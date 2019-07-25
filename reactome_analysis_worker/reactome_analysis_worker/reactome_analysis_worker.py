# -*- coding: utf-8 -*-
"""
Supported environmental variables:
  * REACTOME_VERSION: Version number of the current REACTOME release
"""

import json
import logging
import multiprocessing
import os
import queue

import numpy
import prometheus_client
from reactome_analysis_api.input_deserializer import create_analysis_input_object
from reactome_analysis_api.models.analysis_input import AnalysisInput
from reactome_analysis_api.models.analysis_result import AnalysisResult
from reactome_analysis_api.models.analysis_result_mappings import AnalysisResultMappings
from reactome_analysis_api.models.analysis_status import AnalysisStatus
from reactome_analysis_utils.models import report_request
from reactome_analysis_utils.reactome_mq import ReactomeMQ, REPORT_QUEUE
from reactome_analysis_utils.reactome_storage import ReactomeStorage

from reactome_analysis_worker import result_converter
from reactome_analysis_worker import util
from reactome_analysis_worker.analysers import *
from reactome_analysis_worker.geneset_builder import GeneSet, generate_pathway_filename
from reactome_analysis_worker.models.gene_set_mapping import GeneSetMapping

LOGGER = logging.getLogger(__name__)


RUNNING_ANALYSES = prometheus_client.Gauge("reactome_worker_running_analyses",
                                           "Number of analyses currently running.")
COMPLETED_ANALYSES = prometheus_client.Counter("reactome_worker_completed_analyses",
                                               "Number of successfully completed analyses.")


class ReactomeAnalysisWorker:
    """
    The ReactomeAnalysisWorker class contains all functions
    required to process an analysis request.
    """
    def __init__(self):
        """
        Initializes basic member variables
        """
        self._mq = None
        self._storage = None
        self._performed_analyses = 0
        self.debug = bool(os.getenv("REACTOME_WORKER_DEBUG", False))

    def _get_mq(self):
        """
        Return the current connection to the message queue
        :return: A ReactomeMQ object
        """
        if not self._mq:
            try:
                self._mq = ReactomeMQ()
            except Exception as e:
                LOGGER.error("Failed to connect to MQ service: " + str(e))
                raise Exception("Failed to connect to MQ service.", e)

        return self._mq

    def _get_storage(self):
        """
        Returns the current connection to the reactome storage
        :return: A ReactomeStorage object
        """
        if not self._storage:
            try:
                self._storage = ReactomeStorage()
            except Exception as e:
                LOGGER.error("Failed to connect to storage service: " + str(e))
                raise("Failed to connect to storage service", e)

        return self._storage

    def start_analyses(self):
        """
        Connects to the analysis queue and starts listening for new analysis
        requests. This is a blocking function that will only exit on failure.
        """
        mq = self._get_mq()

        mq.process_analysis(self._on_new_analysis)

    def process_single_message(self):
        """
        Processes a single message from the queue and returns the result.
        This function is only intended for testing purposes
        :return: An analysis result
        """
        mq = self._get_mq()
        mq.process_single_message(self._on_new_analysis)

    def _set_status(self, analysis_id: str, status: str, description: str = None, completed: float = None):
        """
        Set the status of the defined analysis
        :param analysis_id: The analysis's id
        :param status: The status to set ("running", "complete", "failed")
        :param description: (Optional) A more verbose description of the status.
        :param completed: (Optional) The relative completion of the task
        """
        r = self._get_storage()

        # create the status object
        status = AnalysisStatus(analysis_id, status=status, description=description, completed=completed)
        r.set_status(analysis_id, json.dumps(status.to_dict()))

    def _acknowledge_message(self, channel, method):
        """
        Acknowledges a message and thereby marks the analysis as complete.
        :param channel: The channel on which the message was received
        :param method: The method object passed to the analysis
        :return:
        """
        # This function is only here to increase the readability of the code.
        channel.basic_ack(delivery_tag=method.delivery_tag)

        # Decrement the number of running analyses
        RUNNING_ANALYSES.dec()

    def _on_new_analysis(self, ch, method, properties, body):
        """
        Callback function that is triggered whenever a new
        message with an analysis request is received.
        :param ch: The channel the message was received on
        :param method: Method details
        :param properties: Message properties
        :param body: The actual message body (= JSON encoded analysis request)
        """
        LOGGER.debug("Received message.")

        # increment the running analyses
        RUNNING_ANALYSES.inc()

        # create the analysis object
        try:
            LOGGER.debug("Decoding JSON string")
            # decode the JSON information
            body_dict = json.loads(body)
            request = create_analysis_input_object(body_dict)
        except Exception as e:
            # This means that the application has a major problem - this should never happen
            LOGGER.critical("Failed to create analysis request object: " + str(e))
            LOGGER.debug("Error details:", exc_info=1)
            # remove the message from the queue
            self._acknowledge_message(ch, method)

            return

        LOGGER.debug("Received analysis request for " + request.method_name + " (" + request.analysis_id + ")")

        # make sure the request contains datasets
        if len(request.datasets) < 1:
            LOGGER.debug("Analysis request {} does not contain any datasets".format(request.analysis_id))
            # update the status as failed
            self._set_status(request.analysis_id, "failed", description="Request did not contain any datasets",
                             completed=1)
            self._acknowledge_message(ch, method)
            return

        # create the analyser to use
        reactome_analyser = ReactomeAnalyser.get_analyser_for_method(request.method_name.lower())

        # make sure the analyser exists
        if reactome_analyser is None:
            self._set_status(request.analysis_id, status="failed",
                             description="Unsupported method '{}' selected".format(request.method_name),
                             completed=1)
            self._acknowledge_message(ch, method)
            return

        # update the status and mark it as received
        self._set_status(request.analysis_id, status="running", description="Request received", completed=0.05)

        # convert the dataset matrices
        if not self._convert_datasets(request):
            self._acknowledge_message(ch, method)
            return

        # get all identifiers
        all_identifiers = self._extract_identifiers(request.datasets)

        # make sure more than one gene was submitted
        if len(all_identifiers) <= 1:
            LOGGER.debug("Analysis request {} contains an insufficient number of genes ({})".format(
                request.analysis_id, str(len(all_identifiers))))
            self._set_status(request.analysis_id, status="failed", description="Analysis requires >1 genes.",
                             completed=1)
            self._acknowledge_message(ch, method)
            return

        # get the identifier mappings
        self._set_status(request.analysis_id, status="running", description="Mapping identifiers...", completed=0.1)

        try:
            identifier_mappings = util.map_identifiers(all_identifiers, return_all=True)
        except util.MappingException as e:
            LOGGER.debug("Identifier mapping failed", exc_info=1)
            self._set_status(request.analysis_id, status="failed",
                             description="Invalid gene/protein identifiers submitted", completed=1)
            self._acknowledge_message(ch, method)
            return
        except Exception as e:
            LOGGER.error("Failed to connect to mapping service: " + str(e))
            LOGGER.debug("Mapping failed", exc_info=1)
            self._set_status(request.analysis_id, status="failed",
                             description="Failed to contact identifier mapping service. Please try again later.", completed=1)
            self._acknowledge_message(ch, method)
            return

        LOGGER.debug("Mapped {} of {} submitted identifiers".format(
            str(len(identifier_mappings)), str(len(all_identifiers))))

        # make sure that identifiers were mapped
        if len(identifier_mappings) < 1:
            # mark the analysis as failed
            self._set_status(request.analysis_id, status="failed",
                             description="Failed to map any submitted identifiers", completed=1)
            self._acknowledge_message(ch, method)
            return

        # make sure the experimental design matches the number of samples
        if not self._validate_experimental_design(request.datasets, request.analysis_id):
            self._acknowledge_message(ch, method)
            return

        # load the matching gene set
        use_interactors = request.parameter_dict.get("use_interactors", "False").lower() == "true"

        # species is always set at human since we use Reactome's mapping feature "to human"
        gene_set = GeneSet.create_from_file(generate_pathway_filename(resource="reactome",
                                                                      species="Homo sapiens",
                                                                      contains_interactors=use_interactors))

        # filter the datasets
        for dataset in request.datasets:
            dataset.df = self._filter_dataset(dataset.df, identifier_mappings, dataset.design,
                                              max_missing_values=float(request.parameter_dict.get("max_missing_values", 0.5)))

        # get the retained identifiers
        identifiers_after_filter = self._extract_identifiers(datasets=request.datasets)

        # perform the mapping for every dataset
        mappings = dict()
        for dataset in request.datasets:
            LOGGER.debug("Mapping identifiers for dataset {}".format(dataset.name))
            mappings[dataset.name] = GeneSetMapping.create_mapping(gene_set=gene_set,
                                                                   identifiers=dataset.df[:][dataset.df.dtype.names[0]].tolist(),
                                                                   identifier_mapping=identifier_mappings)

        # process the analysis
        self._set_status(request.analysis_id, status="running",
                         description="Performing gene set analysis using {}".format(request.method_name), completed=0.2)

        try:
            # move the analysis to a separate process in order to "stay alive" in the eyes of
            # the queuing system - rpy2 causes python to stop
            is_analysis_complete = multiprocessing.Event()
            result_queue = multiprocessing.Queue()
            status_queue = multiprocessing.Queue()
            analysis_process = AnalysisProcess(analyser=reactome_analyser, request=request, gene_set_mappings=mappings,
                                               gene_set=gene_set,
                                               identifier_mappings=identifier_mappings,
                                               on_complete=is_analysis_complete, result_queue=result_queue,
                                               status_queue=status_queue)
            LOGGER.debug("Launching process to perform the analysis...")

            analysis_process.start()

            # fetch the blueprint (in parallel) for the REACTOME result conversion
            reactome_blueprint = None

            try:
                # only get the blueprint if the option is set
                if request.parameter_dict.get("create_reactome_visualization").lower() == "true":
                    LOGGER.debug("Fetching blueprint for Reactome result conversion")
                    reactome_blueprint = result_converter.perform_reactome_gsa(identifiers=identifiers_after_filter,
                                                                               use_interactors=use_interactors)
            except Exception as e:
                LOGGER.warning("Failed to retrieve Reactome blueprint: " + str(e))

            # wait for completion
            while analysis_process.is_alive() and not is_analysis_complete.is_set():
                # test whether the analysis should be interrupted
                if self._get_mq().get_is_shutdown():
                    LOGGER.debug("Shutdown triggered, terminating analysis process")
                    analysis_process.terminate()
                    analysis_process.join(0.1)
                    return

                if status_queue.qsize() > 0:
                    try:
                        # only use the last update
                        while status_queue.qsize() > 0:
                            status_object = status_queue.get(block=True, timeout=0.5)
                        self._set_status(request.analysis_id, status="running", description=status_object.description,
                                         completed=status_object.completed)
                    except Exception:
                        # this can safely be ignored since it is most commonly caused by the fact that the worker
                        # is too busy and fetching of the message timed out
                        pass

                self._get_mq().sleep(1)

            LOGGER.debug("Analysis process completed. Joining process...")

            # for potential cleanup
            analysis_process.join(1)

            # retrieve the result from the queue
            try:
                gsa_results = result_queue.get(block=True, timeout=0.5)
            except queue.Empty:
                gsa_results = None

            LOGGER.debug("Result received from queue.")

            # make sure a result was received
            if not isinstance(gsa_results, list) or len(gsa_results) < 1:
                LOGGER.error("No analysis result retrieved for {}".format(request.analysis_id))

                # test if an exception is returned instead
                if isinstance(gsa_results, Exception):
                    self._set_status(request.analysis_id, status="failed",
                                     description="{} analysis failed: {}".format(request.method_name, str(gsa_results)),
                                     completed=1)
                else:
                    self._set_status(request.analysis_id, status="failed",
                                     description="{} analysis failed.".format(request.method_name),
                                     completed=1)

                self._acknowledge_message(ch, method)
                return

            # create the AnalysisResult object
            analysis_result = AnalysisResult(release=os.getenv("REACTOME_VERSION", "68"),
                                             results=gsa_results,
                                             mappings=self._convert_mapping_result(identifier_mappings))

            # submit the result to Reactome
            if reactome_blueprint:
                analysis_result.reactome_links = list()
                self._set_status(request.analysis_id, status="running", description="Creating REACTOME visualization",
                                 completed=0.9)

                for reactome_type in reactome_analyser.reactome_result_types:
                    LOGGER.debug("Submitting result for " + reactome_type.name)
                    try:
                        reactome_link = result_converter.submit_result_to_reactome(result=analysis_result,
                                                                                   result_type=reactome_type,
                                                                                   reactome_blueprint=reactome_blueprint,
                                                                                   min_p=0.05)
                        analysis_result.reactome_links.append(reactome_link)
                    except Exception as e:
                        # simply ignore this error
                        LOGGER.warning("Failed to submit result to Reactome: " + str(e))

            # save the result
            storage = self._get_storage()
            storage.set_result(analysis_identifier=request.analysis_id, result=json.dumps(analysis_result.to_dict()))
            # update the status
            self._set_status(request.analysis_id, status="complete", description="Analysis done", completed=1)
            self._acknowledge_message(ch, method)

            # send the request to create the report
            if request.parameter_dict.get("create_reports", "False").lower() == "true" or \
               len(request.parameter_dict.get("email", "")) > 3:
                message_mq = ReactomeMQ(queue_name=REPORT_QUEUE)
                report_request_obj = report_request.ReportRequest(analysis_id=request.analysis_id,
                                                                  user_mail=request.parameter_dict.get("email", None))
                message_mq.post_analysis(analysis=report_request_obj.to_json(), method="report")

            # count the complete analysis
            COMPLETED_ANALYSES.inc()
        except Exception as e:
            self._set_status(request.analysis_id, status="failed", description="Failed to analyse dataset: " + str(e),
                             completed=1)
            self._acknowledge_message(ch, method)
            if self.debug:
                raise e

    def shutdown(self):
        """
        Gracefully shutdown all services
        """
        LOGGER.debug("Shutting down gracefully.")
        if self._mq:
            LOGGER.debug("Closing MQ connection.")
            self._mq.close()

    def _convert_datasets(self, request: AnalysisInput) -> bool:
        """
        Converts the expression matrix of all dataset objects to numpy arrays.
        :param request: The sent request
        :return: A boolean indicating whether the conversion worked.
        """
        # convert matrix for every dataset into an array
        for dataset in request.datasets:
            try:
                dataset.df = util.string_to_array(dataset.data)
            # Mark the analysis as failed if the conversion caused an error.
            except util.ConversionException as e:
                LOGGER.error("Failed to convert dataset '{}' from analysis '{}': {}".format(
                    dataset.name, request.analysis_id, str(e)
                ))

                # mark the analysis as failed
                self._set_status(request.analysis_id, status="failed",
                                 description="Failed to convert dataset '{}'".format(dataset.name), completed=1)

                return False

        return True

    def _extract_identifiers(self, datasets: list) -> set:
        """
        Extract all identifiers from a list of datasets.
        :param datasets: A list of datasets
        :return: A set containing the identifiers
        """
        identifiers = set()

        for dataset in datasets:
            # ndim = 0 only occurs if the dataset contains a single row
            if dataset.df.ndim > 0:
                genes = dataset.df[:][dataset.df.dtype.names[0]].tolist()
            else:
                genes = [dataset.df[dataset.df.dtype.names[0]].item()]

            identifiers.update(genes)

        return identifiers

    def _validate_experimental_design(self, datasets: list, analysis_id: str) -> bool:
        """
        Checks whether the structure of the experimental design matches the structure of
        the expression matrices in every dataset.
        :param datasets: The datasets to test
        :param analysis_id: The analysis' id used to set the status in case it failed
        :return: Boolean to indicate whether the design is valid
        """
        for dataset in datasets:
            # ignore requests that do not contain a design
            if not dataset.design:
                continue

            design_samples = dataset.design.samples
            matrix_samples = dataset.df.dtype.names[1:]

            if not len(design_samples) == len(matrix_samples):
                # mark the analysis as failed
                self._set_status(analysis_id, status="failed",
                                 description="Failed to convert dataset '{}'. Different number of samples in the "
                                             "experimental design ({}) and the expression matrix ({})"
                                 .format(dataset.name, str(len(design_samples)), str(len(matrix_samples))),
                                 completed=1)
                return False

        return True

    def _convert_mapping_result(self, identifier_mappings: dict) -> list:
        """
        Converts the identifier mappings from a dict with lists to a list of IdentifierMappingResult objects
        :param identifier_mappings: The identifier mappings as a dict with the identifier as key and all mappings
               as members of the subsequent list.
        :return: A list of AnalysisResultMappingResult objects
        """
        result = list()

        for identifier in identifier_mappings:
            mapping_obj = AnalysisResultMappings(identifier=identifier, mapped_to=identifier_mappings[identifier])
            result.append(mapping_obj)

        return result

    @staticmethod
    def _filter_dataset(df: numpy.array, identifier_mappings: dict, design, max_missing_values: float):
        """
        Filter the passed expression matrix (a numpy array) removing A) all unmapped genes and
        B) all genes found in less than max_missing_values samples.
        :param df: The numpy array to filter
        :param identifier_mappings: Result of the mapping service
        :param design: The design
        :return: The filtered numpy array
        """
        # remove all unmapped genes by keepingstorage. the mapped ones
        rows_to_keep = list()
        for i in range(df.size):
            gene = df[df.dtype.names[0]][i]
            if gene in identifier_mappings:
                rows_to_keep.append(i)

        LOGGER.debug("Keeping {} mapped identifiers of {}".format(str(len(rows_to_keep)), str(df.size)))
        df = df[rows_to_keep]

        # filter the expression values based on in how many samples they were observed
        if design:
            # get values for group 1
            group1_columns = list()
            group2_columns = list()

            for i in range(len(design.analysis_group)):
                if design.analysis_group[i] == design.comparison.group1:
                    group1_columns.append(df.dtype.names[i + 1])
                elif design.analysis_group[i] == design.comparison.group2:
                    group2_columns.append(df.dtype.names[i + 1])

            # get the proportion of samples a gene as expressed in
            min_missing_values = list()

            for gene_row in df:
                n_group1 = sum([gene_expr == 0 for gene_expr in gene_row[group1_columns].tolist()])
                n_group2 = sum([gene_expr == 0 for gene_expr in gene_row[group2_columns].tolist()])
                # save the maximum relative expression
                min_missing_values.append(min([n_group1 / len(group1_columns), n_group2 / len(group2_columns)]))

            above_expression = [item[0] for item in enumerate(min_missing_values) if item[1] < max_missing_values]
            LOGGER.debug("Keeping {} of {} genes with less than {} missing values".format(str(len(above_expression)),
                                                                                          str(df.size),
                                                                                          str(max_missing_values)))
            df = df[above_expression]
        else:
            # if no design is present, use the relative expression across all samples
            missing_values = list()

            for gene_row in df:
                missing_samples = sum([gene_expr == 0 for gene_expr in gene_row])
                # save the maximum relative expression
                missing_values.append(missing_samples / (len(df.dtype.names) - 1))

            above_expression = [item[0] for item in enumerate(missing_values) if item[1] < max_missing_values]
            LOGGER.debug("Keeping {} of {} genes with less than {} missing values".format(str(len(above_expression)),
                                                                                          str(df.size),
                                                                                          str(max_missing_values)))
            df = df[above_expression]

        return df


class AnalysisProcess(multiprocessing.Process):
    def __init__(self, analyser: ReactomeAnalyser, request: AnalysisInput, gene_set_mappings: dict, gene_set: GeneSet,
                 identifier_mappings: dict, on_complete: multiprocessing.Event, result_queue: multiprocessing.Queue,
                 status_queue: multiprocessing.Queue = None):
        """
        Creates a new AnalysisThread object to perform a gene set analysis
        :param analyser: The analyser to use
        :param request: The AnalysisInput describing the request
        :param gene_set_mappings: The result of the gene set mapping for each dataset
        :param gene_set: The gene set used to create the mapping
        :param identifier_mappings: The identifier mappings for each dataset
        :param on_complete: This event will be triggered once the analysis is complete
        :param result_queue: The result queue to use to return the result
        """
        super().__init__()

        self.analyser = analyser
        self.request = request
        self.gene_set_mappings = gene_set_mappings
        self.gene_set = gene_set
        self.identifier_mappings = identifier_mappings
        self.on_complete = on_complete
        self.result_queue = result_queue
        self.status_queue = status_queue

    def update_status(self, message: str, complete: float):
        if self.status_queue:
            self.status_queue.put(AnalysisStatus(id=self.request.analysis_id, status="running",
                                                 description=message, completed=complete))

    def run(self) -> None:
        try:
            # add the status callback
            LOGGER.debug("Setting status callback for analyser")
            self.analyser.set_status_callback(self.update_status)

            # analyse the request
            LOGGER.debug("Starting analysis....")
            gsa_result = self.analyser.analyse_request(request=self.request, gene_set_mappings=self.gene_set_mappings,
                                                       identifier_mappings=self.identifier_mappings, gene_set=self.gene_set)

            LOGGER.debug("Saving result in queue....")
            self.result_queue.put(gsa_result)
        except Exception as e:
            # put the error message in the queue
            self.result_queue.put(e)
        finally:
            LOGGER.debug("Setting on_complete")
            self.on_complete.set()

