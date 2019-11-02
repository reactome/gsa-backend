"""
A DatasetFetcher to retrieve experiments from EBI's
ExpressionAtlas
"""

from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, DatasetFetcherException
import reactome_analysis_worker
import reactome_analysis_worker.util
from reactome_analysis_utils import reactome_mq
from pkg_resources import resource_string
import re
import urllib3
import json
import enum
import multiprocessing
import tempfile
import rpy2.rinterface as ri
import rpy2.robjects as ro
import rpy2.robjects.packages
import logging
import os

LOGGER = logging.getLogger(__name__)

# initialize R
ri.initr()

# disable R messages
#ri.set_writeconsole_warnerror(None)
#ri.set_writeconsole_regular(None)


class ExpressionAtlasTypes(enum.Enum):
    R_RNA_SEQ = 1
    R_MICROARRAY = 2
    PROTEOMICS = 3


class ExpressionAtlasFetcher(DatasetFetcher):
    """
    A DatasetFetcher to retrieve experiments from EBI's
    ExpressionAtlas
    """
    # URL to retrieve JSON information about RNA-seq experiments:
    # https://www.ebi.ac.uk/fg/rnaseq/api/json/getStudy/E-GEUV-1

    # URL to retrieve JSON formatted metadata
    # https://www.ebi.ac.uk/fg/rnaseq/api/json/getSampleAttributesPerRunByStudy/E-GEUV-1

    # ExpressionAtlas Bioconductor package can retrieve other data as well
    # Microarray example: 

    # Proteomics data is only available through a "direct" download
    # Example: E-PROT-[1-10]
    # 

    # RNA-seq examples:
    # E-MTAB-970

    # Direct download links:
    # All available files: https://www.ebi.ac.uk/gxa/json/experiments/E-ATMX-20/resources/DATA

    def load_dataset(self, identifier: str, reactome_mq: reactome_mq.ReactomeMQ) -> (str, ExternalData):
        """
        Load the specified ExpressionAtlas experiment
        :param identifier: The ExpressionAtlas identifier
        :param reactome_mq: The MQ used to process messages.
        :returns: (data, summary)
        """
        # make sure the identifier matches the pattern
        if not identifier[0:2] == "E-":
            raise DatasetFetcherException("{} is not a valid ExpressionAtlas identifier".format(identifier))

        # get the available files for download
        download_files = self.fetch_available_files(identifier)

        # try to retrieve the type of the passed identifier
        dataset_type = self.get_dataset_type(download_files)

        # load the data based on the type
        if dataset_type == ExpressionAtlasTypes.R_MICROARRAY or dataset_type == ExpressionAtlasTypes.R_RNA_SEQ:
            loaded_data = self.load_r_data(download_files, reactome_mq)
        elif dataset_type == ExpressionAtlasTypes.PROTEOMICS:
            loaded_data = self.load_generic_data(download_files, reactome_mq)
            loaded_data["data_type"] = "proteomics_int"

        # convert the metadata into a summary object
        summary = self._create_summary(identifier, loaded_data["data_type"], loaded_data["metadata"])

        return (loaded_data["expression_values"], summary)

    def _create_summary(self, identifier: str, data_type: str, metadata: str) -> str:
        """
        Create a ExternalData object based on
        the passed metadata
        :param identifier: The ExpressionAtlas identifier.
        :param data_type: The data type of the experiment
        :param metadata: The tab-delimited metadata table.
        :return: The JSON encoded summary object
        """
        # initialize the summary object
        summary = {"type": data_type, "id": identifier, 
                   "title": "ExpressionAtlas dataset {}".format(identifier),
                   "description": "External dataset loaded from ExpressionAtlas"}

        # add the sample specific data
        metadata_array = reactome_analysis_worker.util.string_to_array(metadata)

        # get the sample ids
        summary["sample_ids"] = metadata_array[metadata_array.dtype.names[0]].tolist()

        # add every column as one property
        summary["sample_metadata"] = list()

        for property_name in metadata_array.dtype.names[1:]:
            summary["sample_metadata"].append(
                {
                    "name": property_name,
                    "values": metadata_array[property_name].tolist()
                }
            )

        # normalize RNA-seq count data by default
        if data_type == "rnaseq_counts":
            summary["default_parameters"] = [
                {"name": "discrete_norm_function", "value": "TMM"}
            ]

        return ExternalData.from_dict(summary)

    def load_generic_data(self, download_files: list, reactome_mq: reactome_mq.ReactomeMQ):
        """
        Load the generic ExpressionAtlas data using the provided TSV files.
        :param download_files: The list of files available for download
        :param reactome_mq: The MQ used to process messages.
        :return A dict containing the 'metadata', and 'expression_values'
        """
        # get the metadata file
        design_file = None

        for file in download_files:
            if "experiment-design" in file["url"]:
                design_file = file["url"]
                break

        if design_file is None:
            LOGGER.error("Dataset does not contain an experimental design file")
            raise DatasetFetcherException("Failed to retrieve experimental design from ExpressionAtlas.")

        design_data = self._download_atlas_file(design_file).decode()

        # get the expression values file
        expression_file = None

        for file in download_files:
            if "expression values" in file["description"].lower():
                expression_file = file["url"]
                break

        if expression_file is None:
            LOGGER.error("Dataset does not contain an expression file")
            raise DatasetFetcherException("Failed to retrieve expression values from ExpressionAtlas")

        expression_data = self._download_atlas_file(expression_file).decode()

        # process the metadata file
        clean_metadata = self._filter_metadata(design_data)

        # process the expression data
        clean_expression_data = self._filter_expression_data(expression_data)
        
        return {"metadata": clean_metadata, "expression_values": clean_expression_data}

    def _filter_expression_data(self, expression_data: str) -> str:
        """
        Processes an ExpressionAtlas' expression data file and removes all
        comment lines as well as unwanted columns
        :param expression_data: The expression data as a single string
        :return: The filtered string
        """
        # remove all comment lines
        org_lines = [line for line in expression_data.split("\n") if len(line) > 0 and line[0] != "#"]

        # filter unwanted columns (only the Gene Name)
        header_fields = org_lines[0].split("\t")
        unwanted_columns = dict()

        if header_fields[0] == "Gene ID" and header_fields[1] == "Gene Name":
            unwanted_columns[1] = 1

        # remove the unwanted columns
        filtered_lines = list()

        for line in org_lines:
            fields = line.split("\t")
            filtered_fields = list()

            for field_index in range(0, len(fields)):
                if field_index not in unwanted_columns:
                    filtered_fields.append(fields[field_index])

            filtered_lines.append("\t".join(filtered_fields))

        # return as one string
        return "\n".join(filtered_lines)


    def _filter_metadata(self, metadata_string: str) -> str:
        """
        Processes an ExpressionAtlas experimental design file and
        removes all unwanted columns (such as ontology terms)
        :param metadata_string: The content of the metadata file to process
        :return: The adapted metadata file as a string
        """
        # process as lines
        org_lines = metadata_string.split("\n")
        
        # get the irrelevant columns
        columns_to_ignore = dict()
        header_fields = org_lines[0].split("\t")

        for header_index in range(0, len(header_fields)):
            if "Ontology Term" in header_fields[header_index]:
                columns_to_ignore[header_index] = 1
            # ignore the "Analysed" column
            if header_fields[header_index] == "Analysed":
                columns_to_ignore[header_index] = 1
            # ignore all factors
            if "Factor Value" in header_fields[header_index]:
                columns_to_ignore[header_index] = 1

        # add all relevant fields to a new list
        filtered_lines = list()

        for org_line in org_lines:
            fields = org_line.split("\t")
            relevant_fields = list()

            # only retain fields that are not in "columns_to_ignore"
            for field_index in range(0, len(fields)):
                if field_index not in columns_to_ignore:
                    relevant_fields.append(fields[field_index])

            # combine back into a single string
            filtered_lines.append("\t".join(relevant_fields))

        # adapt the header
        filtered_header_fields = filtered_lines[0].split("\t")
        clean_header_fields = list()

        for header_field in filtered_header_fields:
            if "Sample Characteristic" in header_field:
                clean_header_fields.append(header_field[22:len(header_field) - 1])
            else:
                clean_header_fields.append(header_field)

        filtered_lines[0] = "\t".join(clean_header_fields)

        # combine back into a single string
        return "\n".join(filtered_lines)

    def _download_atlas_file(self, file_url: str):
        """
        Download the specified file from ExpressionAtlas.
        :param file_url: The file's relative URL
        :return: The file's content as a binary
        """
        # download the file
        http = urllib3.PoolManager()
        file_url = "https://www.ebi.ac.uk/gxa/" + file_url
        request = http.request("GET", file_url)

        if request.status != 200:
            LOGGER.error("Failed to download file from ExpressionAtlas ({})".format(str(request.status)))
            raise DatasetFetcherException("Failed to download file from ExpressionAtlas")

        return request.data


    def load_r_data(self, download_files: list, reactome_mq: reactome_mq.ReactomeMQ):
        """
        Load the R-based expression data for the specified dataset
        :param download_files: The list of files available for download
        :param reactome_mq: The MQ used to process messages.
        :return A dict containing the 'data_type', 'metadata', and 'expression_values'
        """
        # get the R file
        r_file = None

        for file in download_files:
            if file["type"] == "icon-Rdata":
                r_file = file["url"]
                break

        if not r_file:
            raise DatasetFetcherException("Failed to locate R file for ExpressionAtlas dataset")

        # download the file
        r_file_content = self._download_atlas_file(r_file)

        # save the file
        try:
            stored_r_file = tempfile.NamedTemporaryFile(suffix=".RData", delete=False)
            stored_r_file.write(r_file_content)
        except Exception as e:
            LOGGER.error("Failed to store ExpressionAtlas file: " + str(e))
            raise DatasetFetcherException("Failed to retrieve ExpressionAtlasFile")
        finally:
            stored_r_file.close()

        # Use the RLoadingProcess to load the file
        file_loaded_event = multiprocessing.Event()
        file_loading_queue = multiprocessing.Queue()

        loading_process = RLoadingProcess(r_file_path=stored_r_file.name, on_complete=file_loaded_event, result_queue=file_loading_queue)
        loading_process.start()

        # wait until the process is complete and get the data
        while loading_process.is_alive and not file_loaded_event.is_set():
            # test whether the fetching was stopped
            if reactome_mq.get_is_shutdown():
                LOGGER.debug("Shutdown triggered, terminating fetching process")
                loading_process.terminate()
                loading_process.join(0.1)
                return

            # sleep for 1 sec
            reactome_mq.sleep(1)

        # Return the required matrices
        loading_process.join(0.5)

        # delete the file
        os.remove(stored_r_file.name)

        try:
            loaded_data = file_loading_queue.get(block=True, timeout=0.5)
        except multiprocessing.queues.Empty:
            loaded_data = None

        # make sure a result was received
        if not isinstance(loaded_data, dict) or len(loaded_data) != 3:
            LOGGER.error("Failed to process R data")

            # test if an Exception was returned
            if isinstance(loaded_data, Exception):
                raise DatasetFetcherException(loaded_data)
            else:
                raise DatasetFetcherException("Failed to process R data")

        return loaded_data

    def get_dataset_type(self, available_files: list) -> ExpressionAtlasTypes:
        """
        Try to guess the type of dataset based on the available files
        :param available_files: The list of available files
        :return: The type of ExpressionAtlas dataset or None if the type cannot be determined
        """
        # check if an R File is available
        has_r_file = False
        for file in available_files:
            if file["type"] == "icon-Rdata":
                has_r_file = True
                break

        # check if it's a microarray experiment
        is_microarray = False
        is_rnaseq = False
        is_proteomics = False

        for file in available_files:
            if "Microarray" in file["url"]:
                is_microarray = True
                break

            if "RnaSeq" in file["url"]:
                is_rnaseq = True
                break

            if "E-PROT-" in file["url"]:
                is_proteomics = True
                break

        if is_rnaseq and has_r_file:
            return ExpressionAtlasTypes.R_RNA_SEQ
        if is_microarray and has_r_file:
            return ExpressionAtlasTypes.R_MICROARRAY
        if is_proteomics:
            return ExpressionAtlasTypes.PROTEOMICS

        return None


    def fetch_available_files(self, identifier: str) -> list:
        """
        Retrieves the available files for the passed ExpressionAtlas experiment
        :param identifier: The ExpressionAtlas identifier
        :returns: A list containing dicts with the respective fields of the available files
        """
        # get the JSON string of available files
        json_data = self._download_atlas_file("json/experiments/{}/resources/DATA".format(identifier))
        
        # decode the JSON string
        file_list = json.loads(json_data)

        # create the result object
        available_files = list()

        for file_item in file_list:
            if "type" not in file_item or "url" not in file_item or "description" not in file_item:
                raise DatasetFetcherException("Invalid file list retrieved from ExpressionAtlas for {}".format(identifier))

            available_files.append({"type": file_item["type"], "url": file_item["url"],
                                    "description": file_item["description"]})

        return available_files


class RLoadingProcess(multiprocessing.Process):
    """
    Class representing a process to load
    the R file and return the required fields
    """
    def __init__(self, r_file_path: str, on_complete: multiprocessing.Event, result_queue: multiprocessing.Queue):
        super().__init__()

        self.r_file_path = r_file_path
        self.on_complete = on_complete
        self.result_queue = result_queue

        # load the r_code
        LOGGER.debug("Loading required R preprocessing functions")
        r_code = resource_string("reactome_analysis_worker.resources.r_code", "preprocessing_functions.R").decode()
        self.preprocessing_functions = ro.packages.SignatureTranslatedAnonymousPackage(r_code, "reactome_preprocessing")

    def run(self) -> None:
        try:
            # inject the r_file_path into the R session
            ri.globalenv["filename"] = ri.StrSexpVector([self.r_file_path])

            # process the file using the R session
            LOGGER.debug("Extracting expression data from R file")
            ro.reval("""
                # load the file
                load(filename)

                # make sure it only contains 1 experiment
                if (length(experimentSummary) != 1) {
                    stop("Error: Unexpected number of experiments")
                }

                # get the metadata
                metadata <- NA
                data_type <- NA
                expression_values <- NA

                # test if it's RNA-seq or microarray
                if (is(experimentSummary[[1]], "RangedSummarizedExperiment")) {
                    # load the required library
                    library(SummarizedExperiment)

                    data_type <- "rnaseq_counts"
                    metadata <- colData(experimentSummary[[1]])

                    expression_values <- assays(experimentSummary[[1]])$counts
                } else if (is(experimentSummary[[1]], "ExpressionSet")) {
                    # load the required library
                    library(Biobase)

                    # TODO: test for microarray normalisation
                    data_type <- "microarray_norm"
                    metadata <- pData(experimentSummary[[1]])

                    expression_values <- data.frame(exprs(experimentSummary[[1]]))
                } else {
                    stop("Error: Unknown assay type encountered.")
                }

                # add a "sample.id" column for the metadata
                metadata$sample.id <- rownames(metadata)
                metadata <- metadata[, c(ncol(metadata), 1:(ncol(metadata)-1))]
            """)

            # convert the R objects to python objects
            LOGGER.debug("Converting R objects to python")
            data_type = str(ri.globalenv["data_type"][0])
            metadata_string = str(self.preprocessing_functions.data_frame_as_string(ri.globalenv["metadata"])[0])
            expression_value_string = str(self.preprocessing_functions.data_frame_as_string(ri.globalenv["expression_values"])[0])

            # save the result and mark the on_complete event
            LOGGER.debug("Returning results through the queue")
            self.result_queue.put({'data_type': data_type, 'metadata': metadata_string, 'expression_values': expression_value_string})
            self.on_complete.set()

        except Exception as e:
            # put the error message in the queue
            LOGGER.error("Error during loading of R file: " + str(e))
            self.result_queue.put(e)
        finally:
            LOGGER.debug("Setting on_complete")
            self.on_complete.set()