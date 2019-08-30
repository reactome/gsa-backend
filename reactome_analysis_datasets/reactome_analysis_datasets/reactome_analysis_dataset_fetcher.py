"""
Class that is listening to dataset loading requests and
loads the actual datasets
"""

import json
import logging

from reactome_analysis_utils import reactome_mq, reactome_storage
from reactome_analysis_utils.models import dataset_request
from reactome_analysis_api.models.analysis_status import AnalysisStatus
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher
from reactome_analysis_datasets.dataset_fetchers.example_fetcher import ExampleDatasetFetcher


LOGGER = logging.getLogger(__name__)


class ReactomeAnalysisDatasetFetcher:
    def __init__(self):
        """
        Initializes basic member variables
        """
        self._mq = None
        self._storage = None
        self._performed_analyses = 0

    def _get_mq(self):
        """
        Return the current connection to the message queue
        :return: A ReactomeMQ object
        """
        if not self._mq:
            try:
                self._mq = reactome_mq.ReactomeMQ(queue_name=reactome_mq.DATASET_QUEUE)
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
                self._storage = reactome_storage.ReactomeStorage()
            except Exception as e:
                LOGGER.error("Failed to connect to storage service: " + str(e))
                raise Exception("Failed to connect to storage service", e)

        return self._storage

    def shutdown(self):
        """
        Gracefully shutdown all services
        """
        LOGGER.debug("Shutting down gracefully.")
        if self._mq:
            LOGGER.debug("Closing MQ connection.")
            self._mq.close()

    def start_listening(self):
        """
        Connects to the report queue and starts listening for new report creation
        requests. This is a blocking function that will only exit on failure.
        """
        mq = self._get_mq()

        LOGGER.debug("Listening for messages...")
        mq.process_analysis(self._on_new_request)

    def process_single_message(self):
        """
        Processes a single message from the queue and returns the result.
        This function is only intended for testing purposes
        :return: An analysis result
        """
        mq = self._get_mq()
        mq.process_single_message(self._on_new_request)
    
    def _acknowledge_message(self, channel, method):
        """
        Acknowledges a message and thereby marks the report as complete.
        :param channel: The channel on which the message was received
        :param method: The method object passed to the analysis
        :return:
        """
        # This function is only here to increase the readability of the code.
        channel.basic_ack(delivery_tag=method.delivery_tag)

    def _set_status(self, request_id: str, status: str, completion: float = 0,
                    description: str = None, reports: list = None):
        """
        Set the status for the dataset loading request
        :param request_id: The request's id
        :param status: Current status ("running", "complete", or "failed")
        :param completion: Relative amount of completion
        :param description: Text describing the current status
        """
        self._get_storage().set_status(
            analysis_identifier=request_id, data_type="dataset",
            status=json.dumps(AnalysisStatus(
                id=request_id, status=status, completed=completion, description=description)
                              .to_dict()))

    def _on_new_request(self, ch, method, properties, body):
        """
        Triggered by a new request for the report generation
        :param ch: The channel the message was received on
        :param method: Method details
        :param properties: Message properties
        :param body: The actual message body (= dataset id)
        """
        LOGGER.debug("Received message.")

        # create the request object
        try:
            LOGGER.debug("Decoding JSON string")
            # decode the JSON information
            request = dataset_request.from_json(body)
        except Exception as e:
            # This means that the application has a major problem - this should never happen
            LOGGER.critical("Failed to create report request object: " + str(e))
            LOGGER.debug("Error details:", exc_info=1)
            # remove the message from the queue
            self._acknowledge_message(ch, method)

            return

        # connect to the storage system
        try:
            storage = self._get_storage()
        except Exception:
            LOGGER.critical("Failed to connect to storage system")
            self._set_status(request_id=request.loading_id, status="failed",
                             description="Failed to connect to storage system")
            self._acknowledge_message(ch, method)
            return

        # test if the dataset already exists
        if storage.request_token_exists(request.dataset_id) and storage.request_data_summary_exists(request.dataset_id):
            self._set_status(request_id=request.loading_id, status="complete", completion=1,
                            description="Dataset {} available.".format(request.dataset_id))
            self._acknowledge_message(ch, method)
            return

        # update the status that it's being loaded
        self._set_status(request_id=request.loading_id, status="running", completion=0.1,
                         description="Dataset {} is being loaded".format(request.dataset_id))

        # get the loading class for the identifier type
        dataset_fetcher = self._get_dataset_fetcher_for_identifier(request.dataset_id)

        if not dataset_request:
            self._set_status(request_id=request.loading_id, status="failed",
                            description="Failed to resolve identifier '{}'.".format(request.dataset_id))
            self._acknowledge_message(ch, method)
            return

        # try to load the dataset
        try:
            (data, summary) = dataset_fetcher.load_dataset(request.dataset_id)

            if data is None:
                raise Exception("Failed to retrieve data.")
            if summary is None:
                raise Exception("Failed to retrieve dataset summary.")

            # save the data
            storage.set_request_data(token=request.dataset_id, data=data, expire=60*60*6)

            # save the summary
            storage.set_request_data_summary(token=request.dataset_id, data=json.dumps(summary.to_dict()))

            # update the status
            self._set_status(request_id=request.loading_id, status="complete", completion=1,
                            description="Dataset {} available.".format(request.dataset_id))

            # acknowledge the message
            self._acknowledge_message(ch, method)
        except Exception as e:
            self._set_status(request_id=request.loading_id, status="failed",
                            description="Failed to load dataset: {}".format(str(e)))
            self._acknowledge_message(ch, method)

    def _get_dataset_fetcher_for_identifier(self, identifier: str) -> DatasetFetcher:
        """
        Returns the matching DatasetFetcher for the passed identifier type or None
        in case the identifier does not match any known format.
        :param identifier: The identifier to get the DatasetFetcher for.
        :return: The matching DatasetFetcher or None if it does not match any known format.
        """
        if len(identifier) > 8 and identifier[0:8] == "EXAMPLE_":
            return ExampleDatasetFetcher()

        return None