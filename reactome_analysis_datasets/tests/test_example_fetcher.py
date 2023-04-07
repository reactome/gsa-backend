import unittest
import os
import logging
import json
from reactome_analysis_datasets import reactome_analysis_dataset_fetcher
from reactome_analysis_utils import reactome_mq, models, reactome_storage
from reactome_analysis_utils.reactome_storage import redis


class DatasetFetcherTest(unittest.TestCase):
    def setUp(self):
        os.environ["REDIS_HOST"] = "192.168.99.100"
        os.environ["REDIS_PORT"] = "32725"
        os.environ["REDIS_PASSWORD"] = "test"
        os.environ["RABBIT_HOST"] = "192.168.99.100"
        os.environ["RABBIT_PORT"] = "30186"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

        logging.basicConfig(level=logging.DEBUG)
        pika_logger = logging.getLogger("pika")
        pika_logger.setLevel(logging.ERROR)

        util_logger = logging.getLogger("reactome_analysis_datasets")
        util_logger.setLevel(logging.DEBUG)

        # delete the datasets from redis
        this_redis = redis.Redis(host=os.getenv("REDIS_HOST"), port=int(os.getenv("REDIS_PORT")),
                                 password=os.getenv("REDIS_PASSWORD"))

        this_redis.delete("request_data:EXAMPLE_1")
        this_redis.delete("request_data:EXAMPLE_MEL_PROT")

    def process_all_messages(self):
        fetcher = reactome_analysis_dataset_fetcher.ReactomeAnalysisDatasetFetcher()
        fetcher.start_listening()

    def test_fetch_dataset(self):
        # set the data directory
        os.environ["EXAMPLE_DIRECTORY"] = os.path.join(os.path.dirname(__file__), "testfiles")

        fetcher = reactome_analysis_dataset_fetcher.ReactomeAnalysisDatasetFetcher()

        # post one analysis request
        mq = reactome_mq.ReactomeMQ(queue_name=reactome_mq.DATASET_QUEUE)
        storage = reactome_storage.ReactomeStorage()

        mq.post_analysis(models.dataset_request.DatasetRequest(
            loading_id="loading_1",
            resource_id="example_datasets",
            parameters=[models.dataset_request.DatasetRequestParameter(name="dataset_id", value="EXAMPLE_1")])
                         .to_json(), method="test")

        # process the message
        fetcher.process_single_message()

        # make sure the status is complete
        status = storage.get_status(analysis_identifier="loading_1", data_type="dataset")

        self.assertIsNotNone(status)
        status_obj = json.loads(status)
        self.assertEqual("complete", status_obj["status"])

        # get the summary
        summary_string = storage.get_request_data_summary("EXAMPLE_1")
        self.assertIsNotNone(summary_string)

        # get the data
        data_string = storage.get_request_data("EXAMPLE_1")
        self.assertIsNotNone(data_string)

    def test_fetch_proteomics(self):
        os.environ["EXAMPLE_DIRECTORY"] = os.path.join(os.path.dirname(__file__), "../example_datasets")
        fetcher = reactome_analysis_dataset_fetcher.ReactomeAnalysisDatasetFetcher()

        # post one analysis request
        mq = reactome_mq.ReactomeMQ(queue_name=reactome_mq.DATASET_QUEUE)
        storage = reactome_storage.ReactomeStorage()

        mq.post_analysis(models.dataset_request.DatasetRequest(
            loading_id="loading_2",
            resource_id="example_datasets",
            parameters=[
                models.dataset_request.DatasetRequestParameter(name="dataset_id", value="EXAMPLE_MEL_PROT")]).to_json(),
                         method="test")

        # process the message
        fetcher.process_single_message()

        # make sure the status is complete
        status = storage.get_status(analysis_identifier="loading_2", data_type="dataset")

        self.assertIsNotNone(status)
        status_obj = json.loads(status)
        self.assertEqual("complete", status_obj["status"])

        # get the summary
        summary_string = storage.get_request_data_summary("EXAMPLE_MEL_PROT")
        self.assertIsNotNone(summary_string)

        # get the data
        data_string = storage.get_request_data("EXAMPLE_MEL_PROT")
        self.assertIsNotNone(data_string)


if __name__ == '__main__':
    unittest.main()
