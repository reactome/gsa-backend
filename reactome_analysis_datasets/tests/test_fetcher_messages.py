import unittest
import os
import logging
import json
from reactome_analysis_datasets import reactome_analysis_dataset_fetcher
from reactome_analysis_utils import reactome_mq, models, reactome_storage
from reactome_analysis_utils.reactome_storage import redis

"""
This testcase tests the complete dataset fetching process
including passing message through the message queue.
Therefore, it needs access to a redis and rabbit mq instance.
"""

class FetcherMessageTest(unittest.TestCase):
    def setUp(self):
        os.environ["REDIS_HOST"] = "localhost"
        os.environ["REDIS_PORT"] = "32088"
        os.environ["REDIS_PASSWORD"] = "test"
        os.environ["RABBIT_HOST"] = "localhost"
        os.environ["RABBIT_PORT"] = "31809"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

        logging.basicConfig(level=logging.DEBUG)
        pika_logger = logging.getLogger("pika")
        pika_logger.setLevel(logging.ERROR)

        util_logger = logging.getLogger("reactome_analysis_datasets")
        util_logger.setLevel(logging.DEBUG)

        # delete the datasets from redis
        this_redis = redis.Redis(host=os.getenv("REDIS_HOST"), port=int(os.getenv("REDIS_PORT")), password=os.getenv("REDIS_PASSWORD"))

        this_redis.delete("request_data:EXAMPLE_1")
        this_redis.delete("request_data:EXAMPLE_MEL_PROT")
        this_redis.delete("request_data:E-GEOD-13316")
        this_redis.delete("request_data:E-MTAB-7078_3")

    def test_process_message(self):
        fetcher = reactome_analysis_dataset_fetcher.ReactomeAnalysisDatasetFetcher()
        fetcher.process_single_message()

    def test_fetch_expression_atlas(self):
        dataset_id = "E-GEOD-13316"
        resource_id = "ebi_gxa"
        loading_id = "loading_gxa_1"

        # set the data directory
        os.environ["EXAMPLE_DIRECTORY"] = os.path.join(os.path.dirname(__file__), "testfiles")

        fetcher = reactome_analysis_dataset_fetcher.ReactomeAnalysisDatasetFetcher()

        # post one analysis request
        mq = reactome_mq.ReactomeMQ(queue_name=reactome_mq.DATASET_QUEUE)
        storage = reactome_storage.ReactomeStorage()

        mq.post_analysis(models.dataset_request.DatasetRequest(
            loading_id=loading_id, 
            resource_id=resource_id, 
            parameters=[models.dataset_request.DatasetRequestParameter(name="dataset_id", value=dataset_id)])
            .to_json(), method="test")

        # process the message
        fetcher.process_single_message()

        # make sure the status is complete
        status = storage.get_status(analysis_identifier=loading_id, data_type="dataset")

        self.assertIsNotNone(status)
        status_obj = json.loads(status)
        self.assertEqual("complete", status_obj["status"])

        # get the summary
        summary_string = storage.get_request_data_summary(dataset_id)
        self.assertIsNotNone(summary_string)

        # get the data
        data_string = storage.get_request_data(dataset_id)
        self.assertIsNotNone(data_string)

    def test_fetch_sc_expression_atlas(self):
        dataset_id = "E-MTAB-7078"
        k = "3"
        resource_id = "ebi_sc_gxa"
        loading_id = "loading_sc_gxa_1"

        # set the data directory
        os.environ["EXAMPLE_DIRECTORY"] = os.path.join(os.path.dirname(__file__), "testfiles")

        fetcher = reactome_analysis_dataset_fetcher.ReactomeAnalysisDatasetFetcher()

        # post one analysis request
        mq = reactome_mq.ReactomeMQ(queue_name=reactome_mq.DATASET_QUEUE)
        storage = reactome_storage.ReactomeStorage()

        mq.post_analysis(models.dataset_request.DatasetRequest(
            loading_id=loading_id, 
            resource_id=resource_id, 
            parameters=[models.dataset_request.DatasetRequestParameter(name="dataset_id", value=dataset_id), 
                        models.dataset_request.DatasetRequestParameter(name="k", value=k)])
            .to_json(), method="test")

        # process the message
        fetcher.process_single_message()

        # make sure the status is complete
        status = storage.get_status(analysis_identifier=loading_id, data_type="dataset")

        self.assertIsNotNone(status)
        status_obj = json.loads(status)
        self.assertEqual("complete", status_obj["status"])
        dataset_id = status_obj["dataset_id"]

        # get the summary
        summary_string = storage.get_request_data_summary(dataset_id)
        self.assertIsNotNone(summary_string)

        # get the data
        data_string = storage.get_request_data(dataset_id)
        self.assertIsNotNone(data_string)
