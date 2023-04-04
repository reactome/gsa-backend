import unittest
import logging
from reactome_analysis_datasets.dataset_fetchers.grein_fetcher import GreinFetcher
from reactome_analysis_utils.models.dataset_request import DatasetRequestParameter
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcherException
import os
from reactome_analysis_utils.reactome_storage import redis


class GreinFetcherTest(unittest.TestCase):
    test_dataset = "GSE112749"

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
    def test_get_identifier(self):
        parameters = [DatasetRequestParameter("dataset_id", self.test_dataset)]
        fetcher = GreinFetcher()
        id_param = fetcher.get_dataset_id(parameters)
        self.assertIsNotNone(id_param)
        self.assertEqual("GSE112749", id_param)
