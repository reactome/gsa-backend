import unittest
import logging
from datetime import time
import os

from reactome_analysis_datasets.dataset_fetchers.grein_fetcher import GreinFetcher
from reactome_analysis_utils.models.dataset_request import DatasetRequestParameter
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import ExternalData



class MockMQ:
    def get_is_shutdown(self):
        return False
    def sleep(self, the_time):
        time.sleep(the_time)
class GreinFetcherTest(unittest.TestCase):
    test_dataset = "GSE112749"
    test_dataset_request = "GSE100007"
    def setUp(self):
       logging.basicConfig(level=logging.DEBUG)

    def test_get_identifier(self):
        parameters = [DatasetRequestParameter("dataset_id", self.test_dataset)]
        fetcher = GreinFetcher()
        id_param = fetcher.get_dataset_id(parameters)
        self.assertIsNotNone(id_param)
        self.assertEqual("GSE112749", id_param)

    def test_load_dataset(self):
        parameters = [DatasetRequestParameter("dataset_id", self.test_dataset_request)]
        fetcher = GreinFetcher()
        external_data = ExternalData()
        external_data = fetcher.load_dataset(parameters, MockMQ)
        self.assertEqual(len(external_data), 2, "Data missing in return")

    def test_load_dataset_2(self):
        os.environ["USE_GREIN_PROXY"] = "True"
        parameters = [DatasetRequestParameter("dataset_id", "GSE100040")]
        fetcher = GreinFetcher()
        external_data = ExternalData()
        external_data = fetcher.load_dataset(parameters, MockMQ)
        self.assertEqual(len(external_data), 2, "Data missing in return")

    def test_overview(self):
        fetcher = GreinFetcher()
        overview = fetcher.get_available_datasets(10)
        overview_len = len(overview)
        self.assertIsNotNone(overview)
        self.assertEqual(overview_len, 10)
