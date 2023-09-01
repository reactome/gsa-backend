import unittest
import logging
from datetime import time

from reactome_analysis_datasets.dataset_fetchers.geo_fetcher import GeoFetcher
from reactome_analysis_utils.models.dataset_request import DatasetRequestParameter
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import ExternalData


class MockMQ:
    def get_is_shutdown(self):
        return False
    def sleep(self, the_time):
        time.sleep(the_time)


class GeoFetcherTest(unittest.TestCase):
    test_dataset = "GSE1563"

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)

    def test_load_dataset(self):
        parameters = [DatasetRequestParameter("dataset_id", self.test_dataset)]
        fetcher = GeoFetcher()
        external_data = ExternalData()
        external_data = fetcher.load_dataset(parameters, MockMQ)
        self.assertEqual(len(external_data), 2, "Data missing in return")
    
