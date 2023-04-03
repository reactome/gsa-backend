import unittest
import logging
from reactome_analysis_datasets.dataset_fetchers.grein_fetcher import GreinFetcher
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcherException


class GreinFetcherTest(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)

    def test_overview(self):
        fetcher = GreinFetcher()
        overview_test = fetcher.load_overview(2)
