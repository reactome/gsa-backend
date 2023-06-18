import logging
import unittest

from reactome_analysis_datasets.dataset_fetchers.public_data_searcher import Generate_search_values, Searcher


class PublicDatatsetSearch(unittest.TestCase):
    keyword ="TFIIH kinase"

    def setUp(self):
       logging.basicConfig(level=logging.DEBUG)

    def test_search(self):
        Generate_search_values.setup_search_events()
        search = Searcher()
        result_dict = search.index_search(['description', 'title'], self.keyword)
        self.assertIsNone(result_dict)
