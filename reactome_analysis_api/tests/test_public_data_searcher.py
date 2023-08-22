import json
import os
import unittest
import logging
import sys

sys.path.insert(0, os.path.dirname(__file__))

from reactome_analysis_api.searcher.public_data_searcher import PublicDatasetSearcher


class PublicDataSearcherTest(unittest.TestCase):
    def test_create_search_index(self):
        # set the logging
        logging.basicConfig(level=logging.DEBUG)

        # use a default path
        path = "/tmp/search_index"

        # create the search
        searcher = PublicDatasetSearcher(path=path)
        searcher.setup_search_events()

        # test get species
        species = searcher.get_species()

        self.assertTrue(len(species) > 20)

        # test search
        search_result = searcher.index_search(keyword="melanoma")

        self.assertIsNotNone(search_result)
