import unittest
import time
import logging
import json
import os
from reactome_analysis_datasets.dataset_fetchers.expression_atlas_fetcher import ExpressionAtlasFetcher, ExpressionAtlasTypes
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcherException
from reactome_analysis_utils.models.dataset_request import DatasetRequestParameter


class MockMQ:
    def get_is_shutdown(self):
        return False
    def sleep(self, the_time):
        time.sleep(the_time)

class ExpressionAtlasFetcherTest(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)

    def test_available_file_fetching(self):
        fetcher = ExpressionAtlasFetcher()

        files_1 = fetcher.fetch_available_files("E-MTAB-970")
        self.assertEqual(7, len(files_1))

        files_2 = fetcher.fetch_available_files("E-PROT-5")
        self.assertEqual(2, len(files_2))

    def test_dataset_type(self):
        fetcher = ExpressionAtlasFetcher()

        files_1 = fetcher.fetch_available_files("E-MTAB-970")
        type_1 = fetcher.get_dataset_type(files_1)
        self.assertEqual(ExpressionAtlasTypes.R_RNA_SEQ, type_1)

        files_2 = fetcher.fetch_available_files("E-PROT-5")
        type_2 = fetcher.get_dataset_type(files_2)
        self.assertEqual(ExpressionAtlasTypes.PROTEOMICS, type_2)

        # E-GEOD-13316
        files_3 = fetcher.fetch_available_files("E-GEOD-13316")
        type_3 = fetcher.get_dataset_type(files_3)
        self.assertEqual(ExpressionAtlasTypes.R_MICROARRAY, type_3)

    def test_load_r_file(self):
        fetcher = ExpressionAtlasFetcher()

        # Test RNA-seq experiments
        files_1 = fetcher.fetch_available_files("E-MTAB-970")
        loaded_data = fetcher.load_r_data(files_1, MockMQ())

        self.assertIsNotNone(loaded_data)
        self.assertEqual(3, len(loaded_data))
        self.assertTrue("data_type" in loaded_data)
        self.assertTrue("metadata" in loaded_data)
        self.assertTrue("expression_values" in loaded_data)

        self.assertEqual("rnaseq_counts", loaded_data["data_type"])

    def test_load_r_microarray(self):
        fetcher = ExpressionAtlasFetcher()
        
        # Test microarray
        files_3 = fetcher.fetch_available_files("E-GEOD-13316")
        loaded_data = fetcher.load_r_data(files_3, MockMQ())

        self.assertIsNotNone(loaded_data)
        self.assertEqual(3, len(loaded_data))
        self.assertTrue("data_type" in loaded_data)
        self.assertTrue("metadata" in loaded_data)
        self.assertTrue("expression_values" in loaded_data)

        self.assertEqual("microarray_norm", loaded_data["data_type"])

    def test_load_generic_data(self):
        fetcher = ExpressionAtlasFetcher()
        files_2 = fetcher.fetch_available_files("E-PROT-5")

        loaded_data = fetcher.load_generic_data(files_2, None)

        self.assertEqual(2, len(loaded_data))
        self.assertTrue("metadata" in loaded_data)
        self.assertTrue("expression_values" in loaded_data)

        # test the metadata
        self.assertEqual(20, len(loaded_data["metadata"].split("\n")))
        self.assertEqual(10934, len(loaded_data["expression_values"].split("\n")))

    def test_complete_loading_process(self):
        fetcher = ExpressionAtlasFetcher()
        (data, summary) = fetcher.load_dataset([DatasetRequestParameter("dataset_id", "E-PROT-5")], MockMQ())

        self.assertIsNotNone(data)
        self.assertIsNotNone(summary)

        self.assertEqual("proteomics_int", summary.type)

    def test_timeout(self):
        # set the timeout to only 1 sec - to trigger a failure
        os.environ["LOADING_MAX_TIMEOUT"] = "1"

        fetcher = ExpressionAtlasFetcher()
        self.assertRaises(DatasetFetcherException, fetcher.load_dataset, [DatasetRequestParameter("dataset_id", "E-MTAB-970")], MockMQ())

    def test_get_identifier(self):
        parameters = [DatasetRequestParameter("dataset_id", "E-MTAB-970")]

        fetcher = ExpressionAtlasFetcher()

        id_param = fetcher.get_dataset_id(parameters)

        self.assertIsNotNone(id_param)
        self.assertEqual("E-MTAB-970", id_param)

    def test_benchmark_loading(self):
        os.environ["LOADING_MAX_TIMEOUT"] = "1000"
        parameters = [DatasetRequestParameter("dataset_id", "E-ENAD-33")]
        fetcher = ExpressionAtlasFetcher()

        fetcher.load_dataset(parameters, MockMQ())
