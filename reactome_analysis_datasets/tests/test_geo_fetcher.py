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

        (count_matrix, metadata_obj) = fetcher.load_dataset(parameters, MockMQ)
        
        self.assertIsNotNone(count_matrix)
        self.assertIsNotNone(metadata_obj)

        self.assertEqual("GSE1563", metadata_obj.id)
        self.assertEqual("microarray_norm", metadata_obj.type)

        # make sure all metadata fields have the correct length
        n_samples = len(metadata_obj.sample_ids)

        self.assertEqual(62, n_samples, "Incorrect number of samples")

        for sample_metadata in metadata_obj.sample_metadata:
            self.assertEqual(len(sample_metadata["values"]), n_samples, 
                             f"Incorrect values for {sample_metadata['name']}. Expected {n_samples} but got {len(sample_metadata['values'])}")


    