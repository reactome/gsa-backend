import unittest
import time
import logging
import json
import os
from reactome_analysis_datasets.dataset_fetchers.sc_expression_atlas_fetcher import ScExpressionAtlasFetcher
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcherException
from reactome_analysis_utils.models.dataset_request import DatasetRequestParameter


class MockMQ:
    def get_is_shutdown(self):
        return False
    def sleep(self, the_time):
        time.sleep(the_time)

class ScExpressionAtlasFetcherTest(unittest.TestCase):
    def setUp(self):
        self.default_parameters = [
            DatasetRequestParameter("dataset_id", "E-CURD-11"),
            DatasetRequestParameter("k", "3")
        ]

        logging.basicConfig(level=logging.DEBUG)

    def test_get_dataset_id(self):
        fetcher = ScExpressionAtlasFetcher()

        dataset_id = fetcher.get_dataset_id(self.default_parameters)

        self.assertIsNotNone(dataset_id)
        self.assertEqual("E-CURD-11_3", dataset_id)

    def test_get_cell_clustering(self):
        fetcher = ScExpressionAtlasFetcher()

        cell_clustering = fetcher._get_cell_clusterings(dataset_id="E-CURD-11", k=3)

        self.assertIsNotNone(cell_clustering)
        self.assertEqual(3, len(cell_clustering))
        self.assertTrue("SRR2049340" in cell_clustering["Cluster 2"])
        self.assertTrue("SRR2049535" in cell_clustering["Cluster 3"])

    def test_download_zip_file(self):
        fetcher = ScExpressionAtlasFetcher()

        file_url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-CURD-11/download/zip?fileType=normalised&accessKey="
        tmp_dir = fetcher._download_zip_file(file_url=file_url)

        self.assertIsNotNone(tmp_dir)
        self.assertTrue(os.path.isdir(tmp_dir.name))
        self.assertTrue(os.path.isfile(os.path.join(tmp_dir.name, "E-CURD-11.aggregated_filtered_normalised_counts.mtx")))

    def test_get_av_cluster_expression(self):
        file_dir = os.path.join(os.path.dirname(__file__), "testfiles")

        class TmpClass:
            pass

        file_dir_obj = TmpClass()
        file_dir_obj.name = file_dir

        fetcher = ScExpressionAtlasFetcher()

        (exp, rows, cols) = fetcher._get_av_cluster_expression(matrix_file_dir=file_dir_obj, 
                                    cell_clusterings={"c1": ["Sample 1"], "c2": ["Sample 2"], "c3": ["Sample 3"]})

        self.assertIsNotNone(exp)
        self.assertIsNotNone(rows)
        self.assertIsNotNone(cols)

        # also test the creation of the final table
        exp_table = fetcher._create_expression_table(exp, rows, cols)

        self.assertEqual("\tc1\tc2\tc3\nGene 1\t1.0\t5.0\t9.0\nGene 2\t2.0\t6.0\t10.0\nGene 3\t3.0\t7.0\t1.1\nGene 4\t4.0\t8.0\t1.2", exp_table)

    def test_extract_factor_values(self):
        test_file = os.path.join(os.path.dirname(__file__), "testfiles", "E-CURD-11.sdrf.txt")

        fetcher = ScExpressionAtlasFetcher()

        cell_factors = fetcher._extract_factor_values_from_sdrf(test_file)

        self.assertEqual(176, len(cell_factors))

    def test_create_summary(self):
        fetcher = ScExpressionAtlasFetcher()

        # get the data
        file_url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-CURD-11/download/zip?fileType=normalised&accessKey="
        tmp_dir = fetcher._download_zip_file(file_url=file_url)

        # get the cell clustering
        cell_clustering = fetcher._get_cell_clusterings(dataset_id="E-CURD-11", k=3)

        (exp, rows, cols) = fetcher._get_av_cluster_expression(matrix_file_dir=tmp_dir, 
                                                               cell_clusterings=cell_clustering)

        # create the summary
        summary = fetcher._create_summary(dataset_id="E-CURD-11", k=3, sample_ids=cols, cell_clusterings=cell_clustering)

        self.assertIsNotNone(summary)
        self.assertEqual(1, len(summary.default_parameters))

    def test_load_exp_design(self):
        dataset_id = "E-HCAD-13"

        fetcher = ScExpressionAtlasFetcher()

        exp_design = fetcher._load_experiment_design_factors(dataset_id)

        self.assertIsNotNone(exp_design)
        self.assertEqual(6263, len(exp_design))
        self.assertTrue("age" in exp_design[list(exp_design.keys())[0]])

    def test_failed_loading(self):
        fetcher = ScExpressionAtlasFetcher()

        failed_experiment = [
            DatasetRequestParameter(name="dataset_id", value="E-HCAD-13"),
            DatasetRequestParameter(name="k", value="12")]

        fetcher.load_dataset(failed_experiment, MockMQ())
        