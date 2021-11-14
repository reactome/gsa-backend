# coding: utf-8

import io
import json
import os
import unittest

import sys
from reactome_analysis_utils.reactome_storage import ReactomeStorage

sys.path.insert(0, os.path.dirname(__file__))

from reactome_analysis_api.__main__ import app


class TestAnalysisController(unittest.TestCase):
    def setUp(self):
        os.environ["REDIS_HOST"] = "localhost"
        os.environ["REDIS_PORT"] = "32606"
        os.environ["REDIS_PASSWORD"] = "test"
        os.environ["RABBIT_HOST"] = "localhost"
        os.environ["RABBIT_PORT"] = "32214"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

        app.app.testing = True

        self.test_tsv = """Sample 1\tSample 2\tSample 3
CD19\t1\t2\t3
CD20\t1\t2\t3
MITF\t1\t2\t3"""

        self.test_csv = """\"Sample 1\";\"Sample 2\";\"Sample 3\"
\"CD19\";1;2;3
\"CD20\";1;2;3
\"MITF\";1;2;3"""

        self.test_csv_gene_header = """\"Gene\";\"Sample 1\";\"Sample 2\";\"Sample 3\"
\"CD19\";1;2;3
\"CD20\";1;2;3
\"MITF\";1;2;3"""

        self.test_incorrect_header = """Gene;\"Sample 1\";\"Sample 2\";\"Sample 3\";Sample 4
\"CD19\";1;2;3
\"CD20\";1;2;3
\"MITF\";1;2;3"""

    def test_failed_get(self):
        with app.app.test_client() as client:
            response = client.get("/upload")

            self.assertEqual(405, response.status_code)

    def test_no_file(self):
        with app.app.test_client() as client:
            response = client.post("/upload")

            self.assertEqual(400, response.status_code)
            error_obj = json.loads(response.data.decode())
            self.assertEqual("Incorrect number of uploaded files. Function requires exactly one file.",
                             error_obj["detail"])

    def test_stored_tsv(self):
        with app.app.test_client() as client:
            response = client.post("/upload", data={"file": (io.BytesIO(self.test_tsv.encode("UTF-8")), "test.tsv")})

            self.assertEqual(200, response.status_code)

            result_obj = json.loads(response.data.decode())

            # make sure the samples are correct
            self.assertEqual("Sample 1:Sample 2:Sample 3", ":".join(result_obj["sample_names"]))
            self.assertEqual(4, result_obj["n_lines"])
            self.assertEqual("CD19:CD20:MITF", ":".join(result_obj["top_identifiers"]))

            # make sure the data was stored correctly
            self.assertIsNotNone(result_obj["data_token"])
            token = result_obj["data_token"]

            # create a new redis instance
            storage = ReactomeStorage()

            self.assertTrue(storage.request_token_exists(token))

            stored_obj = storage.get_request_data(token)

            self.assertEqual("\t" + self.test_tsv, stored_obj.decode("UTF-8"))

    def test_retrieve_tsv(self):
        with app.app.test_client() as client:
            response = client.post("/upload?store=false",
                                   data={"file": (io.BytesIO(self.test_tsv.encode("UTF-8")), "test.tsv")})

            self.assertEqual(200, response.status_code)

            result_obj = json.loads(response.data.decode())

            # make sure the samples are correct
            self.assertEqual("Sample 1:Sample 2:Sample 3", ":".join(result_obj["sample_names"]))
            self.assertEqual(4, result_obj["n_lines"])
            self.assertEqual("CD19:CD20:MITF", ":".join(result_obj["top_identifiers"]))

            # make sure the data was stored correctly
            self.assertFalse("data_token" in result_obj)
            self.assertTrue("data" in result_obj)

            self.assertEqual("\t" + self.test_tsv, result_obj["data"])

    def test_retrieve_csv(self):
        with app.app.test_client() as client:
            response = client.post("/upload?store=false",
                                   data={"file": (io.BytesIO(self.test_csv.encode("UTF-8")), "test.csv")})

            self.assertEqual(200, response.status_code)

            result_obj = json.loads(response.data.decode())

            # make sure the samples are correct
            self.assertEqual("Sample 1:Sample 2:Sample 3", ":".join(result_obj["sample_names"]))
            self.assertEqual(4, result_obj["n_lines"])
            self.assertEqual("CD19:CD20:MITF", ":".join(result_obj["top_identifiers"]))

            # make sure the data was stored correctly
            self.assertFalse("data_token" in result_obj)
            self.assertTrue("data" in result_obj)

            # CSV is converted to tsv in the application and quotes removed
            self.assertEqual("\t" + self.test_tsv, result_obj["data"])

    def test_retrieve_csv_gene_header(self):
        with app.app.test_client() as client:
            response = client.post("/upload?store=false",
                                   data={"file": (io.BytesIO(self.test_csv_gene_header.encode("UTF-8")), "test.csv")})

            self.assertEqual(200, response.status_code)

            result_obj = json.loads(response.data.decode())

            # make sure the samples are correct
            self.assertEqual("Sample 1:Sample 2:Sample 3", ":".join(result_obj["sample_names"]))
            self.assertEqual(4, result_obj["n_lines"])
            self.assertEqual("CD19:CD20:MITF", ":".join(result_obj["top_identifiers"]))

            # make sure the data was stored correctly
            self.assertFalse("data_token" in result_obj)
            self.assertTrue("data" in result_obj)

            # CSV is converted to tsv in the application and quotes removed
            self.assertEqual("Gene\t" + self.test_tsv, result_obj["data"])

    def test_incorrect_header(self):
        with app.app.test_client() as client:
            response = client.post("/upload?store=false",
                                   data={"file": (io.BytesIO(self.test_incorrect_header.encode("UTF-8")), "test.csv")})

            self.assertEqual(400, response.status_code)
            result_obj = json.loads(response.data)

            self.assertEqual("Different number of column names than entries in row 1: header contains 5 fields, "
                             "first line contains 4 fields", result_obj["detail"])

    def test_analysis(self):
        with app.app.test_client() as client:
            response = client.post("/upload?store=true",
                                   data={"file": (io.BytesIO(self.test_tsv.encode("UTF-8")), "test.csv")})
            self.assertEqual(200, response.status_code)
            token = json.loads(response.data)["data_token"]
            self.assertIsNotNone(token)

            request_json_string = '{"methodName": "Camera", "datasets": [{' \
                                  '"data":"' + token + '",' \
                                  '"design": {' \
                                  '"analysisGroup": ["A", "A", "B"], ' \
                                  '"comparison": {"group1": "A", "group2": "B"},' \
                                  '"samples": ["Sample 1", "Sample 2", "Sample 3"]' \
                                  '},' \
                                  '"name": "storedResult", "type": "rnaseq_counts"' \
                                  '}]}' \

            # create the request object
            analysis_response = client.post("/0.1/analysis", data=request_json_string,
                                            content_type="application/json")

            self.assertEqual(200, analysis_response.status_code)


if __name__ == '__main__':
    unittest.main()
