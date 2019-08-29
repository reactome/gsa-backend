import io
import json
import os
import unittest
import time

import sys
from reactome_analysis_utils.reactome_storage import ReactomeStorage

sys.path.insert(0, os.path.dirname(__file__))

from reactome_analysis_api.__main__ import app


class ExternalDataTest(unittest.TestCase):
    def setUp(self):
        os.environ["REDIS_HOST"] = "192.168.99.100"
        os.environ["REDIS_PORT"] = "31297"
        os.environ["REDIS_PASSWORD"] = "test"
        os.environ["RABBIT_HOST"] = "192.168.99.100"
        os.environ["RABBIT_PORT"] = "31715"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

        app.app.testing = True

    def test_wrong_identifier(self):
        token = "EXAMPLE_MEL_PROT"

        with app.app.test_client() as client:
            # load the dataset
            load_response = client.get("/0.1/data/load/" + token)
            self.assertEqual(200, load_response.status_code)

            # keep using the dataset as a token
            load_status = client.get("/0.1/data/status/" + token)
            self.assertEqual(404, load_status.status_code)

    def test_external_dataset(self):
        token = "EXAMPLE_MEL_PROT"

        with app.app.test_client() as client:
            # load the dataset
            load_response = client.get("/0.1/data/load/" + token)
            self.assertEqual(200, load_response.status_code)

            load_token = load_response.data.decode()

            # get the status
            load_status = client.get("/0.1/data/status/" + load_token)
            self.assertEqual(200, load_status.status_code)

            # assess the status
            status_obj = json.loads(load_status.data)
            # status updates happen but sometimes do not appear in the test
            self.assertNotEqual("failed", status_obj["status"])

            # get the summary
            summary_response = client.get("/0.1/data/summary/" + token)
            self.assertEqual(200, summary_response.status_code)

            # make sure it's well formatted
            summary_obj = json.loads(summary_response.data)
            self.assertIsNotNone(summary_obj)

            # submit an analysis request
            request_json_string = '{"methodName": "Camera", "datasets": [{' \
                                    '"data":"' + token + '",' \
                                    '"design": {' \
                                    '"analysisGroup": ["A", "A", "B"], ' \
                                    '"comparison": {"group1": "A", "group2": "B"},' \
                                    '"samples": ["Sample 1", "Sample 2", "Sample 3"]' \
                                    '},' \
                                    '"name": "storedResult", "type": "rnaseq_counts"' \
                                    '}]}' \

            analysis_response = client.open('/0.1/analysis',
                                            method='POST',
                                            data=request_json_string,
                                            content_type='application/json')

            self.assertEqual(200, analysis_response.status_code)

    def test_rna_dataset(self):
        token = "EXAMPLE_MEL_RNA"

        with app.app.test_client() as client:
            # load the dataset
            load_response = client.get("/0.1/data/load/" + token)
            self.assertEqual(200, load_response.status_code)

            load_token = load_response.data.decode()

            # get the status
            load_status = client.get("/0.1/data/status/" + load_token)
            self.assertEqual(200, load_status.status_code)

            # assess the status
            status_obj = json.loads(load_status.data)
            # status updates happen but sometimes do not appear in the test
            self.assertNotEqual("failed", status_obj["status"])

            # get the summary
            summary_response = client.get("/0.1/data/summary/" + token)
            self.assertEqual(200, summary_response.status_code)

            # make sure it's well formatted
            summary_obj = json.loads(summary_response.data)
            self.assertIsNotNone(summary_obj)

            # submit an analysis request
            request_json_string = '{"methodName": "Camera", "datasets": [{' \
                                    '"data":"' + token + '",' \
                                    '"design": {' \
                                    '"analysisGroup": ["A", "A", "B"], ' \
                                    '"comparison": {"group1": "A", "group2": "B"},' \
                                    '"samples": ["Sample 1", "Sample 2", "Sample 3"]' \
                                    '},' \
                                    '"name": "storedResult", "type": "rnaseq_counts"' \
                                    '}]}' \

            analysis_response = client.open('/0.1/analysis',
                                            method='POST',
                                            data=request_json_string,
                                            content_type='application/json')

            self.assertEqual(200, analysis_response.status_code)