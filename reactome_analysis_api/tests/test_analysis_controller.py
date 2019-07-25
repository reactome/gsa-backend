# coding: utf-8

import json
import os
import unittest

import sys

sys.path.insert(0, os.path.dirname(__file__))

from reactome_analysis_api.controllers import analysis_controller
from base_test_case import BaseTestCase


class TestAnalysisController(BaseTestCase):
    def setUp(self):
        os.environ["REDIS_HOST"] = "192.168.99.100"
        os.environ["REDIS_PORT"] = "31727"
        os.environ["REDIS_PASSWORD"] = "test"
        os.environ["RABBIT_HOST"] = "192.168.99.100"
        os.environ["RABBIT_PORT"] = "32528"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

        self.request_json = """
        {
  "analysisId": "analysisId",
  "datasets": [
    {
      "data": "\\tSample 1\\tSample2\\tSample 3\\nGene 1\\t10\\t20\\t2\\n",
      "design": {
        "analysisGroup": [
          "Treatment",
          "Control",
          "Treatment"
        ],
        "comparison": {
          "group1": "Control",
          "group2": "Treatment"
        },
        "samples": [
          "Sample 1",
          "Sample 2",
          "Sample 3"
        ],
        "patient": [
          "Patient 1",
          "Patient 2",
          "Patient 3"
       ]
      },
      "name": "First experiment",
      "type": "rnaseq"
    }
  ],
  "methodName": "Camera",
  "parameters": [
    {
      "name": "permutations",
      "value": "10"
    },
    {
      "name": "permutations",
      "value": "10"
    }
  ]
}
"""

    def test_invalid_json(self):
        response = self.client.open(
            '/0.1/analysis',
            method='POST',
            data="Some bad string",
            content_type='application/json')

        self.assertEqual(415, response.status_code)
        self.assertEqual('415 UNSUPPORTED MEDIA TYPE', response.status)

    def test_invalid_input(self):
        invalid_json = self.request_json.replace("data", "asda")

        response = self.client.open(
            '/0.1/analysis',
            method='POST',
            data=invalid_json,
            content_type='application/json')

        self.assertEqual(400, response.status_code)

    def test_invalid_method(self):
        request_dict = json.loads(self.request_json)
        request_dict["methodName"] = "NonExisting"

        response = self.client.open(
            '/0.1/analysis',
            method='POST',
            data=json.dumps(request_dict),
            content_type='application/json')

        self.assert404(response, 'Response body is : ' + response.data.decode('utf-8'))

    def test_failed_queuing_system(self):
        org_address = os.environ["RABBIT_HOST"]
        os.environ["RABBIT_HOST"] = "not.here"

        response = self.client.open(
            '/0.1/analysis',
            method='POST',
            data=self.request_json,
            content_type='application/json')

        self.assertEqual(503, response.status_code)

        os.environ["RABBIT_HOST"] = org_address

    def test_failed_storage(self):
        org_address = os.environ["REDIS_HOST"]
        os.environ["REDIS_HOST"] = "redis.not.here"

        response = self.client.open(
            '/0.1/analysis',
            method='POST',
            data=self.request_json,
            content_type='application/json')

        self.assertEqual(503, response.status_code)
        os.environ["REDIS_HOST"] = org_address

    def test_fix_additional_properties(self):
        request_object = json.loads(self.request_json)

        self.assertEqual(1, len(request_object["datasets"]))
        self.assertTrue("patient" in request_object["datasets"][0]["design"])

        input_object = analysis_controller.create_analysis_input_object(request_object)

        self.assertTrue(hasattr(input_object.datasets[0].design, "additional_properties"))

    def test_duplicate_dataset_names(self):
        request_json = """
                {
          "analysisId": "analysisId",
          "datasets": [
            {
              "data": "\\tSample 1\\tSample2\\tSample 3\\nGene 1\\t10\\t20\\t2\\n",
              "design": {
                "analysisGroup": [
                  "Treatment",
                  "Control",
                  "Treatment"
                ],
                "comparison": {
                  "group1": "Control",
                  "group2": "Treatment"
                },
                "samples": [
                  "Sample 1",
                  "Sample 2",
                  "Sample 3"
                ],
                "patient": [
                  "Patient 1",
                  "Patient 2",
                  "Patient 3"
               ]
              },
              "name": "First experiment",
              "type": "rnaseq"
            },
            {
              "data": "\\tSample 1\\tSample2\\tSample 3\\nGene 1\\t10\\t20\\t2\\n",
              "design": {
                "analysisGroup": [
                  "Treatment",
                  "Control",
                  "Treatment"
                ],
                "comparison": {
                  "group1": "Control",
                  "group2": "Treatment"
                },
                "samples": [
                  "Sample 1",
                  "Sample 2",
                  "Sample 3"
                ],
                "patient": [
                  "Patient 1",
                  "Patient 2",
                  "Patient 3"
               ]
              },
              "name": "First experiment",
              "type": "rnaseq"
            }
          ],
          "methodName": "Camera",
          "parameters": [
            {
              "name": "permutations",
              "value": "10"
            },
            {
              "name": "permutations",
              "value": "10"
            }
          ]
        }
        """
        request_object = json.loads(request_json)

        self.assertEqual(2, len(request_object["datasets"]))

        response = self.client.open(
            '/0.1/analysis',
            method='POST',
            data=request_json,
            content_type='application/json')

        self.assertEqual(406, response.status_code)

    def test_create_input_object_default_params(self):
        input_dict = {
            "datasets": [],
            "methodName": "camera",
            "parameters": [
                {"name": "max_missing_values", "value": "0.1"}
            ]
        }

        input_object = analysis_controller.create_analysis_input_object(input_dict)

        self.assertIsNotNone(input_object)
        self.assertEqual(4, len(input_object.parameter_dict))
        self.assertEqual("0.1", input_object.parameter_dict["max_missing_values"])

        input_dict = {
            "datasets": [],
            "methodName": "camera",
            "parameters": [ ]
        }

        input_object = analysis_controller.create_analysis_input_object(input_dict)

        self.assertIsNotNone(input_object)
        self.assertEqual(4, len(input_object.parameter_dict))
        self.assertEqual("0.5", input_object.parameter_dict["max_missing_values"])


if __name__ == '__main__':
    unittest.main()
