import copy
import json
import unittest

from reactome_analysis_api import input_deserializer


class InputDeserializerTest(unittest.TestCase):
    def setUp(self) -> None:
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
              "type": "rnaseq_counts"
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
        self.request_obj = json.loads(self.request_json)

    def testNoDatasetParams(self):
        input_object = input_deserializer.create_analysis_input_object(self.request_obj)

        self.assertIsNotNone(input_object)

        # make sure the dataset-level parameters are there
        self.assertIsNotNone(input_object.datasets[0].parameter_dict)
        self.assertEqual(3, len(input_object.datasets[0].parameter_dict))

        for param_name in ["max_missing_values", "discrete_norm_function", "continuous_norm_function"]:
            self.assertTrue(param_name in input_object.datasets[0].parameter_dict)

    def testDatasetLevelParam(self):
        request_obj = copy.deepcopy(self.request_obj)
        request_obj["datasets"][0]["parameters"] = [{"name": "max_missing_values", "value": "-1"}]

        input_object = input_deserializer.create_analysis_input_object(request_obj)

        self.assertIsNotNone(input_object)

        # make sure the dataset-level parameters are there
        self.assertIsNotNone(input_object.datasets[0].parameter_dict)
        self.assertEqual(3, len(input_object.datasets[0].parameter_dict))

        for param_name in ["max_missing_values", "discrete_norm_function", "continuous_norm_function"]:
            self.assertTrue(param_name in input_object.datasets[0].parameter_dict)

        self.assertEqual("-1", input_object.datasets[0].parameter_dict["max_missing_values"])
