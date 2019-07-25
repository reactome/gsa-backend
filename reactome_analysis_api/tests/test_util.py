import unittest
import json
from reactome_analysis_api.models.analysis_input import AnalysisInput
from reactome_analysis_api.models.analysis_result import AnalysisResult
from reactome_analysis_api.models.analysis_result_results import AnalysisResultResults
from reactome_analysis_api.models.analysis_result_mappings import AnalysisResultMappings
import logging


LOGGER = logging.getLogger(__name__)


class TestUtil(unittest.TestCase):
    def setUp(self):
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
        self.result_json = """
        
        """
        logging.basicConfig(level=logging.DEBUG)

    def test_analysis_request(self):
        LOGGER.debug("Starting tests...")
        request = AnalysisInput.from_dict(json.loads(self.request_json))

        self.assertIsNotNone(request)

    def test_analysis_result(self):
        data_result = AnalysisResultResults(name="Test", pathways="Name\\tFDR\\tFC\\n", fold_changes="MAS2\\tasda\\n")
        mapping_result = AnalysisResultMappings(identifier="MS4A1", mapped_to=["1", "2", "3"])

        result = AnalysisResult(release="1", results=[data_result], mappings=[mapping_result])

        dict = result.to_dict()

        # recreate the object
        new_result = AnalysisResult.from_dict(dict)

        self.assertIsNotNone(new_result)
        self.assertIsNotNone(new_result.results, "Results missing")
        self.assertIsNotNone(new_result.mappings, "Mapping results missing")
        self.assertEqual("1", new_result.release)

        self.assertEqual(1, len(new_result.results))
        self.assertEqual("Name\\tFDR\\tFC\\n", new_result.results[0].pathways)
