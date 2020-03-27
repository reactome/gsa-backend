import unittest
from reactome_analysis_utils.models import analysis_request


class AnalysisRequestTest(unittest.TestCase):
    def test_create_convert(self):
        org_request = analysis_request.AnalysisRequest(request_id="request_1")

        json_string = org_request.to_json()

        new_request = analysis_request.from_json(json_string)

        self.assertEqual(org_request.request_id, new_request.request_id)
