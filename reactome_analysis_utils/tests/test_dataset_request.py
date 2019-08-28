import unittest
from reactome_analysis_utils.models.dataset_request import DatasetRequest, from_json


class DatasetRequestTest(unittest.TestCase):
    def test_create_object(self):
        request_obj = DatasetRequest(dataset_id="Analaysis_123", loading_id="L1")

        self.assertEqual("Analaysis_123", request_obj.dataset_id)

        # turn into JSON
        json_str = request_obj.to_json()

        # get a new object
        new_obj = from_json(json_str)

        self.assertEqual(request_obj.dataset_id, new_obj.dataset_id)
        self.assertEqual(request_obj.loading_id, new_obj.loading_id)
