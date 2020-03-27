import unittest
from reactome_analysis_utils.models.dataset_request import DatasetRequest, from_json, DatasetRequestParameter


class DatasetRequestTest(unittest.TestCase):
    def test_create_object(self):
        request_obj = DatasetRequest(resource_id="Example_1", loading_id="L1", 
                                     parameters=[DatasetRequestParameter("p1", "v1"), DatasetRequestParameter("p2", "v2")])

        self.assertEqual("Example_1", request_obj.resource_id)

        # turn into JSON
        json_str = request_obj.to_json()

        # get a new object
        new_obj = from_json(json_str)

        self.assertEqual(request_obj.resource_id, new_obj.resource_id)
        self.assertEqual(request_obj.loading_id, new_obj.loading_id)
        self.assertEqual(2, len(request_obj.parameters))
        self.assertEqual("p1", request_obj.parameters[0].name)
