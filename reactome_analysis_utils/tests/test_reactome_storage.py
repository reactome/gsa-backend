import unittest
import os

from reactome_analysis_utils import reactome_storage


class ReactomeStorageTest(unittest.TestCase):
    def setUp(self):
        os.environ["REDIS_HOST"] = "localhost"
        os.environ["REDIS_PORT"] = "32088"
        os.environ["REDIS_PASSWORD"] = "test"

    def test_connection(self):
        storage = reactome_storage.ReactomeStorage()

        # just tests whether the connection works
        exists = storage.analysis_exists("ABC")

        self.assertFalse(exists)

    def test_compression(self):
        reactome_storage.ReactomeStorage.USE_COMPRSSSION = True

        storage = reactome_storage.ReactomeStorage()

        test_data = '{"name": "Johannes", "age": "38"}'
        test_token = "TEST1"

        # request data - uploaded by the user
        storage.set_request_data(token=test_token, data=test_data)
        fetched_data = storage.get_request_data(token=test_token)

        self.assertEqual(test_data, fetched_data, msg="Request data does not match")

        # result
        storage.set_result(analysis_identifier=test_token, result=test_data, data_type="analysis")
        fetched_data = storage.get_result(analysis_identifier=test_token, data_type="analysis")

        self.assertEqual(test_data, fetched_data, msg="Result data does not match")

        # request data summary
        storage.set_request_data_summary(token=test_token, data=test_data)
        fetched_data = storage.get_request_data_summary(token=test_token)

        self.assertEqual(test_data, fetched_data, msg="Request data summary does not match")

        # analysis request data - the request object
        storage.set_analysis_request_data(token=test_token, data=test_data)
        fetched_data = storage.get_analysis_request_data(token=test_token)

        self.assertEqual(test_data, fetched_data, msg="Analysis request data does not match")

        # result data
        storage.set_result
