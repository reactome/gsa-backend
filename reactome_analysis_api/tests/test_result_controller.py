# coding: utf-8

from __future__ import absolute_import

from . import base_test_case


class TestResultController(base_test_case.BaseTestCase):
    """ResultController integration tests stubs"""

    def test_get_result(self):
        """Test case for get_result

        Retrieves the result for the completed analysis task
        """
        response = self.client.open(
            '/api/result/{analysisId}'.format(analysisId='analysisId_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_get_status(self):
        """Test case for get_status

        Retrieves the status for the specified analysis.
        """
        response = self.client.open(
            '/api/status/{analysisId}'.format(analysisId='analysisId_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
