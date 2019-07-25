# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from reactome_analysis_api.models.analysis_result import AnalysisResult  # noqa: E501
from reactome_analysis_api.models.analysis_status import AnalysisStatus  # noqa: E501
from reactome_analysis_api.test import BaseTestCase


class TestResultController(BaseTestCase):
    """ResultController integration test stubs"""

    def test_get_result(self):
        """Test case for get_result

        Retrieves the result for the completed analysis task
        """
        response = self.client.open(
            '/0.1/result/{analysisId}'.format(analysisId='analysisId_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_get_status(self):
        """Test case for get_status

        Retrieves the status for the specified analysis.
        """
        response = self.client.open(
            '/0.1/status/{analysisId}'.format(analysisId='analysisId_example'),
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
