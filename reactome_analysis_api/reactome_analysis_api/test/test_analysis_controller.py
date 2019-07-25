# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from reactome_analysis_api.models.analysis_input import AnalysisInput  # noqa: E501
from reactome_analysis_api.models.data_type import DataType  # noqa: E501
from reactome_analysis_api.models.method import Method  # noqa: E501
from reactome_analysis_api.test import BaseTestCase


class TestAnalysisController(BaseTestCase):
    """AnalysisController integration test stubs"""

    def test_list_methods(self):
        """Test case for list_methods

        Lists the available analysis methods
        """
        response = self.client.open(
            '/0.1/methods',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_list_types(self):
        """Test case for list_types

        Lists the supported data types
        """
        response = self.client.open(
            '/0.1/types',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_start_analysis(self):
        """Test case for start_analysis

        Performs the specified gene set analysis
        """
        body = AnalysisInput()
        response = self.client.open(
            '/0.1/analysis',
            method='POST',
            data=json.dumps(body),
            content_type='application/json')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
