# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from reactome_analysis_api.models.external_data import ExternalData  # noqa: E501
from reactome_analysis_api.test import BaseTestCase


class TestDatasetsController(BaseTestCase):
    """DatasetsController integration test stubs"""

    def test_get_examples(self):
        """Test case for get_examples

        Lists the available example datasets
        """
        response = self.client.open(
            '/0.1/data/examples',
            method='GET')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
