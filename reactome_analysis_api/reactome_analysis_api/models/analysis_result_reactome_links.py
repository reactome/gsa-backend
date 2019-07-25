# coding: utf-8

from __future__ import absolute_import

from reactome_analysis_api import util
from reactome_analysis_api.models.base_model_ import Model


class AnalysisResultReactomeLinks(Model):
    """NOTE: This class is auto generated by the swagger code generator program.

    Do not edit the class manually.
    """

    def __init__(self, url: str=None, name: str=None, token: str=None, description: str=None):  # noqa: E501
        """AnalysisResultReactomeLinks - a model defined in Swagger

        :param url: The url of this AnalysisResultReactomeLinks.  # noqa: E501
        :type url: str
        :param name: The name of this AnalysisResultReactomeLinks.  # noqa: E501
        :type name: str
        :param token: The token of this AnalysisResultReactomeLinks.  # noqa: E501
        :type token: str
        :param description: The description of this AnalysisResultReactomeLinks.  # noqa: E501
        :type description: str
        """
        self.swagger_types = {
            'url': str,
            'name': str,
            'token': str,
            'description': str
        }

        self.attribute_map = {
            'url': 'url',
            'name': 'name',
            'token': 'token',
            'description': 'description'
        }

        self._url = url
        self._name = name
        self._token = token
        self._description = description

    @classmethod
    def from_dict(cls, dikt) -> 'AnalysisResultReactomeLinks':
        """Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The AnalysisResult_reactome_links of this AnalysisResultReactomeLinks.  # noqa: E501
        :rtype: AnalysisResultReactomeLinks
        """
        return util.deserialize_model(dikt, cls)

    @property
    def url(self) -> str:
        """Gets the url of this AnalysisResultReactomeLinks.

        Link to the result visualization in the Reactome pathway browser  # noqa: E501

        :return: The url of this AnalysisResultReactomeLinks.
        :rtype: str
        """
        return self._url

    @url.setter
    def url(self, url: str):
        """Sets the url of this AnalysisResultReactomeLinks.

        Link to the result visualization in the Reactome pathway browser  # noqa: E501

        :param url: The url of this AnalysisResultReactomeLinks.
        :type url: str
        """
        if url is None:
            raise ValueError("Invalid value for `url`, must not be `None`")  # noqa: E501

        self._url = url

    @property
    def name(self) -> str:
        """Gets the name of this AnalysisResultReactomeLinks.

        Short name of the type of visualization  # noqa: E501

        :return: The name of this AnalysisResultReactomeLinks.
        :rtype: str
        """
        return self._name

    @name.setter
    def name(self, name: str):
        """Sets the name of this AnalysisResultReactomeLinks.

        Short name of the type of visualization  # noqa: E501

        :param name: The name of this AnalysisResultReactomeLinks.
        :type name: str
        """
        if name is None:
            raise ValueError("Invalid value for `name`, must not be `None`")  # noqa: E501

        self._name = name

    @property
    def token(self) -> str:
        """Gets the token of this AnalysisResultReactomeLinks.

        The token of the Reactome analysis  # noqa: E501

        :return: The token of this AnalysisResultReactomeLinks.
        :rtype: str
        """
        return self._token

    @token.setter
    def token(self, token: str):
        """Sets the token of this AnalysisResultReactomeLinks.

        The token of the Reactome analysis  # noqa: E501

        :param token: The token of this AnalysisResultReactomeLinks.
        :type token: str
        """
        if token is None:
            raise ValueError("Invalid value for `token`, must not be `None`")  # noqa: E501

        self._token = token

    @property
    def description(self) -> str:
        """Gets the description of this AnalysisResultReactomeLinks.

        A description of the visualization type.  # noqa: E501

        :return: The description of this AnalysisResultReactomeLinks.
        :rtype: str
        """
        return self._description

    @description.setter
    def description(self, description: str):
        """Sets the description of this AnalysisResultReactomeLinks.

        A description of the visualization type.  # noqa: E501

        :param description: The description of this AnalysisResultReactomeLinks.
        :type description: str
        """

        self._description = description
