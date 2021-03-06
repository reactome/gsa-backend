# coding: utf-8

from __future__ import absolute_import
from datetime import date, datetime  # noqa: F401

from typing import List, Dict  # noqa: F401

from reactome_analysis_api.models.base_model_ import Model
from reactome_analysis_api.models.design import Design  # noqa: F401,E501
from reactome_analysis_api.models.parameter import Parameter  # noqa: F401,E501
from reactome_analysis_api import util


class Dataset(Model):
    """NOTE: This class is auto generated by the swagger code generator program.

    Do not edit the class manually.
    """

    def __init__(self, name: str=None, type: str=None, data: str=None, design: Design=None, parameters: List[Parameter]=None):  # noqa: E501
        """Dataset - a model defined in Swagger

        :param name: The name of this Dataset.  # noqa: E501
        :type name: str
        :param type: The type of this Dataset.  # noqa: E501
        :type type: str
        :param data: The data of this Dataset.  # noqa: E501
        :type data: str
        :param design: The design of this Dataset.  # noqa: E501
        :type design: Design
        :param parameters: The parameters of this Dataset.  # noqa: E501
        :type parameters: List[Parameter]
        """
        self.swagger_types = {
            'name': str,
            'type': str,
            'data': str,
            'design': Design,
            'parameters': List[Parameter]
        }

        self.attribute_map = {
            'name': 'name',
            'type': 'type',
            'data': 'data',
            'design': 'design',
            'parameters': 'parameters'
        }

        self._name = name
        self._type = type
        self._data = data
        self._design = design
        self._parameters = parameters

    @classmethod
    def from_dict(cls, dikt) -> 'Dataset':
        """Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The Dataset of this Dataset.  # noqa: E501
        :rtype: Dataset
        """
        return util.deserialize_model(dikt, cls)

    @property
    def name(self) -> str:
        """Gets the name of this Dataset.

        Every dataset must have a unique name.  # noqa: E501

        :return: The name of this Dataset.
        :rtype: str
        """
        return self._name

    @name.setter
    def name(self, name: str):
        """Sets the name of this Dataset.

        Every dataset must have a unique name.  # noqa: E501

        :param name: The name of this Dataset.
        :type name: str
        """
        if name is None:
            raise ValueError("Invalid value for `name`, must not be `None`")  # noqa: E501

        self._name = name

    @property
    def type(self) -> str:
        """Gets the type of this Dataset.

        Specifies the type of dataset. Currently supported types are RNA-seq (raw read counts), intensity-based proteomics quantification (proteomics-int), raw proteomics spectral counts (proteomics-sc), and microarray data.  # noqa: E501

        :return: The type of this Dataset.
        :rtype: str
        """
        return self._type

    @type.setter
    def type(self, type: str):
        """Sets the type of this Dataset.

        Specifies the type of dataset. Currently supported types are RNA-seq (raw read counts), intensity-based proteomics quantification (proteomics-int), raw proteomics spectral counts (proteomics-sc), and microarray data.  # noqa: E501

        :param type: The type of this Dataset.
        :type type: str
        """
        allowed_values = ["rnaseq_counts", "rnaseq_norm", "proteomics_int", "proteomics_sc", "microarray_norm"]  # noqa: E501
        if type not in allowed_values:
            raise ValueError(
                "Invalid value for `type` ({0}), must be one of {1}"
                .format(type, allowed_values)
            )

        self._type = type

    @property
    def data(self) -> str:
        """Gets the data of this Dataset.

        Tab-delimited expression matrix with the first column containing gene / protein identifiers, the first row containing the sample labels and each subsequent row corresponding to the expression of one gene in all samples. The 'tab' character must be escaped using '\\t' and new-lines must be escaped using '\\n'. If multiple datasets are submitted, shared samples between the datasets must contain identical labels.  # noqa: E501

        :return: The data of this Dataset.
        :rtype: str
        """
        return self._data

    @data.setter
    def data(self, data: str):
        """Sets the data of this Dataset.

        Tab-delimited expression matrix with the first column containing gene / protein identifiers, the first row containing the sample labels and each subsequent row corresponding to the expression of one gene in all samples. The 'tab' character must be escaped using '\\t' and new-lines must be escaped using '\\n'. If multiple datasets are submitted, shared samples between the datasets must contain identical labels.  # noqa: E501

        :param data: The data of this Dataset.
        :type data: str
        """
        if data is None:
            raise ValueError("Invalid value for `data`, must not be `None`")  # noqa: E501

        self._data = data

    @property
    def design(self) -> Design:
        """Gets the design of this Dataset.


        :return: The design of this Dataset.
        :rtype: Design
        """
        return self._design

    @design.setter
    def design(self, design: Design):
        """Sets the design of this Dataset.


        :param design: The design of this Dataset.
        :type design: Design
        """

        self._design = design

    @property
    def parameters(self) -> List[Parameter]:
        """Gets the parameters of this Dataset.


        :return: The parameters of this Dataset.
        :rtype: List[Parameter]
        """
        return self._parameters

    @parameters.setter
    def parameters(self, parameters: List[Parameter]):
        """Sets the parameters of this Dataset.


        :param parameters: The parameters of this Dataset.
        :type parameters: List[Parameter]
        """

        self._parameters = parameters
