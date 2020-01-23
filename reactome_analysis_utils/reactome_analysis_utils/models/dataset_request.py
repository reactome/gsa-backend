"""
Simple class to encapsulate a request to retrieve
a new dataset
"""

import json


class DatasetRequest:
    """
    Request to retrieve a new dataset
    """
    def __init__(self, loading_id: str, resource_id: str, parameters: list):
        """
        Initializes a new DatasetRequest object.
        :param loading_id: The identifier of the loading process.
        :param resource_id: The identifier of the resource to load
        :param parameters: A list of DatasetRequestParameter objects
        """
        self.resource_id = resource_id
        self.loading_id = loading_id
        self.parameters = parameters

    def to_json(self) -> str:
        """
        Serialize to JSON
        :return: The JSON string
        """
        dict_obj = {'resource_id': self.resource_id, 'loading_id': self.loading_id, 'parameters': []}

        for parameter in self.parameters:
            dict_obj['parameters'].append({'name': parameter.name, 'value': parameter.value})

        json_str = json.dumps(dict_obj)

        return json_str


def from_json(json_string: str) -> DatasetRequest:
    """
    Creates a DatasetRequest object from the passed JSON-encoded string.
    :param json_string: The JSON encoded string
    :return: The generated DatasetRequest object.
    """
    dict_obj = json.loads(json_string)

    # make sure the required parameters are present
    required_fields = ["resource_id", "loading_id"]

    for field in required_fields:
        if field not in dict_obj:
            raise Exception("JSON string does not represent a DatasetRequest object. Missing " + field)

    # create the parameter object
    parameters = list()

    for dict_param in dict_obj["parameters"]:
        parameters.append(DatasetRequestParameter(name=dict_param["name"], value=dict_param["value"]))

    # create the object
    request_obj = DatasetRequest(resource_id=dict_obj["resource_id"], loading_id=dict_obj["loading_id"], parameters=parameters)

    return request_obj


class DatasetRequestParameter:
    def __init__(self, name: str, value: str):
        self.name = name
        self.value = value
