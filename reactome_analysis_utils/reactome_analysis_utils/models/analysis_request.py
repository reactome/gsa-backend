"""
Simple class to represent a dataset analysis request
message on RabbitMQ
"""

import json


class AnalysisRequest:
    """
    Represents an analysis request as sent through RabbitMQ.
    """
    def __init__(self, request_id: str):
        """
        Creates a new AnalysisRequest object
        :param request_id: The analysis request id
        """
        self.request_id = request_id

    def to_json(self) -> str:
        """
        Serialize to JSON
        :return: The JSON string
        """
        dict_obj = {'request_id': self.request_id}

        json_str = json.dumps(dict_obj)

        return json_str


def from_json(json_string: str) -> AnalysisRequest:
    """
    Creates an AnalysisRequest object from the passed JSON-encoded string.
    :param json_string: The JSON encoded string
    :return: The generated AnalysisRequest object.
    """
    dict_obj = json.loads(json_string)

    # make sure the required parameters are present
    required_fields = ["request_id"]

    for field in required_fields:
        if field not in dict_obj:
            raise Exception("JSON string does not represent a DatasetRequest object. Missing " + field)

    # create the object
    request_obj = AnalysisRequest(request_id=dict_obj["request_id"])

    return request_obj
