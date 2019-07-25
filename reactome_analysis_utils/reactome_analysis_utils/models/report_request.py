"""
Class containing a report object
"""


import json


class ReportRequest:
    def __init__(self, analysis_id: str, user_mail: str = None):
        """
        Creates a new ReportRequest instance
        :param analysis_id: Analysis id
        :param user_mail: The user's e-mail address (if available)
        """
        self.analysis_id = analysis_id
        self.user_mail = user_mail

    def to_json(self) -> str:
        """
        Serialize to JSON
        :return: The JSON string
        """
        dict_obj = {'analysis_id': self.analysis_id}

        if self.user_mail:
            dict_obj["user_mail"] = self.user_mail

        json_str = json.dumps(dict_obj)

        return json_str


def from_json(json_string: str) -> ReportRequest:
    """
    Creates a ReportRequest object from the passed JSON-encoded string.
    :param json_string: The JSON encoded string
    :return: The generated ReportRequest object.
    """
    dict_obj = json.loads(json_string)

    # make sure the required parameters are present
    required_fields = ["analysis_id"]

    for field in required_fields:
        if field not in dict_obj:
            raise Exception("JSON string does not represent a ReportRequest object. Missing " + field)

    # create the object
    request_obj = ReportRequest(analysis_id=dict_obj["analysis_id"], user_mail=dict_obj.get("user_mail", None))

    return request_obj
