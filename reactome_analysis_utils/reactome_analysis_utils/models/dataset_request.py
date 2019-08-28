"""
Simple class to encapsulate a request to retrieve
a new dataset
"""

import json


class DatasetRequest:
    """
    Request to retrieve a new dataset
    """
    def __init__(self, loading_id: str, dataset_id: str):
        """
        Initializes a new DatasetRequest object.
        :param loading_id: The identifier of the loading process.
        :param dataset_id: The identifier of the dataset to load.
        """
        self.dataset_id = dataset_id
        self.loading_id = loading_id

    def to_json(self) -> str:
        """
        Serialize to JSON
        :return: The JSON string
        """
        dict_obj = {'dataset_id': self.dataset_id, 'loading_id': self.loading_id}

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
    required_fields = ["dataset_id", "loading_id"]

    for field in required_fields:
        if field not in dict_obj:
            raise Exception("JSON string does not represent a DatasetRequest object. Missing " + field)

    # create the object
    request_obj = DatasetRequest(dataset_id=dict_obj["dataset_id"], loading_id=dict_obj["loading_id"])

    return request_obj