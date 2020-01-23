from reactome_analysis_api.models.external_data import ExternalData
from reactome_analysis_utils import reactome_mq


class DatasetFetcher:
    """
    Abstract base class of all dataset fetchers
    """
    def load_dataset(self, parameters: list, reactome_mq: reactome_mq.ReactomeMQ) -> (str, ExternalData):
        """
        Loads the specified dataset.
        :param parameters: A list of DatasetRequestParameter objects.
        :param reactome_mq: The MQ used to process messages.
        :returns: (data, summary) The data as a tab-delimited formatted string, and the summary
                  as a ExternalDataset object
        """
        raise NotImplementedError

    def _get_parameter(self, name: str, parameters: list) -> str:
        """
        Retrieves the parameter with the specified name from the list
        of passed parameters. Returns None if the parameter does not exist.
        :param name: The parameter's name
        :param parameters: A list of DatasetRequestParameter objects
        :returns: The parameter's value or None if it doesn't exist
        """
        for param in parameters:
            if param.name == name:
                return param.value

        return None

    def get_dataset_id(self, parameters: list) -> str:
        """
        Return a unique identifier for this dataset to check whether the dataset
        already exists in the storage.
        :param parameters: A list of DatasetRequestParameter objects.
        :returns: The datsaet identifier
        """
        raise NotImplementedError

class DatasetFetcherException(Exception):
    """
    Container for DatasetFetcher related exceptions
    """
    pass
