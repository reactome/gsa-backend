from reactome_analysis_api.models.external_data import ExternalData
from reactome_analysis_utils import reactome_mq


class DatasetFetcher:
    def __init__(self):
        self._status_callback = None

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

    def set_status_callback(self, status_callback) -> None:
        """
        Set the function to use to update status messages.

        The function must follow the signature (progress: float (0-1), message: str)
        :param status_callback: The function to use for status updates
        """
        self._status_callback = status_callback

    def _update_status(self, progress: float, message: str) -> None:
        """
        Update the current status.

        :param progress: The progress as a float between 0 - 1
        :param message: The status message to show
        """
        if self._status_callback:
            self._status_callback(progress=progress, message=message)

    def get_available_datasets(self, no_datasets: int) -> list:
        """
        Loads overview of datasets based on number of datasets
        :param no_datasets: number of datasets requested
        :returns: list of datasets
        """
        raise NotImplementedError


class DatasetFetcherException(Exception):
    """
    Container for DatasetFetcher related exceptions
    """
    pass
