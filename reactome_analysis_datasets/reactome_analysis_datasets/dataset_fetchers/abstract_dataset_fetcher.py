import typing
from reactome_analysis_api.models.external_data import ExternalData
from reactome_analysis_utils import reactome_mq


class DatasetFetcher:
    def __init__(self):
        self._status_callback = None

    """
    Abstract base class of all dataset fetchers
    """

    def load_dataset(self, parameters: list, reactome_mq: reactome_mq.ReactomeMQ) -> typing.Tuple[str, ExternalData]:
        """
        Loads the specified dataset.
        :param parameters: A list of DatasetRequestParameter objects.
        :param reactome_mq: The MQ used to process messages.
        :returns: (data, summary) The data as a tab-delimited formatted string, and the summary
                  as a ExternalDataset object
        """
        raise NotImplementedError
    
    def get_available_datasets(self, n_max_datasets: int = 0) -> list:
        """Get all available datasets from the dataset fetcher. This is only available
           for certain resources. Resources that do not support this return an empty list.

        :param n_max_datasets: Maximum number of datasets to return. If set to 0 (default) 
                               all available datasets are returned, defaults to 0
        :type n_max_datasets: int, optional
        :return: The list of available datasets as a list of dicts. Thest must contain the
                 following keys: "id", "description", "species", "n_samples", 
                 "loading_params" (this is a JSON encoded string with "name"="value" pairs 
                 containing all parameters to load the dataset from the resource)
                 Other keys are allowed but specific to the respective resource.
        :rtype: list
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

        The return value is a list of dicts. Each dict represents one available
        dataset. Required keys are: 
          * id
          * title
          * study_summary (maybe an empty string)
          * species
          * no_samples
          * technology

        :param no_datasets: number of datasets requested
        :returns: list of datasets
        """
        raise NotImplementedError


class DatasetFetcherException(Exception):
    """
    Container for DatasetFetcher related exceptions
    """
    pass
