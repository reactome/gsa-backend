from reactome_analysis_api.models.external_data import ExternalData
from reactome_analysis_utils import reactome_mq


class DatasetFetcher:
    """
    Abstract base class of all dataset fetchers
    """
    def load_dataset(self, identifier: str, reactome_mq: reactome_mq.ReactomeMQ) -> (str, ExternalData):
        """
        Loads the specified dataset.
        :param identifier: The dataset's identifier to load.
        :param reactome_mq: The MQ used to process messages.
        :returns: (data, summary) The data as a tab-delimited formatted string, and the summary
                  as a ExternalDataset object
        """
        raise NotImplementedError


class DatasetFetcherException(Exception):
    """
    Container for DatasetFetcher related exceptions
    """
    pass
