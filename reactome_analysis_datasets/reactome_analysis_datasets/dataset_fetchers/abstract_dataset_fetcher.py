from reactome_analysis_api.models.external_data import ExternalData


class DatasetFetcher:
    """
    Abstract base class of all dataset fetchers
    """
    def load_dataset(self, identifier: str) -> (str, ExternalData):
        """
        Loads the specified dataset.
        :param identifier: The dataset's identifier to load.
        :returns: (data, summary) The data as a tab-delimited formatted string, and the summary
                  as a ExternalDataset object
        """
        raise NotImplementedError
