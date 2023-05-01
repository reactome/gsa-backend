import grein_loader
import logging
import os

from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, DatasetFetcherException



LOGGER = logging.getLogger(__name__)


class GreinFetcher(DatasetFetcher):
    """
    A DatasetFetcher to retrieve data from GREIN
    """

    # URL to retrieve data from GREIN
    # http://www.ilincs.org/apps/grein/?gse=

    # metadata, description, count matrix are provided via streaming within the grein_loader package
    # description, metadata, count_matrix = grein_loader.load_dataset(geo_accession)
    # Example for geo_accesssion GSE112749

    def __int__(self):
        # constructor of abstract super class
        super().__init__()
        self.max_timeout = int(os.getenv("LOADING_MAX_TIMEOUT", 60))

    def get_dataset_id(self, parameters: list) -> str:
        """
        Returns the dataset identifier
        :param parameters: A list of DatasetRequestParameter objects.
        :returns: The dataset identifier
        """
        return self._get_parameter(name="dataset_id", parameters=parameters)

    def load_dataset(self, parameters: list, reactome_mq) -> (
            str, ExternalData):
        """
        Load the specified GREIN dataset based on GSE id
        :param parameters: GSE id from GREIN
        :param reactome_mq: Not used
        :returns: (metadata, count_matrix)
        """

        # get id parameter
        identifier = self._get_parameter(name="dataset_id", parameters=parameters)

        # check for identifier
        if not identifier:
            raise DatasetFetcherException(
                "Missing required parameter 'dataset_id' to load the example dataset.")

        # check for correct format in geo_accession
        if not identifier[0:3] == "GSE":
            raise DatasetFetcherException(
                f"{identifier} is not a valid geo_accession for GREIN")

        # load the data
        self._update_status(progress=0.2, message="Requesting data from GREIN")

        try:
            description, metadata, count_matrix = grein_loader.load_dataset(identifier)
        except Exception:
            raise DatasetFetcherException("Failed to load data for {}".format(identifier))

        self._update_status(progress=0.7, message="Converting metadata")

        try:
            metadata_obj = ExternalData.from_dict(metadata)
        except Exception:
            raise DatasetFetcherException(
                "Failed to load a valid summary for {}".format(identifier))

        self._update_status(progress=0.8, message="Converting count matrix")
        try:
            count_matrix_tsv = count_matrix.to_csv(sep="\t")
        except Exception:
            raise DatasetFetcherException(
                "Failed to load a valid summary for {}".format(identifier))

        self._update_status(progress=0.8, message="Creating summary data")

        # return data
        return (count_matrix_tsv, metadata_obj)

    def get_available_datasets(self, no_datasets: int) -> list:
        """
        Loads overview of GREIN datasets
        :param no_datasets: Number of datasets loading from GREIN
        :returns: list of datasets with description
        """
        return grein_loader.load_overview(no_datasets)
