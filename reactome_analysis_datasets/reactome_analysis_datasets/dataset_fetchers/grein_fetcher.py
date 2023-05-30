import grein_loader
import logging
import os
from typing import Tuple

from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, \
    DatasetFetcherException
from reactome_analysis_api.models.external_data_sample_metadata import ExternalDataSampleMetadata

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

    def load_dataset(self, parameters: list, reactome_mq) -> Tuple[str, ExternalData]:
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
            metadata_obj = self._create_metadata(metadata, description=description, identifier=identifier)
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

    def _create_metadata(self, metadata, description, identifier: str) -> ExternalData:
        """
        fetches the data in ExternalData object
        :param metadata loaded by the GREIN plugin
        :param description Dict returned by the GREIN plugin
        :param identifier The dataset's original identifier
        :returns ExternalData object
        """
        # initialize the summary
        summary = {"type": "rnaseq_counts", "id": identifier, "title": "Public data from GREIN",
                   "description": "Public dataset from Grein",
                   "sample_ids": list()
                   }
        
        # change to a nice title if available
        if "Title" in description and description["Title"]:
            summary["title"] = "GREIN dataset " + identifier + ": " + description["Title"]
        if "Summary" in description and description["Summary"]:
            summary["description"] = description["Summary"]

        samples = self._get_sample_ids(metadata)
        summary['sample_ids'].append(samples)
        metadata_obj = ExternalData.from_dict(summary)
        list_metadata = self._get_metadata(metadata)
        metadata_obj.sample_metadata = list_metadata  # adds metadata via setter in the object
        return metadata_obj

    def _get_metadata(self, list_metadata):
        """
        gets metadata for each sample, as dictionary with values and a list of the metadat for each sample
        :param list_metadata list of the requested metadata provided
        :returns list of ExternalDataSampleMetadata
        """
        merged_dict = {}
        for dictionary in list_metadata.values():
            for key, value in dictionary.items():
                merged_dict.setdefault(key, []).append(value)

        list_new_metadata = [{'name': key, 'values': values} for key, values in merged_dict.items()]
        list_new_metadata = [ExternalDataSampleMetadata.from_dict(metadata) for metadata in list_new_metadata]
        return list_new_metadata

    def _get_sample_ids(self, list_metadata) -> list:
        """
        gets sample ids of the data
        :param list_metadata list of the requested metadata provided
        :returns list of the sample id as a string
        """
        sample_id_list = []
        for key, data in list_metadata.items():
            sample_id_list.append(key)
        return sample_id_list

    def get_available_datasets(self, no_datasets: int) -> list:
        """
        Loads overview of GREIN datasets
        :param no_datasets: Number of datasets loading from GREIN
        :returns: list of datasets with description
        """
        return grein_loader.load_overview(no_datasets)