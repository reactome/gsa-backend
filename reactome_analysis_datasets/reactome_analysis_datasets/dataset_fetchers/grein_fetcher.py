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
        identifier = self._get_parameter(name="dataset_id", parameters=parameters).strip()

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
            LOGGER.error("Failed to load data for {}".format(identifier))
            raise DatasetFetcherException("Dataset '{}' is not available on GREIN".format(identifier))

        self._update_status(progress=0.7, message="Converting metadata")

        if metadata != '':  # in case metadata is not defined by GREIN
            try:
                metadata_obj = self._create_metadata(metadata, description=description, identifier=identifier)
            except Exception:
                raise DatasetFetcherException(
                    "Failed to load a valid summary for {}".format(identifier))
        else:
            LOGGER.info("No metadata available from GREIN")
            self._update_status(progress=0.75, message="No metadata found")
            metadata_obj = None   # in case metadata is not defined count matrix still will be returned

        self._update_status(progress=0.8, message="Converting count matrix")
        count_matrix = count_matrix.drop("gene_symbol", axis=1)
        try:
            count_matrix_tsv = count_matrix.to_csv(sep="\t", index=False)
        except Exception:
            LOGGER.error("No expression values available on GREIN for {}".format(identifier))

            raise DatasetFetcherException(
                "No expression values available on GREIN for {}".format(identifier))

        self._update_status(progress=1, message="Completed loading  {}".format(identifier))

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

        samples = self._get_sample_ids(metadata)  # gets sample ids for dictionary
        summary['sample_ids'] = samples
        metadata_obj = ExternalData.from_dict(summary)  # converts th
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
        filtered_metadata = self._format_metadata(list_new_metadata)  # formats list based on the values
        filtered_metadata = [ExternalDataSampleMetadata.from_dict(metadata) for metadata in filtered_metadata]
        return filtered_metadata

    def _format_metadata(self, list_metadata):
        """
        formats and filters metadata
        :param list of metadata dictionary
        :return filtered list of dictionaries
        """
        filtered_list = []

        for item in list_metadata:
            entry_value = item['name'].split("_")
            if entry_value[0] == "characteristics":
                value = item['values'][0]
                if value is not None:
                    list_value_data = value.split(": ")
                    item['name'] = list_value_data[0]
                    original_values = item['values']
                    item['values'] = [item.split(": ")[1] for item in original_values]
                    filtered_list.append(item)

        return filtered_list

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

    def get_available_datasets(self, no_datasets: int = None) -> list:
        """
        Loads overview of GREIN datasets
        :param no_datasets: Number of datasets loading from GREIN
        :returns: list of datasets with description
        """
        return grein_loader.load_overview(no_datasets)
