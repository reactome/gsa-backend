import json

import grein_loader
import requests
import pandas
import io
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
            # if set, try to load the dataset from the grein_proxy
            use_grein_proxy = os.getenv("USE_GREIN_PROXY", "False").lower() == "true"

            if use_grein_proxy:
                description, metadata, count_matrix = self._load_dataset_from_proxy(identifier=identifier)

            # in case loading didn't work, try the grein_loader instead
            if not use_grein_proxy or description is None:
                description, metadata, count_matrix = grein_loader.load_dataset(identifier)
        except Exception as e:
            LOGGER.error("Failed to load data for {}".format(identifier))
            LOGGER.exception(e)
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
            metadata_obj = None  # in case metadata is not defined count matrix still will be returned

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

    def _load_dataset_from_proxy(self, identifier: str) -> Tuple:
        """Load the dataset from the internal grein_proxy

        :param identifier: The dataset to load
        :type identifier: str
        :return: A tuple containing the dataset's description, metadata, and counts. All None if
                 the dataset could not be loaded from grein_loader
        :rtype: Tuple
        """
        try:
            LOGGER.debug("Loading dataset %s from grein-proxy...", identifier)

            # description should contain Title and Summary
            status_res = requests.get(f"http://grein-proxy/{ identifier }/status")
            status_res.raise_for_status()
            status = status_res.json()

            if status["status"] != 1:
                raise Exception("Unknown dataset")

            # get the metadata
            meta_res = requests.get(f"http://grein-proxy/{ identifier }/metadata.json")
            meta_res.raise_for_status()

            metadata = meta_res.json()

            # finally, get the count data
            count_res = requests.get(f"http://grein-proxy/{ identifier }/raw_counts.tsv")
            count_res.raise_for_status()

            # create the pandas data.frame
            counts = pandas.read_csv(io.StringIO(count_res.content.decode()), sep="\t")

            LOGGER.debug("Dataset %s successfully loaded from grein-proxy", identifier)

            return ({"Title": status["title"]}, metadata, counts)
        except Exception as e:
            # essentially, ignore all issues
            LOGGER.debug("Failed to load dataset from grein_proxy: %s", str(e))
            return (None, None, None)

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
                    if '' not in list_value_data:
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
