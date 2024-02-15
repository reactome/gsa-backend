import logging
import GEOparse as geoparser
import os
from typing import Tuple

import requests
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, \
    DatasetFetcherException
from reactome_analysis_api.models.external_data_sample_metadata import ExternalDataSampleMetadata
import pandas as pd

LOGGER = logging.getLogger(__name__)


class GeoFetcher(DatasetFetcher):
    """
    A DatasetFetcher to retrieve data from GEO
    """

    def __init__(self):
        # constructor of abstract super class
        super().__init__()

    def load_dataset(self, parameters: list, reactome_mq) -> Tuple[str, ExternalData]:
        """
        Load the specified GEO dataset based on GSE id
        :param parameters: GSE id from GEO
        :returns empty string insted of count matrix if exists, metadata of all datasets
        """

        identifier = self._get_parameter(name="dataset_id", parameters=parameters)
        if not identifier:
            raise DatasetFetcherException(
                "Missing required parameter 'dataset_id' to load the example dataset.")
        if not identifier[0:3] == "GSE":
            raise DatasetFetcherException(
                f"{identifier} is not a valid geo_accession for GREIN")

        # load the data
        LOGGER.info(f"Loading dataset {identifier} from GEO")
        try:
            gse = geoparser.get_GEO(identifier)  # fetching data from Geo via geo parser
        except Exception as e:
            LOGGER.error(f"Error loading dataset {identifier} from GEO: {e}")
            raise DatasetFetcherException(f"Error loading dataset {identifier} from GEO: {e}")

        # create metadata
        if not gse.metadata:
            LOGGER.error(f"Error loading dataset {identifier} from GEO: No metadata found")
            raise DatasetFetcherException(f"Error loading dataset {identifier} from GEO: No metadata found")
        else:
            LOGGER.info(f"Creating metadata for dataset {identifier}")
            experiment_type = self._get_data_type(gse.metadata)
            sample_metadata_list = self._create_sample_metadata(gse.metadata["sample_id"])
            metadata_obj = ExternalData(id=identifier,
                                        title=gse.metadata["title"],
                                        type=experiment_type,
                                        description=gse.metadata["summary"],
                                        group=None,
                                        sample_ids=gse.metadata["sample_id"],
                                        sample_metadata=sample_metadata_list,
                                        default_parameters=None)

        os.remove(identifier + "_family.soft.gz")  # clean up supplementary

        self._clean_up_samples(sample_metadata_list[1].values)  # clean up of downloaded files

        # create matrix if parameter exists
        filename = self._get_parameter(name="filename", parameters=parameters)
        if not filename:
            LOGGER.warning("No Matrix file give only metadata will be returned")
            return "", metadata_obj
        else:
            geo_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={identifier}&format=file&file={filename}"
            matrix_response = requests.get(geo_url)
            if matrix_response.status_code != 200:
                LOGGER.error(f"Error loading dataset {identifier} from GEO")
                raise DatasetFetcherException(f"Error loading matrix file {filename} from GEO")
            else:
                LOGGER.info(f"Expression values with filename: {filename} fetched from GEO")
                matrix = pd.read_csv(geo_url)
                matrix = matrix.to_string()

        return matrix, metadata_obj

    def _create_sample_metadata(self, sample_list) -> list[ExternalDataSampleMetadata]:
        """
        Create sample metadata for each sample in the dataset
        :param sample_list: list of samples in the dataset
        :returns: list of sample metadata formatted as ExternalDataSampleMetadata
        """
        metadata_request_list = [geoparser.get_GEO(sample).metadata for sample in sample_list]

        metadata_dict = {}

        for sample_metadata in metadata_request_list:
            for key, values in sample_metadata.items():
                if key not in metadata_dict:
                    metadata_dict[key] = []
                metadata_dict[key].extend(values)

        formatted_metadata_list = [
            {"name": key, "values": metadata_dict[key]} for key in metadata_dict
        ]
        filtered_metadata = [
            ExternalDataSampleMetadata.from_dict(metadata) for metadata in formatted_metadata_list
        ]
        return filtered_metadata

    def _get_data_type(self, metadata_obj) -> str:
        """
        filters method type from metadata
        :param metadata_obj: metadata of the dataset
        returns: method type as string
        """
        if metadata_obj["type"][0] == "Expression profiling by array":
            LOGGER.info("Dataset is microarray")
            return "microarray"
        elif metadata_obj["type"][0] == "Expression profiling by high throughput sequencing":
            LOGGER.info("Dataset is RNA-seq")
            return "rnaseq"
        else:
            LOGGER.warning(f"Unknown {metadata_obj['type'][0]}")
            return "unknown"

    def _clean_up_samples(self, file_list: list = None):
        """ 
        removes downloaded files 
        :param file_list: list of files to be removed
        """
        for file in file_list:
            file = file + ".txt"
            try:
                os.remove(file)
            except OSError:
                LOGGER.warning(f"Error removing file {file}")
                pass


