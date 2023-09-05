import logging
import GEOparse as geoparser
from typing import Tuple
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, DatasetFetcherException
from reactome_analysis_api.models.external_data_sample_metadata import ExternalDataSampleMetadata   

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
        :returns empty string insteda of count matrix, metadata of all datasets
        """
        identifier = self._get_parameter(name="dataset_id", parameters=parameters).strip()
        if not identifier:
            raise DatasetFetcherException(
                "Missing required parameter 'dataset_id' to load the example dataset.")
        if not identifier[0:3] == "GSE":
            raise DatasetFetcherException(
                f"{identifier} is not a valid geo_accession for GREIN")

        # load the data
        LOGGER.info(f"Loading dataset {identifier} from GEO")  
        try:
            gse = geoparser.get_GEO(identifier)   # fetching data from Geo via geo parser
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
                                        sample_ids= gse.metadata["sample_id"], 
                                        sample_metadata=sample_metadata_list,   
                                        default_parameters=None)
        
        return ("", metadata_obj)
    
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
    
