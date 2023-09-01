import logging
import GEOparse as geoparser
from typing import Tuple
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, DatasetFetcherException
from reactome_analysis_api.models.external_data_sample_metadata import ExternalDataSampleMetadata   


import time

LOGGER = logging.getLogger(__name__)

class GeoFetcher(DatasetFetcher):

    def __init__(self):
        super().__init__()

    def load_dataset(self, parameters: list, reactome_mq) -> Tuple[str, ExternalData]:
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
            gse = geoparser.get_GEO(identifier, destdir=".")   # fetching data from Geo via geo parser
        except Exception as e:
            raise DatasetFetcherException(f"Error loading dataset {identifier} from GEO: {e}")

        # create metadata
        experiment_type = self._get_data_type(gse.metadata)
        sample_metadata_list = self._create_sample_metadata(gse.metadata["sample_id"])
        metadata_obj = ExternalData(id=identifier, 
                                    title=gse.metadata["title"],
                                    type=experiment_type, 
                                    description=gse.metadata["summary"],
                                    group=None,
                                    sample_ids= gse.metadata["sample_id"], 
                                    sample_metadata=sample_metadata_list,       #TODO  
                                    default_parameters=None)
        return ("", metadata_obj)
    
    def _create_sample_metadata(self, sample_list) -> list[ExternalDataSampleMetadata]:    
        list_of_sample_metadata = []
        for sample in sample_list:
            sample_metadata = geoparser.get_GEO(sample, destdir=".")   
            LOGGER.info(f"Loading sample metadata {sample} from GEO")
            #metadata_list = [str(sample_metadata.metadata["title"][0])] #TODO change to tuple
            metadata_list = [str(sample_metadata.metadata)]
            external_metadata_dict = {
                "name": str(sample),
                "values": metadata_list
                }
            #print(external_metadata_dict)
            external_metadata_object = ExternalDataSampleMetadata.from_dict(external_metadata_dict)
            #print(external_metadata_object)
            list_of_sample_metadata.append(external_metadata_object)
        return list_of_sample_metadata            


    def _get_data_type(self, metadata_obj) -> str:
        #print(metadata_obj["type"])
        if metadata_obj["type"][0] == "Expression profiling by array":  
            return "microarray"
        elif metadata_obj["type"][0] == "Expression profiling by high throughput sequencing":
            return "rnaseq"
        else:
            LOGGER.warning(f"Unknown {metadata_obj['type'][0]}")
            return "unknown"
        