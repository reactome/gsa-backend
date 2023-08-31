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

        print(gse.metadata)

        #metadata_obj = ExternalData(id=identifier, title=gse.metadata["title"], type=gse.metadata["type"], description=gse.metadata["summary"])
        return None
    
    def _get_data_type(self, metadata_obj) -> str:
        if metadata_obj["type"] == "Expression profiling by array":   # make this more generic TODO
            return "microarray"
        elif metadata_obj["type"] == "Expression profiling by high throughput sequencing":
            return "rnaseq"
        else:
            return "unknown"
        

class MockMQ:
    def get_is_shutdown(self):
        return False
    def sleep(self, the_time):
        time.sleep(the_time)
from reactome_analysis_utils.models.dataset_request import DatasetRequestParameter

test_dataset = "GSE219121"
parameters = [DatasetRequestParameter("dataset_id", test_dataset)]
fetcher = GeoFetcher()
external_data = ExternalData()
external_data = fetcher.load_dataset(parameters, MockMQ)