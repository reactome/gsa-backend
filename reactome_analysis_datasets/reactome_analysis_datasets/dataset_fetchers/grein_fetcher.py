import grein_loader
import pandas as pd
import abstract_dataset_fetcher
import logging
import os
LOGGER = logging.getLogger(__name__)


class GreinFetcher(abstract_dataset_fetcher.DatasetFetcher):

    def load_dataset(self, parameters: list, reactome_mq: abstract_dataset_fetcher.ReactomeMQ) -> (
            str, abstract_dataset_fetcher.ExternalData):

        identifier = self._get_parameter("dataset_id", parameters)

        if not identifier:
            raise abstract_dataset_fetcher.DatasetFetcherException("Missing required parameter 'dataset_id' to load the example dataset.")

        # prevent the injection of "mean" characters
        identifier = identifier.replace("/", "_")
        identifier = identifier.replace(".", "_")
        identifier = identifier.replace("$", "_")

        # build the path
        data_dir = os.getenv("GREIN", "/data/grein_data")

        # load the data
        self._update_status(progress=0.2, message="Requesting data from GREIN")

        try:
            description, metadata, count_matrix = grein_loader.load_dataset(identifier)
        except Exception:
            raise abstract_dataset_fetcher.DatasetFetcherException("Failed to load data for {}".format(identifier))

        self._update_status(progress=0.7, message="Converting metadata")

        try:
            metadata_obj = abstract_dataset_fetcher.ExternalData.from_dict(metadata)
        except Exception:
            raise abstract_dataset_fetcher.DatasetFetcherException("Failed to load a valid summary for {}".format(identifier))

        self._update_status(progress=0.8, message="Converting countmatrix")
        try:
            count_matrix_tsv = count_matrix.to_csv(data_dir, sep="\t")
        except Exception:
            raise abstract_dataset_fetcher.DatasetFetcherException("Failed to load a valid summary for {}".format(identifier))


        #return data
        return (metadata_obj, count_matrix_tsv)

    def load_overview(self, no_datasets: int):
        return grein_loader.load_overview(no_datasets)
