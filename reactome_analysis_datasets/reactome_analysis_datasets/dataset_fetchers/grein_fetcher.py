import grein_loader
import pandas as pd
import abstract_dataset_fetcher


class GreinFetcher(abstract_dataset_fetcher.DatasetFetcher):

    def load_dataset(self, parameters: list, reactome_mq: abstract_dataset_fetcher.ReactomeMQ) -> (
            str, abstract_dataset_fetcher.ExternalData):
        description, metadata, count_matrix = grein_loader.load_dataset(parameters[0])
        metadata = abstract_dataset_fetcher.ExternalData(metadata)
        count_matrix_tsv = count_matrix.to_csv(sep="\t")
        return metadata, count_matrix_tsv

    def load_overview(self, no_datasets: int):
        return grein_loader.load_overview(no_datasets)
