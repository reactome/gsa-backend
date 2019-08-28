"""
Dataset fetcher the retrieves the stored example datasets

This class expects that the datasets are stored in EXAMPLE_DIRECTORY and
follow the pattern [dataset_id].data (text file) and [dataset_id].summary
(JSON encoded text file)
"""

import os
import json
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData


class ExampleDatasetFetcher(DatasetFetcher):
    def load_dataset(self, identifier: str) -> (str, ExternalData):
        """
        Loads the defined example dataset
        :params identifier: The dataset identifier
        :returns: (data, summary)
        """
        # prevent the injection of "mean" characters
        identifier = identifier.replace("/", "_")
        identifier = identifier.replace(".", "_")
        identifier = identifier.replace("$", "_")

        # build the path
        data_dir = os.getenv("EXAMPLE_DIRECTORY", "/data/examples")
        data_file = os.path.join(data_dir, "{}.data".format(identifier))
        summary_file = os.path.join(data_dir, "{}.summary".format(identifier))

        if not os.path.isfile(data_file) or not os.path.isfile(summary_file):
            raise Exception("Unknown example data identifier {}".format(identifier))

        # load the data
        try:
            example_data = self.load_textfile(data_file)
        except Exception:
            raise Exception("Failed to load data for {}".format(identifier))

        try:
            summary_data = self.load_textfile(summary_file)
        except Exception:
            raise Exception("Failed to load summary data for {}".format(identifier))

        # convert the summary into an object
        try:
            summary_dict = json.loads(summary_data)
            summary_obj = ExternalData.from_dict(summary_dict)
        except Exception:
            raise Exception("Failed to load a valid summary for {}".format(identifier))

        # return the data
        return (example_data, summary_obj)

    @staticmethod
    def load_textfile(filename: str) -> str:
        with open(filename, "r") as reader:
            text = reader.read()

        return text