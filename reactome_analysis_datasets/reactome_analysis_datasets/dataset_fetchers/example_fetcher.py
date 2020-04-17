"""
Dataset fetcher the retrieves the stored example datasets

This class expects that the datasets are stored in EXAMPLE_DIRECTORY and
follow the pattern [dataset_id].data (text file) and [dataset_id].summary
(JSON encoded text file)
"""

import os
import json
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, DatasetFetcherException


class ExampleDatasetFetcher(DatasetFetcher):
    def load_dataset(self, parameters: list, reactome_mq) -> (str, ExternalData):
        """
        Loads the defined example dataset
        :param parameters: A list of DatasetRequestParameter objects.
        :param reactome_mq: Not used
        :returns: (data, summary)
        """
        # get the id parameter
        identifier = self._get_parameter("dataset_id", parameters)

        if not identifier:
            raise DatasetFetcherException("Missing required parameter 'dataset_id' to load the example dataset.")

        # prevent the injection of "mean" characters
        identifier = identifier.replace("/", "_")
        identifier = identifier.replace(".", "_")
        identifier = identifier.replace("$", "_")

        # build the path
        data_dir = os.getenv("EXAMPLE_DIRECTORY", "/data/examples")
        data_file = os.path.join(data_dir, "{}.data".format(identifier))
        summary_file = os.path.join(data_dir, "{}.summary".format(identifier))

        if not os.path.isfile(data_file) or not os.path.isfile(summary_file):
            raise DatasetFetcherException("Unknown example data identifier {}".format(identifier))

        # load the data
        self._update_status(progress=0.2, message="Loading expression data")

        try:
            example_data = self.load_textfile(data_file)
        except Exception:
            raise DatasetFetcherException("Failed to load data for {}".format(identifier))

        self._update_status(progress=0.7, message="Loading summary data")

        try:
            summary_data = self.load_textfile(summary_file)
        except Exception:
            raise DatasetFetcherException("Failed to load summary data for {}".format(identifier))

        # convert the summary into an object
        self._update_status(progress=0.8, message="Converting loaded data")
        try:
            summary_dict = json.loads(summary_data)
            summary_obj = ExternalData.from_dict(summary_dict)
        except Exception:
            raise DatasetFetcherException("Failed to load a valid summary for {}".format(identifier))

        # return the data
        return (example_data, summary_obj)

    def get_dataset_id(self, parameters: list) -> str:
        """
        Returns the dataset identifier
        :param parameters: A list of DatasetRequestParameter objects.
        :returns: The dataset identifier
        """
        return self._get_parameter("dataset_id", parameters)

    @staticmethod
    def load_textfile(filename: str) -> str:
        with open(filename, "r") as reader:
            text = reader.read()

        return text
