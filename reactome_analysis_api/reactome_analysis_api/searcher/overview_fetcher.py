import grein_loader
import requests
import json
import logging

LOGGER = logging.getLogger(__name__)


class fetcher():

    @staticmethod
    def get_available_datasets_grein(no_datasets: int = None) -> list:
        """
        Returns the available datasets, is used exclusively during the index build for the keyword searcher
        :param no_datasets: number of datasets to retrieve
        :returns: datasets in ExternalData format
        """
        grein_datasets = grein_loader.load_overview(no_datasets)
        list_overview = []
        for dataset in grein_datasets:
            overview_dict = {
                "id": dataset["geo_accession"],
                "title": dataset["title"],
                "study_summary": dataset["study_summary"],
                "species": dataset["species"],
                "no_samples": dataset["no_samples"],
                "technology": "",
                "resource_id": dataset["geo_accession"],
                "loading_parameters": json.dumps({"id": dataset["geo_accession"]})
            }
            list_overview.append(overview_dict)
        return list_overview


    @staticmethod
    def get_available_datasets_expression_atlas(no_datasets: int = None) -> list:
        """
        Returns the available datasets, is used exclusively during the index build for the keyword searcher
        :param no_datasets: number of datasets to retrieve
        :returns: datasets in ExternalData format
        """

        if no_datasets is None:
            no_datasets = 10000

        experiments_external_data_list = list()
        try:
            experiments_url = "https://www.ebi.ac.uk/gxa/json/experiments"
            response = requests.get(experiments_url)
            response.raise_for_status()
            json_response = response.json()
            experiments_list = json_response['experiments'][0:no_datasets]

            for experiment in experiments_list:
                experiment_data_dict = {
                    "id": experiment['experimentAccession'],
                    "title": experiment['experimentDescription'],
                    "study_summary": "",
                    "species": experiment['species'],
                    "no_samples": experiment['numberOfAssays'],
                    "technology": experiment['technologyType'],
                    "resource_id": experiment['experimentAccession'],
                    "loading_parameters": json.dumps({"id": experiment['experimentAccession']})
                }
                experiments_external_data_list.append(experiment_data_dict)
        except requests.exceptions.RequestException:
            LOGGER.error("Response not available")
        return experiments_external_data_list