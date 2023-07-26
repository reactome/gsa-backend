import logging
import os

from whoosh.fields import Schema, TEXT, KEYWORD, NUMERIC
from whoosh.index import create_in
from whoosh import index
from whoosh import qparser
import pickle
from dataclasses import dataclass

from reactome_analysis_api.reactome_analysis_api.searcher.overview_fetcher import Fetcher

LOGGER = logging.getLogger(__name__)


@dataclass
class Species:
    SPECIES_DICT = {
    }


class PublicDatasetSearcher():
    """
    performs searching based on keyword and species, in previous created index
    the path for index creation is defined in the constructor!
    """

    _path = ""
    schema = Schema(data_source=KEYWORD, id=TEXT(stored=True), title=TEXT(stored=True), species=TEXT(stored=True),
                    description=TEXT(stored=True), no_samples=NUMERIC(stored=True), technology=TEXT(stored=True),
                    resource_id=TEXT(stored=True), loading_parameters=TEXT(stored=True))

    def __init__(self, path):
        self._path = path

    def setup_search_events(self):  # top level funktion hinzufpügen create_search_index
        """
        sets up the index for later search process based on schema, sets up species for later filtering
        """
        LOGGER.info("Creating index for searching")
        if not os.path.exists(path=self._path):
            os.mkdir(self._path)

        ix = create_in(self._path, self.schema)
        LOGGER.info("Created index: ", self._path)
        writer = ix.writer()
        LOGGER.info("Fetching available datasets from GREIN")
        grein_datasets = Fetcher.get_available_datasets_grein()
        for dataset in grein_datasets:
            writer.add_document(data_source="grein", id=str(dataset['id']), title=str(dataset['title']),
                                species=str(dataset['species']),
                                description=str(dataset['study_summary']), no_samples=str(dataset['no_samples']),
                                technology=str(dataset['technology']),
                                resource_id=str(dataset['resource_id']),
                                loading_parameters=str(dataset['loading_parameters']))

        LOGGER.info("Fetching available datasets from ExpressionAtlas")
        expression_atlas_datasets = Fetcher.get_available_datasets_expression_atlas()
        for dataset in expression_atlas_datasets:
            writer.add_document(data_source="ebi_gxa", id=str(dataset['id']), title=str(dataset['title']),
                                species=str(dataset['species']),
                                description=str(''), no_samples=dataset['no_samples'],
                                technology=str(dataset['technology']),
                                resource_id=str(dataset['resource_id']),
                                loading_parameters=str(dataset['loading_parameters']))
        writer.commit()
        list_data = grein_datasets + expression_atlas_datasets
        species_in_datasets = self._get_species(list_data)  # gets species based on public datasets
        Species.SPECIES_DICT = {item.replace(' ', '_'): item for item in species_in_datasets}
        with open('species.pickle', 'wb') as f:
            pickle.dump(Species.SPECIES_DICT, f, pickle.HIGHEST_PROTOCOL)

    def _get_species(self, datasets) -> set:
        """
        :param datasets: list of dictionaries from public datasets
        :return species_values: list of species in public datasets
        """
        values = set()
        for dictionary in datasets:
            if 'species' in dictionary:
                values.add(dictionary['species'])
        species_values = sorted(values)
        return species_values

    def get_species(self) -> set:
        """
        returns species stored in binary file
        """
        with open('species.pickle', 'rb') as f:
            return pickle.load(f)

    def index_search(self, keyword: str, species: str) -> dict:
        """
        :param keyword, species: searches in title and description, species is based on the dictionary defined, searches only in
        species of the schema
        :return dictionary of the search results
        """
        LOGGER.info("Searching keyword: ", keyword, "species: ", species)
        ix = index.open_dir(self._path)

        with ix.searcher() as searcher:
            description_parser = qparser.QueryParser("description", self.schema)
            title_parser = qparser.QueryParser("title", self.schema)
            species_parser = qparser.QueryParser("species", self.schema)

            description_query = description_parser.parse(keyword)
            title_query = title_parser.parse(keyword)
            species_query = species_parser.parse(species)
            combined_query = (description_query | title_query) & species_query
            results = searcher.search(combined_query, limit=100)
            result_dict = {}
            for result in results:
                result_dict[result["id"]] = {
                    "description": result["description"],
                    "title": result["title"],
                    "species": result["species"],
                    "resource_id": result["resource_id"],  # change this to grein or expresion atlas
                    "loading_parameters": result["loading_parameters"]
                }

            return result_dict
