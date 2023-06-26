import logging
import os

from whoosh.fields import Schema, TEXT, KEYWORD, NUMERIC
from whoosh.index import create_in
from whoosh import index
from whoosh import qparser
from itertools import groupby
from dataclasses import dataclass

from reactome_analysis_datasets.reactome_analysis_datasets.dataset_fetchers.grein_fetcher import GreinFetcher
from reactome_analysis_datasets.reactome_analysis_datasets.dataset_fetchers.expression_atlas_fetcher import \
    ExpressionAtlasFetcher

schema = Schema(data_source=KEYWORD, id=TEXT(stored=True), title=TEXT(stored=True), species=TEXT(stored=True),
                description=TEXT(stored=True), no_samples=NUMERIC(stored=True))

LOGGER = logging.getLogger(__name__)

@dataclass
class Species:
    SPECIES_DICT = {
    }


class Generate_search_values():
    """
    Generator to build the search index based on the overview parameters
    """
    PATH= "../dataset_fetchers/index"

    def setup_search_events(self):
        """
        sets up the index for later search process based on schema, sets up species for later filtering
        """
        LOGGER.info("Creating index for searching")
        if not os.path.exists(path=self.PATH):
            os.mkdir(self.PATH)

        ix = create_in("index", schema)
        LOGGER.info("Created index: ", self.PATH)
        writer = ix.writer()
        fetcher_grein = GreinFetcher()
        LOGGER.info("Fetching available datasets from GREIN")
        grein_datasets = fetcher_grein.get_available_datasets()
        for dataset in grein_datasets:
            writer.add_document(data_source="grein", id=str(dataset['geo_accession']), title=str(dataset['title']),
                                species=str(dataset['species']),
                                description=str(dataset['study_summary']), no_samples=dataset['no_samples'])

        expression_atlas_fetcher = ExpressionAtlasFetcher()
        LOGGER.info("Fetching available datasets from ExpressionAtlas")
        expression_atlas_datasets = expression_atlas_fetcher.get_available_datasets()
        for dataset in expression_atlas_datasets:
            writer.add_document(data_source="ebi_gxa", id=str(dataset['id']), title=str(dataset['title']),
                                species=str(dataset['species']),
                                description=str(''), no_samples=dataset['no_samples'])
        writer.commit()
        list_data = grein_datasets + expression_atlas_datasets
        species_in_datasets = self._get_species(list_data)  # gets species based on public datasets
        Species.SPECIES_DICT = {item.replace(' ', '_'): item for item in species_in_datasets}

    def _get_species(self, datasets) -> list:
        """
        :param datasets: list of dictionaries from public datasets
        :return species_values: list of species in public datasets
        """
        values = []
        for dictionary in datasets:
            if 'species' in dictionary:
                values.append(dictionary['species'])
        species_values = [next(group) for key, group in groupby(sorted(values))]
        species_values.sort()
        return species_values


class Searcher():
    """
    performs searching based on keyword and species, in previous created index
    """

    dirname = ""

    def __init__(self):
        self.dirname = "index"

    def index_search(self, keyword: str, species: str) -> dict:
        """
        :param keyword, species: searches in title and description, species is based on the dictionary defined, searches only in
        species of the schema
        :return dictionary of the search results
        """
        LOGGER.info("Searching keyword: ", keyword, "species: ", keyword)
        ix = index.open_dir(self.dirname)

        with ix.searcher() as searcher:
            description_parser = qparser.QueryParser("description", schema)
            title_parser = qparser.QueryParser("title", schema)
            species_parser = qparser.QueryParser("species", schema)

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
                    "species": result["species"]
                }

            return result_dict


"""example script"""
generate_search_values = Generate_search_values()
generate_search_values.setup_search_events()
searcher = Searcher()
search_result = searcher.index_search("gene", "Homo_sapiens")
print(search_result)