import logging
import os
import pickle
import click

from whoosh.fields import Schema, TEXT, KEYWORD, NUMERIC
from whoosh.index import create_in
from whoosh import index
from whoosh import qparser
from whoosh.qparser import MultifieldParser

from reactome_analysis_api.searcher.overview_fetcher import PublicDataFetcher

LOGGER = logging.getLogger(__name__)


class PublicDatasetSearcher():
    """
    performs searching based on keyword and species, in previous created index
    the path for index creation is defined in the constructor!
    """

    _path = ""
    _species_list = None
    _ix = None
    _path_whitelist = ""

    schema = Schema(data_source=TEXT(stored=True), id=TEXT(stored=True), title=TEXT(stored=True),
                    species=TEXT(stored=True),
                    description=TEXT(stored=True), no_samples=NUMERIC(stored=True), technology=TEXT(stored=True),
                    resource_id=TEXT(stored=True), loading_parameters=TEXT(stored=True), link=TEXT(stored=True))

    def __init__(self, path: str, path_whitelist: str):
        """Initialize the public data searcher

        :param path: Path to where the index is stored
        :type path: str
        :param path_whitelist: Path to where the withelist is stored
        """
        self._path = path
        self._path_whitelist = path_whitelist

    def setup_search_events(self):
        """
        sets up the index for later search process based on schema, sets up species for later filtering. Data
        is stored in the path defined in the constructor.
        """
        LOGGER.info("Creating index for searching")
        if not os.path.exists(path=self._path):
            os.mkdir(self._path)

        ix = create_in(self._path, self.schema)

        LOGGER.debug("Created index: %s", self._path)
        writer = ix.writer()

        LOGGER.debug("Fetching available datasets")

        whitelist = self._get_whitelist_datasets(self._path_whitelist)
        datasets = PublicDataFetcher.get_available_datasets()  # all available public datasets

        for dataset in datasets:
            # ignore datasets without an id (happens sometimes in GREIN)
            if not "id" in dataset or type(dataset["id"]) != str or len(dataset["id"].strip()) < 3:
                continue
            # only valid datasets should be writen in the index
            # check if grein data is valid, assumes that ebi data is always valid
            if dataset['id'] in whitelist or dataset['resource_id'] == "ebi_gxa":
                writer.add_document(data_source=str(dataset['resource_id_str']), id=str(dataset['id']),
                                    title=str(dataset['title']),
                                    species=str(dataset['species']),
                                    description=str(dataset['study_summary']), no_samples=str(dataset['no_samples']),
                                    technology=str(dataset['technology']),
                                    resource_id=str(dataset['resource_id']),
                                    loading_parameters=str(dataset['loading_parameters']))
        writer.commit()

        # gets species based on public datasets
        species_in_datasets = self._get_species(datasets=datasets)
        with open(self._path + 'species.pickle', 'wb') as f:
            pickle.dump(species_in_datasets, f, pickle.HIGHEST_PROTOCOL)

    def _get_whitelist_datasets(self, path) -> list:
        """ 
        based on the whitelist datasets are fetched from the public datasets
        :param: path to whitelist 
        :return: list of datasets in whitelist used to build the index
        """
        whitelist = []
        list_reader = open(path, "r")
        for line in list_reader.readlines():
            whitelist.append(line.strip())
        list_reader.close()
        return whitelist

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
        if "character(0)" in species_values: species_values.remove("character(0)")
        return species_values

    def get_species(self) -> list:
        """
        returns species stored in binary file
        """
        # only try to load if the species list wasn't loaded before
        if self._species_list:
            return self._species_list

        # try to load the list from file
        try:
            if os.path.exists(self._path + 'species.pickle'):
                with open(self._path + 'species.pickle', 'rb') as f:
                    self._species_list = pickle.load(f)
        except FileNotFoundError:
            logging.error(f"File not found: {self._path}species.pickle")
            self._species_list = []

        return self._species_list

    def index_search(self, keyword: list, species: str = None, search_in_description: bool = False) -> list:
        """
        :param keyword, species: searches in title and description, species is based on the dictionary defined, searches only in
        species of the schema, search_in_description: boolean to switch of description searching
        :return dictionary of the search results
        """
        LOGGER.info("Searching keyword: %s, species: %s", keyword, species)

        if not self._ix:
            self._ix = index.open_dir(self._path)

        if species is None:
            species = "Homo sapiens"
            LOGGER.debug("Default species is set with: %s", species)

        # if there's only one keyword passed as string, add it as single item to a list
        if type(keyword) == str:
            keyword = [keyword]

        with self._ix.searcher() as searcher:

            if search_in_description == True:
                description_parser = MultifieldParser(["description", "title"], self.schema)
            else:
                description_parser = MultifieldParser(["title"], self.schema)
            species_parser = qparser.QueryParser("species", self.schema)

            query_string = " AND ".join(keyword)
            description_title_query = description_parser.parse(query_string)
            species_query = species_parser.parse(species)
            combined_query = description_title_query & species_query
            results = searcher.search(combined_query, limit=100)
            results_list = list()
            for result in results:
                if result["id"] != '':
                    results_list.append({
                        "id": result["id"],
                        "description": result["description"],
                        "title": result["title"],
                        "species": result["species"],
                        "resource_id": result["resource_id"],
                        "loading_parameters": result["loading_parameters"],
                        "data_source": result["data_source"],
                        "web_link": result["link"]
                    })
            return results_list


@click.command()
@click.option('--path', default=None,
              help="If set, the path to store the search index in. Otherwise the environment variable 'SEARCH_INDEX_PATH' is used.")
def create_search_index(path, path_whitelist):
    """Create the initial search index.
    """
    # set the logging
    logging.basicConfig(level=logging.DEBUG)

    # use the environment variable if no arg passed
    if not path:
        path = os.getenv("SEARCH_INDEX_PATH", "../")

    if not path_whitelist:
        path_whitelist = os.getenv("SEARCH_INDEX_WHITELIST", "../")

    # create the search
    searcher = PublicDatasetSearcher(path=path, path_whitelist=path_whitelist)
    searcher.setup_search_events()
