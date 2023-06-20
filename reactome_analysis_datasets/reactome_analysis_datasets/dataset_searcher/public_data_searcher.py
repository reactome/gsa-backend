import os
from whoosh.fields import Schema, TEXT, KEYWORD, NUMERIC
from whoosh.index import create_in
from whoosh import index
from whoosh import qparser
from whoosh.query import Term
from itertools import groupby
from enum import Enum
from dataclasses import dataclass

from reactome_analysis_datasets.reactome_analysis_datasets.dataset_fetchers.grein_fetcher import GreinFetcher
from reactome_analysis_datasets.reactome_analysis_datasets.dataset_fetchers.expression_atlas_fetcher import ExpressionAtlasFetcher

schema = Schema(data_source=KEYWORD, id=TEXT(stored=True), title=TEXT(stored=True), species=KEYWORD,
                description=TEXT(stored=True), no_samples=NUMERIC(stored=True))

@dataclass
class Species:
    SPECIES_DICT = {
    }

class Generate_search_values():

    def setup_search_events(self):

        if not os.path.exists("../dataset_fetchers/index"):
            os.mkdir("../dataset_fetchers/index")

        ix = create_in("index", schema)

        writer = ix.writer()
        fetcher_grein = GreinFetcher()
        grein_datasets = fetcher_grein.get_available_datasets()
        for dataset in grein_datasets:
            writer.add_document(data_source="grein", id=str(dataset['geo_accession']), title=str(dataset['title']), species=str(dataset['species']),
                                description=str(dataset['study_summary']), no_samples=dataset['no_samples'])

        expression_atlas_fetcher = ExpressionAtlasFetcher()
        expression_atlas_datasets = expression_atlas_fetcher.get_available_datasets()
        for dataset in expression_atlas_datasets:
            writer.add_document(data_source="ebi_gxa", id=str(dataset['id']), title=str(dataset['title']),
                                species=str(dataset['species']),
                                description=str(''), no_samples=dataset['no_samples'])
        writer.commit()
        list_data = grein_datasets + expression_atlas_datasets
        species_in_datasets = self._get_species(list_data)
        print(species_in_datasets)
        Species.SPECIES_DICT = {item.replace(' ', '_'): item for item in species_in_datasets}

    def _get_species(self, datasets):
        values = []
        for dictionary in datasets:
            if 'species' in dictionary:
                values.append(dictionary['species'])
        species_values = [next(group) for key, group in groupby(sorted(values))]
        return species_values

class Searcher():

    dirname =""
    def __init__(self):
        self.dirname= "index"

    def index_search(self, search_fields: list, keyword: str, species: str) -> dict:
        ix = index.open_dir(self.dirname)
        schema = ix.schema

        search_field_description = "description"
        search_field_title = "title"
        search_field_species = "species"

        og = qparser.OrGroup.factory(0.9)
        mp = qparser.MultifieldParser(search_fields, schema, group=og)

        description_query = mp.parse(search_field_description + keyword)
        title_query = mp.parse(search_field_title + keyword)

        species_query = Term(search_field_species, species)
        combined_query = description_query | title_query | species_query

        # keyword_parser = mp.parse(keyword)

        with ix.searcher() as s:
            results = s.search(combined_query, terms=True, limit=100)
            result_dict = {}
            for result in results:
                result_dict[result['id']] = result.fields()
            return result_dict


keyword = "TFIIH kinase"
generator = Generate_search_values()
generator.setup_search_events()
search = Searcher()
result_dict = search.index_search(['description', 'title', 'species'], keyword, Species.SPECIES_DICT["Mus_musculus"])
print(result_dict)