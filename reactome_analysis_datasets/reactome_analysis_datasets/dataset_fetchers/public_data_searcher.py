import os
from whoosh.fields import Schema, TEXT, KEYWORD, NUMERIC
from whoosh.index import create_in
from whoosh import index
from whoosh import qparser

from grein_fetcher import GreinFetcher
from expression_atlas_fetcher import ExpressionAtlasFetcher

schema = Schema(data_source=KEYWORD, id=TEXT(stored=True), title=TEXT(stored=True), species=KEYWORD,
                description=TEXT(stored=True), no_samples=NUMERIC(stored=True))

class Generate_search_values():

    @staticmethod
    def setup_search_events():

        if not os.path.exists("index"):
            os.mkdir("index")

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


class Searcher():

    dirname =""
    def __init__(self):
        self.dirname= "index"

    def index_search(self, search_fields: list, keyword: str) -> dict:
        ix = index.open_dir(self.dirname)
        schema = ix.schema

        og = qparser.OrGroup.factory(0.9)
        mp = qparser.MultifieldParser(search_fields, schema)

        keyword_parser = mp.parse(keyword)

        with ix.searcher() as s:
            results = s.search(keyword_parser, terms=False, limit=100)
            result_dict = {}
            for result in results:
                result_dict[result['id']] = result.fields()
            return result_dict
