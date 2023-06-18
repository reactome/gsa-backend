from whoosh.fields import Schema, TEXT, KEYWORD, NUMERIC
from grein_fetcher import GreinFetcher
from whoosh.index import create_in
import os
from whoosh import index
from whoosh import qparser


schema = Schema(data_source=KEYWORD, id=TEXT(stored=True), title=TEXT(stored=True), species=KEYWORD,
                description=TEXT(stored=True))


class Generate_search_values():

    @staticmethod
    def setup_search_events():

        if not os.path.exists("index"):
            os.mkdir("index")

        ix = create_in("index", schema)

        writer = ix.writer()
        fetcher_grein = GreinFetcher()
        grein_datasets = fetcher_grein.get_available_datasets(10)  # fetches

        for dataset in grein_datasets:
            writer.add_document(data_source="GREIN", id=str(dataset['geo_accession']), title=str(dataset['title']), species=str(dataset['species']),
                                description=str(dataset['study_summary']))
        writer.commit()


class Searcher():

    dirname =""
    def __init__(self):
        self.dirname= "index"

    def index_search(self, search_fields, keyword):
        ix = index.open_dir(self.dirname)
        schema = ix.schema

        og = qparser.OrGroup.factory(0.9)
        mp = qparser.MultifieldParser(search_fields, schema, group=og)

        keyword_parser = mp.parse(keyword)

        with ix.searcher() as s:
            results = s.search(keyword_parser, terms=True, limit=100)
            result_dict = {}
            for result in results:
                result_dict[result['id']] = result.fields()
            return result_dict


Generate_search_values.setup_search_events()
search = Searcher()
result_dict = search.index_search (['description', 'title'], "TFIIH kinase")