import ix as ix
from whoosh.fields import Schema, TEXT, KEYWORD
from whoosh.index import create_in
import os
from grein_fetcher import GreinFetcher
from expression_atlas_fetcher import ExpressionAtlasFetcher



class Generate_search_values():

    @staticmethod
    def setup_search_events():
        schema = Schema(load_dataset_information=KEYWORD, title=TEXT, species=KEYWORD, description=TEXT)

        if not os.path.exists("index"):
            os.mkdir("index")
        ix = create_in("index", schema)

        writer = ix.writer()
        fetcher_grein = GreinFetcher()
        grein_datasets = fetcher_grein.get_available_datasets()

        for dataset in grein_datasets:
            writer.add_document(load_dataset_information="GREIN", title=['title'], species=dataset['species'],
                                description=['study_summary'])

        """
        fetcher_expression_atlas = ExpressionAtlasFetcher()  # not merged yet
        fetcher_expression_atlas.get_available_datasets()
        
        for dataset in grein_datasets:
            writer.add_document(load_dataset_information="Expression Atlas", title=['title'], species=dataset['species'],
                                description=['study_summary'])
        """
        writer.commit()