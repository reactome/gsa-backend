"""
The main entry point for the reactome_analysis_report package.
"""

import logging
import os

from prometheus_client import start_http_server
from reactome_analysis_datasets.reactome_analysis_dataset_fetcher import ReactomeAnalysisDatasetFetcher
from reactome_analysis_utils.reactome_logging import get_default_logging_handlers


LOG_FORMAT = ('%(levelname) -10s %(asctime)s %(name) -30s %(funcName) '
              '-35s %(lineno) -5d: %(message)s')

LOGGER = logging.getLogger(__name__)

# start the prometheus client http server
start_http_server(port=int(os.getenv("PROMETHEUS_PORT", 9000)))


def main():
    logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT, handlers=get_default_logging_handlers())

    # set config for other packages
    pika_logger = logging.getLogger("pika")
    pika_logger.setLevel(level=logging.ERROR)
    logging.getLogger("rediscluster").setLevel(logging.CRITICAL)

    dataset_worker = ReactomeAnalysisDatasetFetcher()

    """try:  # Approach of setting up the search index 
        Generate_search_values.setup_search_events()
    except Exception as e:
        LOGGER.error("Dataset for search failed: "+str(e))"""

    try:
        dataset_worker.start_listening()
    except Exception as e:
        LOGGER.error("Datasets process failed: " + str(e))
        # stack trace on debug
        LOGGER.debug("Error:", exc_info=1)
        dataset_worker.shutdown()

    LOGGER.info("Exiting.")


if __name__ == "__main__":
    main()
