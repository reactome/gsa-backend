"""
The main entry point for the reactome_analysis_report package.
"""

import logging
import os

from prometheus_client import start_http_server

LOG_FORMAT = ('%(levelname) -10s %(asctime)s %(name) -30s %(funcName) '
              '-35s %(lineno) -5d: %(message)s')

LOGGER = logging.getLogger(__name__)

# start the prometheus client http server
start_http_server(port=int(os.getenv("PROMETHEUS_PORT", 9000)))


def main():
    logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)

    # set config for other packages
    pika_logger = logging.getLogger("pika")
    pika_logger.setLevel(level=logging.ERROR)

    report_generator = reactome_analysis_report_generator.ReactomeAnalysisReportGenerator()

    try:
        report_generator.start_listening()
    except Exception as e:
        LOGGER.error("Report process failed: " + str(e))
        # stack trace on debug
        LOGGER.debug("Error:", exc_info=1)
        report_generator.shutdown()

    LOGGER.info("Exiting.")


if __name__ == "__main__":
    main()
