import logging

from flask import abort, Response
from reactome_analysis_api.models.analysis_result import AnalysisResult  # noqa: E501
from reactome_analysis_utils.reactome_storage import ReactomeStorage, ReactomeStorageException

LOGGER = logging.getLogger(__name__)


def get_result(analysisId):  # noqa: E501
    """Retrieves the result for the completed analysis task

     # noqa: E501

    :param analysisId: The analysis identified returned by &#39;/analysis&#39;
    :type analysisId: str

    :rtype: AnalysisResult
    """
    try:
        # check if an extension was present
        extension = None

        if "." in analysisId:
            extension = analysisId[analysisId.find(".") + 1:]
            analysisId = analysisId[:analysisId.find(".")]

        storage = ReactomeStorage()

        if extension == "xlsx":
            xlsx_file = storage.get_result(analysis_identifier=analysisId, data_type="report")

            if xlsx_file is not None:
                return Response(response=xlsx_file, status=200, headers={"content-type": "application/xlsx"})
        elif extension == "pdf":
            pdf_file = storage.get_result(analysis_identifier=analysisId, data_type="pdf_report")

            if pdf_file is not None:
                return Response(response=pdf_file, status=200, headers={"content-type": "application/pdf"})
        elif extension == "r":
            r_file = storage.get_result(analysis_identifier=analysisId, data_type="r_script")

            if r_file is not None:
                return Response(response=r_file, status=200, headers={"content-type": "text/plain", 
                                                                      "content-disposition": "attachment; filename=\"ReactomeGSA_analysis_script.R\""})
        else:
            result = storage.get_result(analysisId)

            if result is not None:
                return Response(response=result, status=200, headers={"content-type": "application/json"})

        # find out why the result doesn't exist
        status = storage.get_status(analysisId)

        if not status:
            LOGGER.debug("Unknown identifier to get_result: " + analysisId)
            abort(404, "Unknown analysis identifier passed.")

        # the identifier is valid, so for some reason the result is not ready (yet)
        abort(406, "Analysis is not complete.")
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(503, "Failed to connect to storage system. Please try again in a few minutes.")


def get_status(analysisId):  # noqa: E501
    """Retrieves the status for the specified analysis.

     # noqa: E501

    :param analysisId: The analysis identifier returned by &#39;/analysis&#39;
    :type analysisId: str

    :rtype: InlineResponse200
    """
    try:
        storage = ReactomeStorage()

        status = storage.get_status(analysisId)

        if status is None:
            LOGGER.debug("Unknown identifier passed to get_status: " + analysisId)
            abort(404, "Unknown identifier")
        else:
            # return a Response object to prevent connexion from
            # de-serializing the object into a JSON object
            return Response(response=status, status=200, headers={"content-type": "application/json"})
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(503, "Failed to connect to storage system. Please try again in a few minutes.")


def get_report_status(analysisId):
    """Retrieves the status for the specified report.

    :param analysisId: The analysis identifier
    :type analysisId: str
    :rtype: InlineResponse200
    """
    try:
        storage = ReactomeStorage()

        status = storage.get_status(analysis_identifier=analysisId, data_type="report")

        if status is None:
            LOGGER.debug("Unknown identifier passed to get_report_status: " + analysisId)
            abort(404, "Unknown identifier")
        else:
            # return a Response object to prevent connexion from
            # de-serializing the object into a JSON object
            return Response(response=status, status=200, headers={"content-type": "application/json"})
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(503, "Failed to connect to storage system. Please try again in a few minutes.")
