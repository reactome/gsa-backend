import logging
import socket
import uuid
import zlib
import json

import connexion
import prometheus_client
from flask import abort
from reactome_analysis_api import methods, input_deserializer
from reactome_analysis_api.encoder import JSONEncoder
from reactome_analysis_api.models.analysis_status import AnalysisStatus
from reactome_analysis_api.models.data_type import DataType
from reactome_analysis_utils.reactome_mq import ReactomeMQ, ReactomeMQException
from reactome_analysis_utils.reactome_storage import ReactomeStorage, ReactomeStorageException
from reactome_analysis_utils.models.analysis_request import AnalysisRequest

LOGGER = logging.getLogger(__name__)
STARTED_ANALYSIS_COUNTER = prometheus_client.Counter('reactome_api_started_analyses',
                                                     'Analysis requests started through the API', ["client"])

MISSING_DATA_TOKEN_COUNTER = prometheus_client.Counter('reactome_api_missing_token',
                                                       'Missing data for data tokens.')


def list_methods():  # noqa: E501
    """Lists the available analysis methods

    Some analysis services may provide different methods to process the data. The available methods and their specification can be retrieved using this function. Most services will only support a single method though. # noqa: E501


    :rtype: List[Method]
    """
    return methods.get_available_methods()


def list_types():
    """
    Lists the supported data types
    :return: List[DataType]
    """
    data_types = [
        DataType(id="rnaseq_counts",
                 name="RNA-seq (raw counts)",
                 description="Raw RNA-seq based read counts per gene (recommended)."),
        DataType(id="rnaseq_norm",
                 name="RNA-seq (normalized)",
                 description="log2 transformed, normalized RNA-seq based read counts per gene (f.e. RPKM, TPM)"),
        DataType(id="proteomics_int",
                 name="Proteomics (intensity)",
                 description="Intensity-based quantitative proteomics data (for example, "
                             "iTRAQ/TMT or intensity-based label-free quantitation). Values "
                             "must be log2 transformed."),
        DataType(id="proteomics_sc",
                 name="Proteomics (spectral counts)",
                 description="Raw spectral-counts of label-free proteomics experiments"),
        DataType(id="microarray_norm",
                 name="Microarray (normalized)",
                 description="Normalized and log2 transformed microarray-based gene expression values."),
        DataType(id="ribo_seq",
                 name="Ribo-seq",
                 description="Ribo seq data analysis using RNA-seq data and Ribo-seq data")
    ]

    return data_types


def start_analysis(body):  # noqa: E501
    """Performs the specified gene set analysis

     # noqa: E501

    :param body: Specification of analysis to perform
    :type body: dict | bytes

    :rtype: str
    """
    user_client = "Unknown"

    # set the basic user client
    if "r-curl" in connexion.request.headers.get("User-Agent", "Unknown"):
        user_client = "ReactomeGSA R"
    elif "reactome.org/PathwayBrowser" in connexion.request.headers.get("Referer", "Unknown"):
        user_client = "PathwayBrowser"

    # get the JSON-encoded dict from the request object
    if connexion.request.is_json:
        analysis_dict = connexion.request.get_json(cache=False)
    # de-compress if it's a gzipped string
    elif connexion.request.content_type == "application/gzip":
        LOGGER.debug("Received gzipped analysis request. Decompressing...")

        decompressed_string = zlib.decompress(connexion.request.data)
        analysis_dict = json.loads(decompressed_string)
        
        # free the memory again
        del decompressed_string
    else:
        LOGGER.debug("Invalid analysis request submitted. Request body does not describe a JSON object.")
        abort(406, "Invalid analysis request submitted. Request body does not describe a JSON object.")
        return
    
    try:
        analysis_request = input_deserializer.create_analysis_input_object(analysis_dict)
    except Exception as e:
        if "Unknown analysis method" in str(e):
            LOGGER.debug("Unknown analysis method submitted: " + analysis_dict["methodName"])
            abort(404, "Unknown analysis method selected.")

        LOGGER.debug("Invalid analysis request: " + str(e))

        # try to find a nice error message
        if "Invalid value for `type`" in str(e):
            abort(400, "Invalid request: " + str(e))
        
        abort(400, "Invalid analysis request submitted")

    # make sure all datasets have unique names
    all_names = [dataset.name for dataset in analysis_request.datasets]

    if len(all_names) != len(set(all_names)):
        LOGGER.debug("Analysis request contains duplicate names")
        abort(406, "Datasets must not have duplicate names")

    # make sure the analysis design is present
    for n_dataset in range(0, len(analysis_request.datasets)):
        if not analysis_request.datasets[n_dataset].design:
            LOGGER.debug("Analysis request misses design")
            abort(406, "Invalid request. Dataset '{name}' misses the required experimental design.".format(name=analysis_request.datasets[n_dataset].name))
        if not analysis_request.datasets[n_dataset].design.comparison:
            LOGGER.debug("Analysis request misses design comparison")
            abort(406, "Invalid request. Dataset '{name}' misses the required comparison specification.".format(name=analysis_request.datasets[n_dataset].name))


    # generate an analysis id
    analysis_id = str(uuid.uuid1())

    try:
        storage = ReactomeStorage()

        # a very basic sanity check to make sure it's unique
        while storage.analysis_exists(analysis_id):
            analysis_id = str(uuid.uuid1())

        # Load request data from storage
        for n_dataset in range(0, len(analysis_dict["datasets"])):
            data = analysis_dict["datasets"][n_dataset]["data"]

            # check if dataset is riboseq data and process metadata accordingly
            if analysis_dict["datasets"][n_dataset]["type"] == "ribo_seq":
                analysis_dict["datasets"][n_dataset]["design"]["analysisGroup"] = analysis_dict["datasets"][n_dataset]["design"]["analysisGroup"] + analysis_dict["datasets"][n_dataset]["design"]["analysisGroup"]
                analysis_dict["datasets"][n_dataset]["design"]["samples"] = analysis_dict["datasets"][n_dataset]["design"]["samples"] + analysis_dict["datasets"][n_dataset]["design"]["samples"]
                
                if len(analysis_dict["datasets"][n_dataset]["design"]["samples"]) != analysis_dict["datasets"][n_dataset]["design"]["analysisGroup"]:
                    abort(500, "Samples and analysis groups do not match.")

            # Update for external datasets
            if data[0:4] == "rqu_" or len(data) < 20:
                # make sure the request data exists
                if not storage.request_token_exists(data):
                    MISSING_DATA_TOKEN_COUNTER.inc()
                    abort(500, "No data available for storage token '{}'".format(data))

                # load the data
                stored_data = storage.get_request_data(data)

                # if a bytes object is returned, this still has to be decoded
                if type(stored_data) == bytes:
                    stored_data = stored_data.decode("UTF-8")

                # update the request object
                analysis_dict["datasets"][n_dataset]["data"] = stored_data

        # Set the initial status
        encoder = JSONEncoder()

        status = AnalysisStatus(id=analysis_id, status="running", completed=0, description="Queued")
        storage.set_status(analysis_id, encoder.encode(status))

        # Save the request data
        analysis_dict["analysisId"] = analysis_id
        storage.set_analysis_request_data(token=analysis_id, data=encoder.encode(analysis_dict))

        try:
            # Submit the request to the queue
            queue = ReactomeMQ()
            queue.post_analysis(AnalysisRequest(request_id=analysis_id).to_json(), analysis_request.method_name)
            LOGGER.debug("Analysis " + analysis_id + " submitted to queue")
            queue.close()

            STARTED_ANALYSIS_COUNTER.labels(client=user_client).inc()

            return analysis_id
        except socket.gaierror as e:
            # update the status
            LOGGER.error("Failed to connect to queuing system: " + str(e))
            status = AnalysisStatus(id=analysis_id, status="failed", completed=0,
                                    description="Failed to connect to queuing system.")
            storage.set_status(analysis_id, encoder.encode(status))

            abort(503, "Failed to connect to queuing system. Please try again in a few seconds.")
        except ReactomeMQException as e:
            LOGGER.error("Failed to post message to queuing system: " + str(e))
            # update the status
            status = AnalysisStatus(id=analysis_id, status="failed", completed=0,
                                    description="Failed to connect to queuing system.")
            storage.set_status(analysis_id, encoder.encode(status))

            abort(503, "The number of analysis requests is currently too high. Please try again in a few minutes.")
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(503, "Failed to connect to storage system. Please try again in a few minutes.")
    except (socket.timeout, socket.gaierror) as e:
        LOGGER.error("Socket timeout connecting to storage or queuing system: " + str(e))
        abort(503, "Failed to connect to downstream system. Please try again in a few minutes.")
