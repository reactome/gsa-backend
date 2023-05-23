import connexion
import six
import logging
import prometheus_client
from flask import abort, Response
import uuid
import socket

from reactome_analysis_api.encoder import JSONEncoder
from reactome_analysis_api.models.dataset_loading_status import DatasetLoadingStatus  # noqa: E501
from reactome_analysis_api.models.external_data import ExternalData  # noqa: E501
from reactome_analysis_api import util
from reactome_analysis_utils.reactome_mq import ReactomeMQ, ReactomeMQException, DATASET_QUEUE
from reactome_analysis_utils.reactome_storage import ReactomeStorage, ReactomeStorageException
from reactome_analysis_utils.models.dataset_request import DatasetRequest, DatasetRequestParameter
from reactome_analysis_api.models.external_datasource import ExternalDatasource, ExternalDatasourceParameters


LOGGER = logging.getLogger(__name__)
DATASET_LOADING_COUNTER = prometheus_client.Counter("reactome_api_loading_datasets",
                                                    "External datasets loaded.")


def get_examples():  # noqa: E501
    """Lists the available example datasets

     # noqa: E501#


    :rtype: ExternalData
    """
    return [
        ExternalData(id="EXAMPLE_MEL_RNA", title="Melanoma RNA-seq example", type="rnaseq_counts",
                     description="RNA-seq analysis of melanoma associated B cells.",
                     group="GRISS_MELANOMA"),
        ExternalData(id="EXAMPLE_MEL_PROT", title="Melanoma proteomics example", type="proteomics_int",
                     description="Quantitative (TMT-labelled) proteomics analysis of melanoma associated B cells.",
                     group="GRISS_MELANOMA"),
        ExternalData(id="EXAMPLE_SC_B_CELLS", title="B cell scRNAseq example", type="rnaseq_counts",
                     description="Single-cell RNA-seq data of B cells extracted from the Jerby-Arnon at al. study (Cell 2018).",
                     group="SC_EXAMPLES"),
        ExternalData(id="GREIN", title="Public data from GREIN", type="rnaseq_counts",
                     description="Public dataset from Grein", group="")
    ]


def get_data_sources():  # noqa: E501
    """Lists the available external data sources

     # noqa: E501


    :rtype: ExternalDatasource
    """
    return [
        ExternalDatasource(id = "example_datasets", name = "Example datasets", 
                           description="Example datasets to quickly test the application.",
                           parameters=[
            ExternalDatasourceParameters(name="dataset_id", display_name="Dataset Id", 
                                         type="string", description="Identifier of the dataset", required=True)
        ]),
        ExternalDatasource(id = "ebi_gxa", name = "Expression Atlas", 
                           description="EBI's Expression Atlas resource for consistently reprocessed 'omics data.",
                           parameters=[
            ExternalDatasourceParameters(name="dataset_id", type="string", display_name="Dataset Id",
                                         description="Identifier of the dataset", required=True)
        ]),
        ExternalDatasource(id="ebi_sc_gxa", name = "Single Cell Expression Atlas",
                           description="EBI's Single Cell Expression Atlas resource for consistently reprocessed scRNA-seq data.",
                           parameters=[
            ExternalDatasourceParameters(name="dataset_id", display_name="Dataset Id",
                                         type="string", description="Identifier of the dataset", required=True),
            ExternalDatasourceParameters(name="k", type="int", display_name="K",
                                         description="Parameter k used to create the cell clusters", required=True),
        ]),
        ExternalDatasource(id="grein_data", name="GREIN Data",
                           description="Public data from GREIN",
                           parameters=[
            ExternalDatasourceParameters(name="dataset_id", display_name="Dataset Id",
                                         type="string", description="Identifier of the dataset", required=True)])
    ]


def get_data_loading_status(loadingId):  # noqa: E501
    """Retrieves the status for the dataset loading process.

     # noqa: E501

    :param loadingId: The loading identifier returned by &#39;/data/load&#39;
    :type loadingId: str

    :rtype: DatasetLoadingStatus
    """
    try:
        storage = ReactomeStorage()

        status = storage.get_status(analysis_identifier=loadingId, data_type="dataset")

        if status is None:
            LOGGER.debug("Unknown identifier passed to get_status: " + loadingId)
            abort(404, "Unknown identifier")
        else:
            # return a Response object to prevent connexion from
            # de-serializing the object into a JSON object
            return Response(response=status, status=200, headers={"content-type": "application/json"})
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(503, "Failed to connect to storage system. Please try again in a few minutes.")


def get_summary(datasetId):  # noqa: E501
    """Retrieves a summary of the loaded data. This function is only available once. The data is fully loaded.

     # noqa: E501

    :param datasetId: The dataset identifier used to trigger the download
    :type loadingId: str

    :rtype: ExternalData
    """
    try:
        storage = ReactomeStorage()

        if not storage.request_data_summary_exists(datasetId):
            abort(404, "Unknown identifier passed.")

        summary_data = storage.get_request_data_summary(datasetId)

        if summary_data is not None:
            return Response(response=summary_data, status=200, headers={"content-type": "application/json"})

        abort(404, "Unknown identifier passed.")
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(503, "Failed to connect to storage system. Please try again in a few minutes.")


def load_data(resourceId, parameters):  # noqa: E501
    """Start the retrieval of an external or example dataset.

     # noqa: E501

    :param resourceId: The identifier of the data source to load from
    :type resourceId: str

    :param parameters: The parameters for the selected resource.

    :rtype: str
    """
    try:
        storage = ReactomeStorage()

        # generate an id for the request
        loading_id = str(uuid.uuid1())

        # Set the initial status
        encoder = JSONEncoder()
        status = DatasetLoadingStatus(id=loading_id, status="running", completed=0, description="Queued")
        storage.set_status(loading_id, encoder.encode(status), data_type="dataset")

        # convert the parameters
        request_parameters = list()


        for dict_param in parameters:
            request_parameters.append(DatasetRequestParameter(name=dict_param["name"], value=dict_param["value"]))

        # create the request
        request = DatasetRequest(loading_id=loading_id, resource_id=resourceId, parameters=request_parameters)

        try:
            queue = ReactomeMQ(queue_name=DATASET_QUEUE)
            queue.post_analysis(analysis=request.to_json(), method="DatasetLoading")
            LOGGER.debug("Dataset process " + loading_id + " submitted to queue")
            queue.close()

            DATASET_LOADING_COUNTER.inc()

            return loading_id
        except socket.gaierror as e:
            # update the status
            LOGGER.error("Failed to connect to queuing system: " + str(e))
            status = DatasetLoadingStatus(id=loading_id, status="failed", completed=0,
                                    description="Failed to connect to queuing system.")
            storage.set_status(loading_id, encoder.encode(status), data_type="dataset")

            abort(503, "Failed to connect to queuing system. Please try again in a few seconds.")
        except ReactomeMQException as e:
            LOGGER.error("Failed to post message to queuing system: " + str(e))
            # update the status
            status = DatasetLoadingStatus(id=loading_id, status="failed", completed=0,
                                    description="Failed to connect to queuing system.")
            storage.set_status(loading_id, encoder.encode(status), data_type="dataset")

            abort(503, "The number of analysis requests is currently too high. Please try again in a few minutes.")
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(503, "Failed to connect to storage system. Please try again in a few minutes.")
    except (socket.timeout, socket.gaierror) as e:
        LOGGER.error("Socket timeout connecting to storage or queuing system: " + str(e))
        abort(503, "Failed to connect to downstream system. Please try again in a few minutes.")
