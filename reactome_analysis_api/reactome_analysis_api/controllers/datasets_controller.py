import logging
import prometheus_client
from flask import abort, Response, current_app
import uuid
import socket
import json

from reactome_analysis_api.encoder import JSONEncoder
from reactome_analysis_api.models.dataset_loading_status import DatasetLoadingStatus  # noqa: E501
from reactome_analysis_api.models.external_data import ExternalData  # noqa: E501
from reactome_analysis_api import util
from reactome_analysis_utils.reactome_mq import ReactomeMQ, ReactomeMQException, DATASET_QUEUE
from reactome_analysis_utils.reactome_storage import ReactomeStorage, ReactomeStorageException
from reactome_analysis_utils.models.dataset_request import DatasetRequest, DatasetRequestParameter
from reactome_analysis_api.models.external_datasource import ExternalDatasource, ExternalDatasourceParameters
from reactome_analysis_api.searcher.public_data_searcher import PublicDatasetSearcher
from reactome_analysis_api.models import data_search_result
from reactome_analysis_api.models import parameter

LOGGER = logging.getLogger(__name__)
DATASET_LOADING_COUNTER = prometheus_client.Counter("reactome_api_loading_datasets",
                                                    "External datasets loaded.", ["resource"])

DATASET_SEARCH_COUNTER = prometheus_client.Counter(
    "reactome_api_dataset_searches", "Number of searches performed.")


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
        ExternalData(id="EXAMPLE_RIBO_SEQ", title="Ribo seq data example", type="riboseq_counts",
                     description="Ribo seq data requires Ribo seq and RNA seq data.",
                     group="RIBO_EXAMPLES")
    ]


def get_data_sources():  # noqa: E501
    """Lists the available external data sources

     # noqa: E501


    :rtype: ExternalDatasource
    """
    return [
        ExternalDatasource(id="example_datasets", name="Example datasets",
                           description="Example datasets to quickly test the application.",
                           url="https://reactome.org/gsa",
                           parameters=[
                               ExternalDatasourceParameters(name="dataset_id", display_name="Dataset Id",
                                                            type="string", description="Identifier of the dataset",
                                                            required=True)
                           ]),
        ExternalDatasource(id="ebi_gxa", name="Expression Atlas",
                           description="EBI's Expression Atlas resource for consistently reprocessed 'omics data.",
                           url="https://www.ebi.ac.uk/gxa/home",
                           parameters=[
                               ExternalDatasourceParameters(name="dataset_id", type="string", display_name="Dataset Id",
                                                            description="Identifier of the dataset", required=True)
                           ]),
        ExternalDatasource(id="ebi_sc_gxa", name="Single Cell Expression Atlas",
                           description="EBI's Single Cell Expression Atlas resource for consistently reprocessed scRNA-seq data.",
                           url="https://www.ebi.ac.uk/gxa/sc/home",
                           parameters=[
                               ExternalDatasourceParameters(name="dataset_id", display_name="Dataset Id",
                                                            type="string", description="Identifier of the dataset",
                                                            required=True),
                               ExternalDatasourceParameters(name="k", type="int", display_name="K",
                                                            description="Parameter k used to create the cell clusters",
                                                            required=True),
                           ]),
        ExternalDatasource(id="grein", name="GREIN Data",
                           description="GREIN is an NCBI project that consistently reprocesses RNA-seq data from GEO.",
                           url="http://www.ilincs.org/apps/grein/?gse=",
                           parameters=[
                               ExternalDatasourceParameters(name="dataset_id", display_name="Dataset Id",
                                                            type="string", description="Identifier of the dataset",
                                                            required=True)]),
        ExternalDatasource(id="geo_microarray", name="GEO query",
                           description="Uses 'GEO query' to load datasets directly from GEO. Primarily supports microarray data.",
                           parameters=[
                               ExternalDatasourceParameters(name="dataset_id", display_name="Dataset Id",
                                                            type="string", description="Identifier of the dataset",
                                                            required=True)])
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

        status = storage.get_status(
            analysis_identifier=loadingId, data_type="dataset")

        if status is None:
            LOGGER.debug(
                "Unknown identifier passed to get_status: " + loadingId)
            abort(404, "Unknown identifier")
        else:
            # return a Response object to prevent connexion from
            # de-serializing the object into a JSON object
            return Response(response=status, status=200, headers={"content-type": "application/json"})
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(
            503, "Failed to connect to storage system. Please try again in a few minutes.")


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
        abort(
            503, "Failed to connect to storage system. Please try again in a few minutes.")


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
        status = DatasetLoadingStatus(
            id=loading_id, status="running", completed=0, description="Queued")
        storage.set_status(loading_id, encoder.encode(
            status), data_type="dataset")

        # convert the parameters
        request_parameters = list()

        for dict_param in parameters:
            request_parameters.append(DatasetRequestParameter(
                name=dict_param["name"], value=dict_param["value"]))

        # create the request
        request = DatasetRequest(
            loading_id=loading_id, resource_id=resourceId, parameters=request_parameters)

        try:
            queue = ReactomeMQ(queue_name=DATASET_QUEUE)
            queue.post_analysis(analysis=request.to_json(),
                                method="DatasetLoading")
            LOGGER.debug("Dataset process " +
                         loading_id + " submitted to queue")
            queue.close()

            DATASET_LOADING_COUNTER.labels(resource=resourceId).inc()

            return loading_id
        except socket.gaierror as e:
            # update the status
            LOGGER.error("Failed to connect to queuing system: " + str(e))
            status = DatasetLoadingStatus(id=loading_id, status="failed", completed=0,
                                          description="Failed to connect to queuing system.")
            storage.set_status(loading_id, encoder.encode(
                status), data_type="dataset")

            abort(
                503, "Failed to connect to queuing system. Please try again in a few seconds.")
        except ReactomeMQException as e:
            LOGGER.error("Failed to post message to queuing system: " + str(e))
            # update the status
            status = DatasetLoadingStatus(id=loading_id, status="failed", completed=0,
                                          description="Failed to connect to queuing system.")
            storage.set_status(loading_id, encoder.encode(
                status), data_type="dataset")

            abort(
                503, "The number of analysis requests is currently too high. Please try again in a few minutes.")
    except ReactomeStorageException as e:
        LOGGER.error("Failed to connect to redis: " + str(e))
        abort(
            503, "Failed to connect to storage system. Please try again in a few minutes.")
    except (socket.timeout, socket.gaierror) as e:
        LOGGER.error(
            "Socket timeout connecting to storage or queuing system: " + str(e))
        abort(503, "Failed to connect to downstream system. Please try again in a few minutes.")

def download_dataset(datasetId, format = None):
    """Download a previously loaded dataset

    :param datasetId: The dataset's id to download
    :type datasetId: string
    :param format: The format in which the data should be presented (xlsx, meta, expr)
    :type format: string
    """
    try:
        # make sure the format is set
        if not format:
            abort(400, "Missing required paramter 'format'")

        storage = ReactomeStorage()

        # check if the dataset exists
        if format == "meta":
            summary = storage.get_request_data_summary(token=datasetId)

            # convert to TSV
            summary = json.loads(summary)

            tsv_string = []

            # first column is the samples
            header_string = "Sample Id\t"
            header_string += "\t".join( [field["name"] for field in summary["sample_metadata"]] )
            tsv_string.append(header_string)

            # add the field data
            for index, sample in enumerate(summary["sample_ids"]):
                sample_string = sample + "\t"

                # add the fields
                sample_string += "\t".join( [str(field["values"][index]) for field in summary["sample_metadata"]] )

                tsv_string.append(sample_string)

            return Response(response="\n".join(tsv_string), status=200, headers={"content-type": "text/plain", 
                                                                      "content-disposition": f"attachment; filename=\"{datasetId}_meta.tsv\""})
        elif format == "expr":
            expression_data = storage.get_request_data(token=datasetId)

            # fix special character encoding issue
            expression_data = expression_data.replace("\\n", "\n")
            expression_data = expression_data.replace("\\t", "\t")

            return Response(response=expression_data, status=200, headers={"content-type": "text/plain", 
                                                                      "content-disposition": f"attachment; filename=\"{datasetId}_expr.tsv\""})
        else:
            abort(404, "Unsupported format passed.")
    except ReactomeStorageException:
        abort(404, "Failed to retrieve dataset. Make sure the dataset was successfully loaded beforehand.")

def get_search_species():  # noqa: E501
    """Returns the available species presented in the available datasets.
     # noqa: E501
    :rtype: list
    """
    try:
        species_list = current_app.public_searcher.get_species()

        return species_list
    except FileNotFoundError as e:
        LOGGER.error(f"Loading species failed.")

        abort(500, "Failed to load available species.")


def search_data(keywords, species=None):  # noqa: E501
    """Key search for public datasets across multiple resources. This function supports all
        resources that are supported by the features to load external datasets.

    :param keywords: The space delimited keywords to search for. Keywords are combined using the logical AND.
    :type keywords: string
    :param species: If set, only samples for this species are being returned.
    :type species: string
    """
    try:
        search_response = current_app.public_searcher.index_search(
            keywords.split(" "), species)
    except Exception as e:
        LOGGER.error(f"Search failed: {keywords}")
        LOGGER.exception(e)
        abort(500, "Internal server error. Search failed.")

    DATASET_SEARCH_COUNTER.inc()

    # convert to a list of search result objects
    search_response_list = list()

    for search_result in search_response:
        # convert the loading parameters
        loading_parameters = list()

        original_parameters = json.loads(search_result["loading_parameters"])
        for param_name in original_parameters.keys():
            loading_parameters.append(parameter.Parameter(name=param_name, value=original_parameters[param_name]))

        # create the search result object
        search_response_result = data_search_result.DataSearchResult(id=search_result["id"],
                                                                     title=search_result["title"],
                                                                     description=search_result["description"],
                                                                     species=search_result["species"],
                                                                     resource_name=search_result["data_source"],
                                                                     resource_loading_id=search_result["resource_id"], 
                                                                     loading_parameters=loading_parameters,
                                                                     web_link=search_result["web_link"])
        search_response_list.append(search_response_result)

    return search_response_list
