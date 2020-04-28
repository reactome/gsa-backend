import unittest
import subprocess
import logging
import requests
import json
import time
import numpy
import click
import os
import sys


"""
Test Specification:
   'tests': [
       {
           'name': Test name,
           'dataset': str,
           'type': ["pathways", "fold_changes"]
           'value': numeric
       }
   ]

"""


logger = logging.getLogger(__name__)


@click.command()
@click.option("-s", "--server", default=None)
@click.argument("filenames", nargs=-1, type=click.Path())
def submit_request(server, filenames):
    logging.basicConfig(level=logging.DEBUG)
    urllib_logger = logging.getLogger("urllib3")
    urllib_logger.setLevel(logging.ERROR)

    for filename in filenames:
        process_file(server, filename)


def process_file(server: str, filename: str) -> None:
    """
    Process a request file
    :param server: Name of the remote server to use.
    :param filename: Path to the file
    """
    # make sure the file exists
    if not os.path.isfile(filename):
        print("Error: {name} does not exist".format(filename))
        sys.exit(1)

    # get the service URL
    service_url = get_url_for_system(server)

    # load the request
    with open(filename, "r") as reader:
        request_data = reader.read()

    request_object = json.loads(request_data)

    # submit the request
    print("Submitting request to ReactomeGSA...")

    analysis_id = run_analysis(request_object, service_url)

    if not analysis_id:
        print("Analysis failed. Exiting...")
        sys.exit(1)

    # get the result object
    result = retrieve_reactome_result(analysis_id, service_url)

    print("Analysis completed. Final id = " + analysis_id)

    print_result_statistics(result)

    if "tests" in request_object:
        run_tests(request_object["tests"], result)


def run_tests(tests: list, result: dict) -> None:
    """
    Run the tests defined in the request object's test structure.
    :param tests: The tests to run in the above defined structure
    :param result: The result object as a dict
    """

    print("\nRunning tests...")

    for the_test in tests:
        if not "type" in the_test or not "name" in the_test:
            logger.error("Invalid test specification found")
            continue

        print("* " + the_test["name"] + "...", end = "")

        if the_test["type"] == "fold_changes":
            # get the result
            dataset = get_dataset(result, the_test["dataset"])
            if not dataset:
                print("Failed.\n  Dataset '{name}' does not exist".format(name=the_test["dataset"]))
                continue
            
            n_fold_changes = len(dataset["fold_changes"].split("\n")) - 1
            if n_fold_changes != int(the_test["value"]):
                print("Failed.\n  Expected {exp} but got {num}"
                .format(exp=str(the_test["value"]), num=str(n_fold_changes)))
                continue
            else:
                print("OK.")

        if the_test["type"] == "pathways":
            # get the result
            dataset = get_dataset(result, the_test["dataset"])
            if not dataset:
                print("Failed.\n  Dataset '{name}' does not exist".format(name=the_test["dataset"]))
                continue
            
            n_pathways = len(dataset["pathways"].split("\n")) - 1
            if n_pathways != int(the_test["value"]):
                print("Failed.\n  Expected {exp} but got {num}"
                .format(exp=str(the_test["value"]), num=str(n_pathways)))
                continue
            else:
                print("OK.")


def get_dataset(result: dict, name: str) -> dict:
    """
    Returns the result for the specified dataset.
    :param result: The analysis result as a dict
    :param name: The dataset's name to extract
    :returns: The extracted datataset
    """
    for dataset in result["results"]:
        if dataset["name"] == name:
            return dataset

    return None


def print_result_statistics(result_obj: dict) -> None:
    """
    Print basic statistics about the retrieved result.
    :param result_obj: The analysis result as a dict object
    """
    print("------------- Analysis Result ------------------")
    print("Reactome Version: " + result_obj["release"])
    print("# of datasets: " + str(len(result_obj["results"])))

    # process all datasets
    for dataset in result_obj["results"]:
        print("\nDataset " + dataset["name"])
        print("  Pathways: " + str(len(dataset["pathways"].split("\n")) - 1))
        print("  Fold-changes: " + str(len(dataset["fold_changes"].split("\n")) - 1))


def retrieve_reactome_result(analysis_id: str, reactome_url: str) -> dict:
    """
    Retrieves the result for the specified analysis from ReactomeGSA.
    :param analysis_id: The analysis id
    :param reactome_url: The URL of the ReactomeGSA system to contact.
    :returns: The result object as a dict
    """
    result_obj = send_request(
        "{base_url}result/{analysis_id}".format(base_url=reactome_url, analysis_id=analysis_id), 
        "Failed to retrieve result.")

    return result_obj


def get_url_for_system(service_name) -> str:
    """
    Loads the services from the sources.json and returns the matching
    URL for the specified service. If not service_name is set, the
    user is prompted to enter one.

    :param service_name: The service's name
    :returns: The url
    """
    # load the available services from the config
    service_file = os.path.join(os.path.dirname(__file__), "sources.json")

    if not os.path.isfile(service_file):
        logger.warn("Failed to find sources.json. Using default list.")

        if service_name == "reactome-prod":
            return "https://gsa.reactome.org/0.1/"
        else:
            print("Error: Failed to find sources.json. Service {name} is not known"
            .format(name=service_name))
            sys.exit(1)

    # load the service list from file
    services = dict()

    with open(service_file, "r") as reader:
        service_list = json.loads(reader.read())

    for service in service_list:
        services[service["name"]] = service["url"]
        if not services[service["name"]].endswith("/"):
            services[service["name"]] += "/"

    # ask the user to choose a service name
    if not service_name:
        service_name = input("Server (" + ", ".join(services.keys()) + "): ")

    # make sure the user supplied service is defined
    if service_name not in services:
        print("Error: Unknown service '{name}'".format(name=service_name))
        sys.exit(1)

    return services[service_name] + "0.1/"


def send_request(url, error_msg, decode_json=True):
    """
    Sends a request to the API. If the request does not respond with a
    200 fail with the respective message. Otherwise, the decoded
    JSON result object is returned.

    :param url: The url to fetch
    :param error_msg: The error message to show if the request fails
    :param decode_json: Decode the result as JSON object

    :returns: The decoded JSON result object
    """
    request_obj = requests.get(url)

    if request_obj.status_code != 200:
        raise Exception(error_msg)

    # decode the JSON data
    if decode_json:
        data = json.loads(request_obj.content)
    else:
        data = request_obj.content.decode()

    return data


def load_expression_atlas(dataset_id: str, service_url: str) -> str:
    """
    Load the specified ExpressionAtlas dataset.
    :param dataset_id: The ExpressionAtlas dataset id
    :param service_url: The ReactomeGSA url
    :returns: The final dataset identifier
    """
    loading_request = requests.post(service_url + "data/load/ebi_gxa", json=[{"name": "dataset_id", "value": dataset_id}])

    if loading_request.status_code != 200:
        raise Exception("Failed to send loading request: {}".format(loading_request.content.decode()))

    loading_token = loading_request.content.decode()

    # get the status
    loading_status = {"status": "running"}

    logger.debug("Waiting for status to change from 'running'...")

    while loading_status["status"] == "running":
        time.sleep(1)
        loading_status = send_request(service_url + "data/status/" + loading_token, "Failed to send status request")

    if loading_status["status"] != "complete":
        raise Exception("Failed to load dataset " + dataset_id)

    return loading_status["dataset_id"]


def run_analysis(request: dict, service_url: str) -> str:
    """
    Run the specified analysis.
    :param request: The request object as a dict
    :param service_url: Reactome's URL
    :returns: The analysis identifier
    """
    # start the analysis
    analysis_request = requests.post(service_url + "analysis", json=request)

    if analysis_request.status_code != 200:
        raise Exception("Failed to submit analysis: " + analysis_request.content.decode())

    # get the status
    analysis_id = analysis_request.content.decode()

    status = send_request(service_url + "status/" + analysis_id, "Failed to retrieve status")

    previous_status = status["description"]

    while status["status"] == "running":
        if status["description"] != previous_status:
            logger.debug(status["description"])
            previous_status = status["description"]

        # wait
        time.sleep(1)

        # update the status
        status = send_request(service_url + "status/" + analysis_id, "Failed to retrieve status")

    if status["status"] != "complete":
        logger.error("Analysis failed: " + status["description"])
        return None

    return analysis_id


if __name__ == "__main__":
    submit_request(None, None)
