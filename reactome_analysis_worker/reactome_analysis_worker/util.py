"""
A collection of basic utility functions
"""

import json
import logging
import os
from io import StringIO

import numpy as np
import urllib3

LOGGER = logging.getLogger(__name__)


class ConversionException(Exception):
    """
    Simple exception class to handle any exceptions
    during the conversion process.
    """
    pass


class MappingException(Exception):
    """
    Simple exception class to handle any exceptions
    during the identifier mapping process.
    """
    pass


def string_to_array(string: str) -> np.ndarray:
    """
    Convert an expression matrix with escaped line-breaks into
    a numpy array
    :param string: The string to convert.
    :param delimiter: The delimiter to use. Tab by default.
    :return: The resulting numpy array.
    """
    # Unescape the tab and new-line delimiters
    formatted_string = string.replace('\\t', '\t').replace('\\n', '\n')

    # If the string starts with a tab, add a "Gene" header
    if formatted_string[0] == "\t":
        formatted_string = "Gene" + formatted_string

    # Convert to a numpy array
    try:
        array = np.genfromtxt(StringIO(formatted_string), names=True, autostrip=True, delimiter="\t", dtype=None,
                              encoding=None)

        return array
    except Exception as e:
        raise ConversionException(e)


def map_identifiers(identifiers: set, return_all: bool = True) -> dict:
    """
    Map the passed identifiers using REACTOME's mapping service but returns the result as a dict
    with the original identifier as a key and the mappings as values.
    :param identifiers: Identifiers to map.
    :param return_all: Indicates whether all mappings or just the first one should be used.
    :return: A dict with the original identifier as key and the mapping result as values of a list
    """
    proxy = os.getenv("PROXY", None)

    if proxy:
        http = urllib3.ProxyManager(proxy_url=proxy)
        LOGGER.debug("Using proxy " + proxy + " for requests")
    else:
        http = urllib3.PoolManager()

    LOGGER.debug("Mapping identifiers using REACTOME...")

    url = "https://dev.reactome.org/AnalysisService/mapping/projection/?interactors=true"
    request = http.request("POST", url, body="\n".join(set(identifiers)), headers={"content-type": "text/plain"},
                           timeout=5)

    if request.status != 200:
        msg = "Failed to retrieve mappings: Invalid identifiers submitted".format(str(request.status))
        LOGGER.error(msg + " (" + str(request.status) + ")")
        raise MappingException(msg)

    # get the mapping data
    reactome_mappings = json.loads(request.data.decode("utf-8"))

    final_mappings = dict()

    for mapping in reactome_mappings:
        # ignore any unmapped identifiers
        if len(mapping["mapsTo"]) < 1:
            continue

        if not return_all:
            final_mappings[mapping["identifier"]] = [mapping["mapsTo"][0]["identifier"]]
        else:
            final_mappings[mapping["identifier"]] = list()
            for mapped_identifier in mapping["mapsTo"]:
                final_mappings[mapping["identifier"]].append(mapped_identifier["identifier"])

    return final_mappings
