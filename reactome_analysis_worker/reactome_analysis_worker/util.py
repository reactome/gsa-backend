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


def string_to_array(string: str, first_column_str: bool = True) -> np.ndarray:
    """
    Convert an expression matrix with escaped line-breaks into
    a numpy array
    :param string: The string to convert.
    :param delimiter: The delimiter to use. Tab by default.
    :param first_column_str: If set, the function ensures that the first column
                             is interpreted as a string.
    :return: The resulting numpy array.
    """
    # Unescape the tab and new-line delimiters
    formatted_string = string.replace('\\t', '\t').replace('\\n', '\n')

    # If the string starts with a tab, add a "Gene" header
    if formatted_string[0] == "\t":
        formatted_string = "Gene" + formatted_string

    # Convert to a numpy array
    
    # Note: Due to some bug, the automatic type detection in numpy no longer works
    #       properly in numpy, therefore, types are detected "manually"

    # get the second line
    lines = formatted_string.splitlines()
    if len(lines) < 2:
        raise ConversionException("Insuffcient number of lines found")
    
    fields = lines[1].split("\t")
    field_types = list()

    for field in fields:
        # test for an integer first
        try:
            value = int(field)
            field_types.append("i8")
            continue
        except Exception:
            pass

        # try a float as second option
        try:
            value = float(field)
            field_types.append("f8")
            continue
        except Exception:
            pass

        # use string as last type
        field_types.append("U30")
    try:
        array = np.genfromtxt(StringIO(formatted_string), names=True, autostrip=True, delimiter="\t", dtype=field_types,
                                encoding=None, missing_values="NA", filling_values=0)

        # ensure that the first (identifier) column is a string
        if first_column_str and not str(array.dtype[0]).startswith("<U"):
            dt = array.dtype.descr
            dt[0] = (dt[0][0], '<U15')
            array = array.astype(dt)

        return array
    except Exception as e:
        raise ConversionException(e)


def map_identifiers(identifiers: set, return_all: bool = True, reactome_server: str = "production") -> dict:
    """
    Map the passed identifiers using REACTOME's mapping service but returns the result as a dict
    with the original identifier as a key and the mappings as values.
    :param identifiers: Identifiers to map.
    :param return_all: Indicates whether all mappings or just the first one should be used.
    :param reactome_server: The Reactome server to use. Available options are 'production', 'dev', and 'release'
    :return: A dict with the original identifier as key and the mapping result as values of a list
    """
    # convert all identifiers to string
    identifiers = [str(value) for value in identifiers]

    # make sure that the identifiers are valid
    _check_valid_identifiers(identifiers)

    proxy = os.getenv("PROXY", None)

    if proxy:
        http = urllib3.ProxyManager(proxy_url=proxy)
        LOGGER.debug("Using proxy " + proxy + " for requests")
    else:
        http = urllib3.PoolManager()

    LOGGER.debug("Mapping identifiers using REACTOME...")

    reactome_url = get_reactome_url(reactome_server)
    url = "https://{}/AnalysisService/mapping/projection/?importableOnly=true&interactors=true".format(reactome_url)
    request = http.request("POST", url, body="\n".join(set(identifiers)), headers={"content-type": "text/plain"},
                           timeout=30)

    if request.status != 200:
        msg = "Failed to retrieve mappings: Invalid identifiers submitted"
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


def _check_valid_identifiers(identifiers: list) -> None:
    """
    Checks whether all passed identifiers are valid - or rather
    tries to detect common problems. If a problem is found, an
    exception is raised.
    :param identifiers: The list of identifiers to check
    """
    # check for protein groups
    for identifier in identifiers:
        if ";" in identifier:
            raise MappingException("Invalid gene/protein identifiers submitted. Identifiers contain ';'. Did you submit protein groups? Only single identifiers are supported.")
        if " " in identifier:
            raise MappingException("Invalid gene/protein identifiers submitted. Identifier '{}' contains a ' '.".format(identifier))
        if "," in identifier:
            raise MappingException("Invalid gene/protein identifiers submitted. Identifier '{}' contains a ','.".format(identifier))



def get_reactome_url(reactome_server: str) -> str:
    """
    Returns the URL for the specified Reactome server.
    :param reactome_server: The Reactome server to use. Available options are 'production', 'dev', and 'release'
    :return str The URL for the server.
    """
    reactome_url = "www.reactome.org"

    if reactome_server == "dev":
        reactome_url = "dev.reactome.org"
    if reactome_server == "release":
        reactome_url = "release.reactome.org"

    return reactome_url