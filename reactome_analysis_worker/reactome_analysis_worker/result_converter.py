"""
Converts a REACTOME Analysis System result object into the
result structure expected by the REACTOME Pathway Browser
"""


import copy
import enum
import gzip
import json
import logging
import os
import typing

import urllib3
from reactome_analysis_api.models.analysis_result import AnalysisResult
from reactome_analysis_api.models.analysis_result_reactome_links import AnalysisResultReactomeLinks
from scipy.stats import combine_pvalues, zscore
from statsmodels.stats.multitest import fdrcorrection

from reactome_analysis_worker import util

LOGGER = logging.getLogger(__name__)


class ConversionException(Exception):
    """
    Class to encapsulate exceptions that occurred during the conversion process.
    """
    pass


class ReactomeResultTypes(enum.Enum):
    gsa = "gsa"
    gsa_p_values = "gsa_p"
    gsva = "gsva"


def perform_reactome_gsa(identifiers: typing.Iterable, use_interactors: bool = False, reactome_server: str = "production",
                         include_disease: bool = True) -> dict:
    """
    Use the REACTOME GSA service to perform a complete overrepresentation analysis
    of all identifiers. This result is then used as a "blueprint" to enhance
    it with the results created from this project's analyses.
    :param identifiers: The identifiers to use for the ORA
    :param use_interactors: Indicates whether interactors should be included in the analysis.
    :param reactome_server: The Reactome server to use. Available options are 'production', 'dev', and 'release'
    :param include_disease: Indicates whether to include disease pathways.
    :return: A dict object representing the JSON encoded result
    """
    # use a proxy if set
    proxy = os.getenv("PROXY", None)

    if proxy:
        http = urllib3.ProxyManager(proxy_url=proxy)
        LOGGER.debug("Using proxy " + proxy + " for requests")
    else:
        http = urllib3.PoolManager()

    # set the reactome server url
    reactome_url = util.get_reactome_url(reactome_server)

    body_text = "#Multi-sample_analysis\n" + "\n".join(set(identifiers))

    url = "https://{reactome_url}/AnalysisService/identifiers/projection?interactors={interactors}" \
          "&pageSize=0&page=1&sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1&includeDisease={disease}" \
          .format(reactome_url=reactome_url,  interactors=str(use_interactors).lower(), disease=str(include_disease).lower())

    LOGGER.debug("Performing Reactome GSA analysis with url = " + url)

    request = http.request("POST", url, body=body_text,
                           headers={"content-type": "text/plain"},
                           timeout=30)

    if request.status != 200:
        msg = "Failed get REACTOME Analysis token: {}".format(str(request.status))
        LOGGER.error(msg + " (" + str(request.status) + ")")
        raise ConversionException(msg)

    # get the result object
    reactome_result = json.loads(request.data.decode("utf-8"))

    # extract the token
    try:
        token = reactome_result["summary"]["token"]
    except Exception:
        msg = "REACTOME Analysis token missing in REACTOME response"
        LOGGER.error(msg)
        raise ConversionException(msg)

    # get the actual JSON file
    file_url = "https://{}/AnalysisService/download/{}/result.json".format(reactome_url, token)
    download_request = http.request("GET", file_url)

    if download_request.status != 200:
        msg = "Failed to download REACTOME JSON blueprint"
        LOGGER.error(msg + " (" + str(request.status) + ")")
        raise ConversionException(msg)

    reactome_blueprint = json.loads(download_request.data.decode("utf-8"))

    return reactome_blueprint


def submit_result_to_reactome(result: AnalysisResult, result_type: ReactomeResultTypes, reactome_blueprint: dict,
                              min_p: float = 0.05, reactome_server: str = "production") -> AnalysisResultReactomeLinks:
    """
    Submits the passed AnalysisResult object to Reactome and retrieves the matching token for the stored result.

    :param result: The AnalysisResult to convert and submit.
    :param result_type: The type of result to be converted. Possible values are "gsa", "gsa_p_values" (p-values of the
                        GSA result), and "gsva".
    :param reactome_blueprint: The Reactome result object of the ORA analysis fetched using the `perform_reactome_gsa`
                               function.
    :param min_p: The minimum p-value to consider a pathway as significantly regulated. This is only used in the "gsa"
                  approach.
    :param reactome_server: The Reactome server to use. Available options are 'production', 'dev', and 'release'
    :return: An `AnalysisResultReactomeLinks` object representing the link
    """
    # convert the result
    if result_type == ReactomeResultTypes.gsa:
        converted_result = _convert_gsa_result(result=result, reactome_blueprint=reactome_blueprint, min_p=min_p,
                                               use_p=False)
    elif result_type == ReactomeResultTypes.gsa_p_values:
        converted_result = _convert_gsa_result(result=result, reactome_blueprint=reactome_blueprint, min_p=min_p,
                                               use_p=True)
    elif result_type == ReactomeResultTypes.gsva:
        converted_result = _convert_gsva_result(result=result, reactome_blueprint=reactome_blueprint)
    else:
        raise Exception("Invalid result_type '{}' passed to submit_result_to_reactome. Valid values are gsa, gsa_p, "
                        "or gsva")

    # upload the result to Reactome
    proxy = os.getenv("PROXY", None)

    if proxy:
        http = urllib3.ProxyManager(proxy_url=proxy)
        LOGGER.debug("Using proxy " + proxy + " for requests")
    else:
        http = urllib3.PoolManager()

    # use the file upload end-point as it supports gziped content
    reactome_url = util.get_reactome_url(reactome_server)

    url = "https://{}/AnalysisService/import/form".format(reactome_url)
    # compress the data - this speeds up the upload dramatically
    compressed_data = gzip.compress(json.dumps(converted_result).encode("UTF-8"))

    request = http.request("POST", url,
                           headers={"accepts": "application/json"},
                           fields={"file": ("pathway.json.gz", compressed_data)},
                           timeout=30)

    # check the return code
    if request.status != 200:
        msg = "Failed upload result to REACTOME: " + request.data.decode("UTF-8")
        LOGGER.error(msg + " (" + str(request.status) + ")")
        raise ConversionException(msg)

    # get the token
    result_object = json.loads(request.data.decode('UTF-8'))

    # add the name and description
    if result_type == ReactomeResultTypes.gsa:
        name = "Gene Set Analysis Summary"
        description = "Overview over all submitted datasets showing significantly and non-significantly up- and down- " \
                      "regulated pathways"
    elif result_type == ReactomeResultTypes.gsa_p_values:
        name = "Gene Set Analysis Statistics"
        description = "P-values of all pathways and datasets. This visualisation does not contain information on the " \
                      "direction of the regulated pathway (up or down)."
    elif result_type == ReactomeResultTypes.gsva:
        name = "Gene Set Variation Summary"
        description = "GSVA-derived pathway-level expression for every sample across all datasets. Protein-/gene-level " \
                      "expression is shown as z-score normalised values (across the respective dataset)."
    else:
        name = ""
        description = ""

    return AnalysisResultReactomeLinks(
        url="https://{}/PathwayBrowser/#/DTAB=AN&ANALYSIS={}".format(reactome_url, result_object["token"]),
        token=result_object["token"],
        name=name,
        description=description)


def _get_identifier_changes(identifier_fcs: list, all_identifiers: set, min_p: float = 0.05, return_p: bool = False) -> dict:
    """
    Creates a dict with the identifier as key and the regulation type of the various
    datasets as values. If an identifier is not observed in a dataset, 0 is added. Identifier
    regulation is encoded as follows:

      * 2: significant up-regulated
      * 1: non-significantly up-regulated
      * 0: not found
      * -1: non-significantly down-regulated
      * -2: significantly down-regulated

    :param identifier_fcs: A list of numpy arrays one for each dataset
    :param all_identifiers: A set containing all identifiers
    :param min_p: The minimum p-value to consider a protein / gene significantly regulated.
    :param return_p: If set, the protein FDR values are returned instead of the logFC
    :return: A dict with the identifier as key and all expression values as a list
    """
    identifier_changes = dict([(identifier, list()) for identifier in all_identifiers])

    for identifier_fc in identifier_fcs:
        processed_identifiers = set()

        try:
            identifier_index = identifier_fc.dtype.names.index("Identifier")
            fc_index = identifier_fc.dtype.names.index("logFC")
            p_index = identifier_fc.dtype.names.index("adjPVal")
        except ValueError:
            raise ConversionException("Missing required field in fold-change table")

        for identifier_row in identifier_fc:
            identifier = identifier_row[identifier_index]

            if return_p:
                value = float(identifier_row[p_index])
            else:
                fold_change = float(identifier_row[fc_index])
                is_sig = float(identifier_row[p_index]) <= min_p

                if fold_change > 0 and is_sig:
                    value = 2
                elif fold_change > 0:
                    value = 1
                elif fold_change < 0 and is_sig:
                    value = -2
                elif fold_change < 0:
                    value = -1
                else:
                    value = 0

            identifier_changes[identifier].append(value)
            processed_identifiers.add(identifier)

        # add 0 for missing identifiers
        for identifier in all_identifiers:
            if identifier not in processed_identifiers:
                if return_p:
                    identifier_changes[identifier].append(1)
                else:
                    identifier_changes[identifier].append(0)

    return identifier_changes


def _get_gsva_pathway_expression(result: AnalysisResult) -> dict:
    """
    Extracts the GSVA expression values per pathway and returns the concatenated values
    across all datasets as a dict with the pathway id as key and the expression values in a list as value.
    :param result: An AnalysisResult object
    :return: A dict with the pathway id as key and expression values as value
    """
    # get all observed pathways
    all_pathways = set()
    pathway_fcs = list()

    for dataset in result.results:
        pathway_fc = util.string_to_array(dataset.pathways)
        pathway_fcs.append(pathway_fc)
        pathway_id = pathway_fc.dtype.names.index("Pathway")
        dataset_pathways = [row[pathway_id] for row in pathway_fc]
        all_pathways.update(dataset_pathways)

    # Use the GSVA score as pathway expression value of all experiments
    gsva_expr_per_pathway = dict([(pathway_id, list()) for pathway_id in all_pathways])

    for pathway_fc in pathway_fcs:
        processed_pathways = set()
        # -2 since the pathway id + name column should be ignored
        n_samples = len(pathway_fc.dtype.names) - 2

        for pathway_row in pathway_fc:
            pathway_id = pathway_row[0]
            gsva_expr_per_pathway[pathway_id] += pathway_row.tolist()[2:]
            processed_pathways.add(pathway_id)

        # add missing values
        for pathway_id in all_pathways:
            if pathway_id not in processed_pathways:
                gsva_expr_per_pathway[pathway_id] += [0] * n_samples

    return gsva_expr_per_pathway


def _get_identifier_zscores(result: AnalysisResult) -> dict:
    """
    Extract the expression values for every identifier and returns them as a dict with the
    identifier as a key and the expression values across all datasets as value (list). Expression
    values are z-score normalised.
    :param result: An AnalysisResult object
    :return: A dict with the pathway as key and the expression values as value
    """

    # get the observed proteins
    all_identifiers = set()
    identifier_fcs = list()

    for dataset in result.results:
        if not dataset.fold_changes:
            raise ConversionException("Fold-change data missing in dataset '{}'".format(dataset.name))

        identifier_fc = util.string_to_array(dataset.fold_changes)
        identifier_fcs.append(identifier_fc)

        identifier_index = identifier_fc.dtype.names.index("Identifier")
        dataset_identifiers = [row[identifier_index] for row in identifier_fc]

        all_identifiers.update(dataset_identifiers)

    # Use the z-score across the genes as identifier expression values
    z_scores_per_identifier = dict([(identifier, list()) for identifier in all_identifiers])

    for identifier_fc in identifier_fcs:
        processed_identifiers = set()
        n_samples = len(identifier_fc.dtype.names) - 1

        for identifier_row in identifier_fc:
            identifier = identifier_row[0]
            expression_values = identifier_row.tolist()[1:]
            z_scores = zscore(expression_values)

            z_scores_per_identifier[identifier] += z_scores.tolist()
            processed_identifiers.add(identifier)

        # add the missing values
        for identifier in all_identifiers:
            if identifier not in processed_identifiers:
                z_scores_per_identifier[identifier] += [0] * n_samples

    return z_scores_per_identifier


def _convert_gsva_result(result: AnalysisResult, reactome_blueprint: dict) -> dict:
    """
    Adds the data of the passed AnalysisResult to the passed reactome_blueprint.

    **Note:** Pathway-level p-values are not changed since GSVA approaches do not provide these.

    :param result: The AnalysisResult to draw the data from
    :param reactome_blueprint: The result retrieved from the REACTOME ORA analysis
    :return: The adapted reactome result as a dict.
    """
    reactome_blueprint = copy.deepcopy(reactome_blueprint)

    # get pathway and identifier expression values
    gsva_expr_per_pathway = _get_gsva_pathway_expression(result)
    z_scores_per_identifier = _get_identifier_zscores(result)

    # get the min and max expression values
    all_expr_values = [value for values in gsva_expr_per_pathway.values() for value in values]
    all_expr_values += [value for values in z_scores_per_identifier.values() for value in values]
    min_expr = min(all_expr_values)
    max_expr = max(all_expr_values)

    # set the type
    reactome_blueprint["summary"]["type"] = "GSVA"
    reactome_blueprint["summary"]["sampleName"] = "Multi-sample analysis"

    # create the column names as dataset name + sample name
    column_names = list()

    for dataset in result.results:
        # get the first line of the pathway values (samples should have the same sort-order)
        pathway_header = dataset.pathways[0:dataset.pathways.index("\n")]
        # ignore the first two columns since they contain the pathway id + name
        sample_names = [dataset.name + "-" + field.strip() for field in pathway_header.split("\t")[2:]]
        column_names += sample_names

    # add the dataset names as column names
    reactome_blueprint["expressionSummary"]["columnNames"] = column_names
    reactome_blueprint["expressionSummary"]["min"] = min_expr
    reactome_blueprint["expressionSummary"]["max"] = max_expr

    # populate the pathway data
    for i in range(0, len(reactome_blueprint["pathways"])):
        pathway_id = reactome_blueprint["pathways"][i]["stId"]

        # add the entity-level expression values
        for identifier_level in ("entities", "interactors"):
            for entity_index in range(0, len(reactome_blueprint["pathways"][i]["data"][identifier_level])):
                org_id = reactome_blueprint["pathways"][i]["data"][identifier_level][entity_index]["id"]

                if org_id not in z_scores_per_identifier:
                    raise ConversionException("Missing expression values for " + org_id)

                reactome_blueprint["pathways"][i]["data"][identifier_level][entity_index]["exp"] = \
                    z_scores_per_identifier[org_id]

        # update the pathway-level expression values
        for resource_index in range(0, len(reactome_blueprint["pathways"][i]["data"]["statistics"])):
            if pathway_id not in gsva_expr_per_pathway:
                raise ConversionException("Missing pathway GSVA information for '" + pathway_id + "'")

            reactome_blueprint["pathways"][i]["data"]["statistics"][resource_index]["exp"] = \
                gsva_expr_per_pathway[pathway_id]

    # populate the "not found" data
    for i in range(0, len(reactome_blueprint["notFound"])):
        identifier = reactome_blueprint["notFound"][i]["id"]

        if identifier not in z_scores_per_identifier:
            raise ConversionException("Missing expression data for " + identifier)

        # add the expression data
        reactome_blueprint["notFound"][i]["exp"] = z_scores_per_identifier[identifier]

    return reactome_blueprint


def _convert_gsa_result(result: AnalysisResult, reactome_blueprint: dict, min_p: float = 0.05, use_p: bool = False) -> dict:
    """
    Adds the data of the passed AnalysisResult to the passed reactome_blueprint.
    :param result: The AnalysisResult to draw the data from
    :param reactome_blueprint: The result retrieved from the REACTOME ORA analysis
    :param min_p: The minimum p-value in order to consider a pathway as significantly regulated.
    :param use_p: If set, p-values instead of fold-changes / pathway direction are used as "expression values"
    :return: The adapted reactome result as a dict.
    """
    reactome_blueprint = copy.deepcopy(reactome_blueprint)

    # get all observed pathways
    all_pathways = set()
    pathway_fcs = list()

    for dataset in result.results:
        pathway_fc = util.string_to_array(dataset.pathways)
        pathway_fcs.append(pathway_fc)
        pathway_id = pathway_fc.dtype.names.index("Pathway")
        dataset_pathways = [row[pathway_id] for row in pathway_fc]
        all_pathways.update(dataset_pathways)

    # get the pathway-level changes
    pathway_expr = _get_pathway_changes(pathway_fcs, all_pathways=all_pathways, min_p=min_p, return_p=use_p)
    pathway_p = _get_pathway_p_values(pathway_fcs)

    # get the observed proteins
    all_identifiers = set()
    identifier_fcs = list()

    for dataset in result.results:
        if not dataset.fold_changes:
            raise ConversionException("Fold-change data missing in dataset '{}'".format(dataset.name))

        identifier_fc = util.string_to_array(dataset.fold_changes)
        identifier_fcs.append(identifier_fc)

        identifier_index = identifier_fc.dtype.names.index("Identifier")
        dataset_identifiers = [row[identifier_index] for row in identifier_fc]

        all_identifiers.update(dataset_identifiers)

    # get the gene-/protein-level changes
    identifier_expr = _get_identifier_changes(identifier_fcs, all_identifiers, return_p=use_p)

    # set the type
    reactome_blueprint["summary"]["type"] = "GSA_STATISTICS" if use_p else "GSA_REGULATION"
    reactome_blueprint["summary"]["sampleName"] = "Multi-sample analysis"

    # add the dataset names as column names
    reactome_blueprint["expressionSummary"]["columnNames"] = [dataset.name for dataset in result.results]
    reactome_blueprint["expressionSummary"]["min"] = 0 if use_p else -2
    reactome_blueprint["expressionSummary"]["max"] = 1 if use_p else 2

    # populate the pathway data
    for i in range(0, len(reactome_blueprint["pathways"])):
        pathway_id = reactome_blueprint["pathways"][i]["stId"]

        # add the entity-level expression values
        for identifier_level in ("entities", "interactors"):
            for entity_index in range(0, len(reactome_blueprint["pathways"][i]["data"][identifier_level])):
                org_id = reactome_blueprint["pathways"][i]["data"][identifier_level][entity_index]["id"]

                if org_id not in identifier_expr:
                    raise ConversionException("Missing expression values for " + org_id)

                reactome_blueprint["pathways"][i]["data"][identifier_level][entity_index]["exp"] = identifier_expr[org_id]

        # update the statistics
        for resource_index in range(0, len(reactome_blueprint["pathways"][i]["data"]["statistics"])):
            if pathway_id not in pathway_expr:
                raise ConversionException("Missing pathway regulation information for '" + pathway_id + "'")
            if pathway_id not in pathway_p:
                raise ConversionException("Missing p-value for pathway '" + pathway_id + "'")

            reactome_blueprint["pathways"][i]["data"]["statistics"][resource_index]["entitiesPValue"] = \
                pathway_p[pathway_id]["p"]
            reactome_blueprint["pathways"][i]["data"]["statistics"][resource_index]["entitiesFDR"] = \
                pathway_p[pathway_id]["fdr"]

            reactome_blueprint["pathways"][i]["data"]["statistics"][resource_index]["exp"] = pathway_expr[pathway_id]

    # populate the "not found" data
    for i in range(0, len(reactome_blueprint["notFound"])):
        identifier = reactome_blueprint["notFound"][i]["id"]

        if identifier not in identifier_expr:
            raise ConversionException("Missing expression data for " + identifier)

        # add the expression data
        reactome_blueprint["notFound"][i]["exp"] = identifier_expr[identifier]

    return reactome_blueprint


def _get_pathway_p_values(pathway_fcs: list) -> dict:
    """
    Creates one combine p-value per pathway. If pathways are regulated in
    the same direction, p-values are joined using Fisher's method. If they
    are not, the higher (= worse) p-value is inverted (1 - x), then Fisher's
    method is used.
    :param pathway_fcs: A list containing the pathway tables.
    :return: A dictionary with the pathway ids as key, and a second dict with the "p"-value and "fdr" as values.
    """
    # get all p-values and directions from the datasets
    pathway_ps = dict()

    for pathway_fc in pathway_fcs:
        pathway_index = pathway_fc.dtype.names.index("Pathway")
        direction_index = pathway_fc.dtype.names.index("Direction")
        p_index = pathway_fc.dtype.names.index("PValue")

        for pathway_row in pathway_fc:
            pathway = pathway_row[pathway_index]
            direction = pathway_row[direction_index]
            p = float(pathway_row[p_index])

            if pathway not in pathway_ps:
                pathway_ps[pathway] = list()

            pathway_ps[pathway].append({"direction": direction, "p": p})

    # combine the p-values
    pathway_ids = list()
    combined_p_values = list()

    for pathway in pathway_ps:
        # save the pathway id to align with the combined_p_values array
        pathway_ids.append(pathway)

        # if there is only one entry, use this p-value
        if len(pathway_ps[pathway]) == 1:
            combined_p_values.append(pathway_ps[pathway][0]["p"])
        else:
            # get the lowest p
            min_p = min([entry["p"] for entry in pathway_ps[pathway]])
            min_entry = [entry for entry in pathway_ps[pathway] if entry["p"] == min_p][0]

            # if the direction is different to the lowest, use the inverted p-value
            p_values = list()

            for entry in pathway_ps[pathway]:
                if entry["direction"] == min_entry["direction"]:
                    p_values.append(entry["p"])
                else:
                    p_values.append(1 - entry["p"])

            # combine the p-values using Fisher's method
            combined_p = combine_pvalues(p_values)[1]
            combined_p_values.append(combined_p)

    # correct the p-values
    corrected_p_values = fdrcorrection(combined_p_values)[1]

    # create the result dict
    pathway_p_values = dict()

    for i in range(0, len(pathway_ids)):
        pathway_p_values[pathway_ids[i]] = {"p": combined_p_values[i], "fdr": corrected_p_values[i]}

    return pathway_p_values


def _get_pathway_changes(pathway_fcs: list, all_pathways: set, min_p: float, return_p: bool = False):
    """
    Creates a dict with the pathway id as key and a list as value that
    holds the pathway expression values for every dataset. In this implementation,
    the pathway expression represents the direction of change:

      * 2: significant up-regulated
      * 1: non-significantly up-regulated
      * 0: not found
      * -1: non-significantly down-regulated
      * -2: significantly down-regulated

    :param pathway_fcs: A list of numpy arrays holding the pathway information.
    :param all_pathways: A set containing all pathway ids
    :param min_p: The minimum FDR to consider a pathway significantly regulated
    :param return_p: If set, the pathway's adjusted p-value is returned instead of the direction of change.
                     Missing pathways get a value of "1".
    :return: A dict with the pathway id as key and the expression values per sample as a list
    """
    # initialize the dict to hold the pathway changes
    pathway_changes = dict([(pathway, list()) for pathway in all_pathways])

    for pathway_fc in pathway_fcs:
        # store the column index for the required columns
        try:
            direction_column = pathway_fc.dtype.names.index("Direction")
            fdr_column = pathway_fc.dtype.names.index("FDR")
            pathway_index = pathway_fc.dtype.names.index("Pathway")
        except ValueError:
            raise ConversionException("Missing required pathway column")

        # process every pathway
        processed_pathways = set()

        for pathway_row in pathway_fc:
            pathway_id = pathway_row[pathway_index]
            direction = pathway_row[direction_column]
            fdr = float(pathway_row[fdr_column])

            # make sure no duplicate pathways exist
            if pathway_id in processed_pathways:
                raise ConversionException("Duplicate entries for pathway '" + pathway_id + "'")

            # add the information whether a pathway is
            if return_p:
                value = fdr
            else:
                if direction == "Up" and fdr <= min_p:
                    value = 2
                elif direction == "Up":
                    value = 1
                elif direction == "Down" and fdr <= min_p:
                    value = -2
                elif direction == "Down":
                    value = -1
                else:
                    value = 0

            pathway_changes[pathway_id].append(value)

            # track all processed pathways
            processed_pathways.add(pathway_id)

        # at 0 expression for non-observed pathways
        for missing_pathway in all_pathways:
            if missing_pathway not in processed_pathways:
                if return_p:
                    pathway_changes[missing_pathway].append(1)
                else:
                    pathway_changes[missing_pathway].append(0)

    return pathway_changes
