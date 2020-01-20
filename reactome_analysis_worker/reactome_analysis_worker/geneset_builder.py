"""
A collection of scripts to create the gene sets used for the analysis based on the
available data from REACTOME.

If this file is called directly, it retrieves REACTOME data for the currently supported
species and stores it as binary files. The scripts evaluates the following environmental variables:
  * `REACTOME_STORAGE_PATH`: Path where the files should be stored (default is the temporary directory).
  * `REACTOME_SOURCE`: Location from where the files should be retrieved (default REACTOME website).
  * `REACTOME_INTERACTOR_FILE`: Location of the file to load the interactors from. By default, the
                                packaged file is being used.
  * `REACTOME_DISEASE_FILE`: Location or URL of the file to load the disease pathways from.

Files to download:

  * Pathway mappings:
    * https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
    * https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt
  * Reaction mappings:
    * https://reactome.org/download/current/UniProt2ReactomeReactions.txt
    * https://reactome.org/download/current/Ensembl2ReactomeReactions.txt
  * Disease pathway ids:
    * https://dev.reactome.org/download/current/HumanDiseasePathways.txt
  * Reaction to pathways:
    * From REACTOME content service - for the specific species
"""

import codecs
import logging
import os
import sys
from tempfile import gettempdir

import urllib3

from reactome_analysis_worker.models.gene_set import GeneSet

LOGGER = logging.getLogger(__name__)


def fetch_reactome_geneset(source: str, species: str) -> GeneSet:
    """
    Creates two GeneSet objects from the REACTOME sources, one for proteomics data (UniProt mappings)
    and one for genomics data (ENSEMBL mappings)
    :param source: May either contain a URL to fetch the files from or a path to load
                   them from directly.
    :param species: The species to filter the pathways for.
    :return: The GeneSet object
    :raise FileNotFoundError: Raised if a directory was specified as source and the pathway file could not be found.
    """
    pathway_files = ("UniProt2Reactome_All_Levels.txt", "Ensembl2Reactome_All_Levels.txt")
    reactome_mappings = ""

    # get the reactome mappings
    for pathway_file in pathway_files:
        if source[0:4] == "http":
            # make sure there is a trailing "/"
            if source[-1] != "/":
                source += "/"

            url = source + pathway_file

            LOGGER.debug("Loading REACTOME mapping from '{}'".format(url))

            http = urllib3.PoolManager()
            request = http.request("GET", url, retries=False)

            LOGGER.debug("Response code = {}".format(request.status))

            if request.status != 200:
                msg = "Failed to retrieve pathway file '{}': Code {}".format(url, request.status)
                LOGGER.error(msg)
                raise FileNotFoundError(msg)

            reactome_mappings += request.data.decode('utf-8')
        else:
            filename = os.path.join(source, pathway_file)

            if not os.path.isfile(filename):
                msg = "Failed to find pathway file '{}'".format(filename)
                LOGGER.error(msg)
                raise FileNotFoundError(msg)

            LOGGER.debug("Loading REACTOME mapping from file '{}'".format(filename))

            with codecs.open(filename=filename, mode="r", encoding='UTF-8') as reader:
                reactome_mappings += reader.read()

    # create the GeneSet object
    try:
        gene_set = GeneSet.create_from_reactome_mapping(mappings=reactome_mappings, species=species)
    except SyntaxError as e:
        msg = str(e)
        LOGGER.error(msg)
        raise e

    return gene_set


def load_reactome_interactors(filename: str) -> dict:
    """
    Loads the interaction data from the REACTOME interaction file
    and returns them as a dict with the key as A interacting with all
    members of the subsequent set. Interactions are filtered based on a
    minimum score of 0.45 (used in Reactome).
    :param filename: File to load the interactors from. By default, the packaged file will be used.
    :return: Dict with molecule id as key and its interactors as value = set
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError("Failed to find interactor file {}".format(filename))

    interactions = dict()

    with open(filename, "r") as reader:
        header = reader.readline()
        header_fields = header.split("\t")

        if len(header_fields) < 2 or header_fields[0] != "A" or header_fields[1] != "B" or header_fields[2] != "Score":
            raise SyntaxError("Invalid interaction file passed. Column 1 must be 'A', column 2 'B', "
                              "and column 3 'Score'")

        for line in reader:
            line = line.strip()

            # ignore empty lines
            if len(line) == 0:
                continue

            fields = line.split("\t")

            score = float(fields[2])

            # ignore every interaction below a score of 0.45
            if score < 0.45:
                continue

            if fields[0] not in interactions:
                interactions[fields[0]] = set()

            interactions[fields[0]].add(fields[1])

    return interactions


def load_reactome_disease_pathways(filename: str) -> list:
    """
    Loads the disease pathways from the specified file.
    :param filename: The path to the disease specification file. May be a URL or local path
    :returns: A list holding the disease pathway ids
    """
    # check whether the file should be downloaded
    if filename[0:4] == "http":
        disease_request = urllib3.PoolManager().request("GET", filename, retries=False)

        if disease_request.status != 200:
            LOGGER.error("Failed to download disease pathway file {} ({})".format(filename, str(disease_request.status)))
            raise FileNotFoundError("Failed to download disease pathway file.")

        disease_string = disease_request.data.decode("utf-8")
    else:
        if not os.path.isfile(filename):
            msg = "Failed to find disease pathway file '{}'".format(filename)
            LOGGER.error(msg)
            raise FileNotFoundError(msg)

        LOGGER.debug("Loading REACTOME diseases from file '{}'".format(filename))

        with codecs.open(filename=filename, mode="r", encoding='UTF-8') as reader:
            disease_string = reader.read()
    
    # process the file line by line - the first column is the pathway id
    pathway_ids = list()

    for line in disease_string.split("\n"):
        fields = line.split("\t")
        pathway_id = fields[0].strip()

        if len(pathway_id) > 0:
            pathway_ids.append(pathway_id)

    return pathway_ids


def generate_pathway_filename(resource: str = "reactome", species: str = "Homo sapiens",
                              contains_interactors: bool = False, contains_disease: bool = True) -> str:
    """
    Generates the filename for the binary file storing the pathway data. The location (directory) for the files
    is determined through the `REACTOME_STORAGE_PATH` environmental variable. By default, the system's temporary
    directory is being used.
    :param resource: Name of the pathway resource
    :param species: Name of the species
    :param contains_interactors: Boolean indicating whether the pathways contain interactors
    :param contains_disease: Boolean indicating whether the pathways contain disease pathways as well.
    :return: A string representing the abosolute path of the file
    """
    storage_location = os.getenv("REACTOME_STORAGE_PATH", gettempdir())

    if contains_disease:
        disease_string = ""
    else:
        disease_string = "_no-dis"

    if contains_interactors:
        interactor_string = "_interactors"
    else:
        interactor_string = ""

    filename = "{resource}_{species}{disease}{interactors}.pkl".format(
        resource=resource.lower().replace(" ", "_"),
        species=species.lower().replace(" ", "_"),
        disease=disease_string,
        interactors = interactor_string
    )

    return os.path.join(os.path.abspath(storage_location), filename)


def main():
    logging.basicConfig(level=logging.DEBUG)

    pathway_source = os.getenv("REACTOME_SOURCE", "https://reactome.org/download/current")
    interactor_file = os.getenv("REACTOME_INTERACTOR_FILE", None)
    disease_pathway_file = os.getenv("REACTOME_DISEASE_FILE", "https://dev.reactome.org/download/current/HumanDiseasePathways.txt")

    if not interactor_file:
        LOGGER.error("Missing required environmental variable REACTOME_INTERACTOR_FILE")
        sys.exit(1)

    # load the interactors
    LOGGER.info("Loading interactors from " + interactor_file)
    interactors = load_reactome_interactors(interactor_file)

    # load the disease pathway ids
    LOGGER.info("Loading disease pathways from " + disease_pathway_file)
    disease_pathways = load_reactome_disease_pathways(disease_pathway_file)

    # get the gene sets
    for species in ["Homo sapiens"]:
        LOGGER.info("Loading gene set from " + pathway_source)
        gene_set = fetch_reactome_geneset(source=pathway_source, species=species)

        # save the file
        target_filename = generate_pathway_filename("reactome", species, contains_interactors=False, contains_disease=True)
        if os.path.exists(target_filename):
            os.unlink(target_filename)

        LOGGER.info("Saving gene set to " + target_filename)
        gene_set.save(target_filename)

        # save the file without disease pathways
        gene_set_no_disease = gene_set.copy()
        gene_set_no_disease.remove_pathways(disease_pathways)

        target_filename = generate_pathway_filename("reactome", species, contains_interactors=False, contains_disease=False)
        if os.path.exists(target_filename):
            os.unlink(target_filename)

        LOGGER.info("Saving gene set to " + target_filename)
        gene_set_no_disease.save(target_filename)

        # add the interactors
        gene_set.add_interactors(interactors)

        target_filename = generate_pathway_filename("reactome", species, contains_interactors=True, contains_disease=True)

        if os.path.isfile(target_filename):
            os.unlink(target_filename)

        LOGGER.info("Saving gene set with interactors to " + target_filename)
        gene_set.save(target_filename)

        # save the file without disease pathways
        gene_set_no_disease = gene_set.copy()
        gene_set_no_disease.remove_pathways(disease_pathways)

        target_filename = generate_pathway_filename("reactome", species, contains_interactors=True, contains_disease=False)
        if os.path.exists(target_filename):
            os.unlink(target_filename)

        LOGGER.info("Saving gene set to " + target_filename)
        gene_set_no_disease.save(target_filename)


if __name__ == "__main__":
    main()
