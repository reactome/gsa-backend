import os
import pickle
from copy import deepcopy


class GeneSet:
    """
    Represents a basic gene set from a source.

    :ivar gene_sets: A dict holding all gene sets with the id as key and the members as values of a list
    :ivar gene_set_names: A dict holding additional descriptions / names per gene set
    :ivar gene_set_size: Size of each gene set as a dict with the gene set id as key and the size as value
    :ivar n_curated: Number of curated genes in the gene set as a dict with the gene set id as key and the number
                   of curated genes as value
    :ivar n_interactors: Number of interactors in the gene set as a dict with the gene set id as key and the number
                       of interactors as value
    :ivar interactors: A dict with the pathway id as key and the interactors added to the pathway as values
    """
    def __init__(self, gene_sets: dict, gene_set_names: dict = None):
        """
        Initializes a GeneSet object with the specified gene_sets.
        :param gene_sets: A dict with the pathway name / identifier as key and the members as values of a list.
        :param gene_set_names: A dict with the pathway identifier as key and the name / description as value.
        """
        self.gene_sets = gene_sets

        self.gene_set_names = gene_set_names if gene_set_names else dict()

        # initialize the lengths
        self.gene_set_size = dict([(gene_set_id, len(self.gene_sets[gene_set_id]))
                                   for gene_set_id in self.gene_sets])

        # number of curated genes - all
        self.n_curated = deepcopy(self.gene_set_size)

        # number of interactors - initialized with 0
        self.n_interactors = dict([(gene_set_id, 0) for gene_set_id in self.gene_sets])

        self.interactors = dict()

    @staticmethod
    def create_from_reactome_mapping(mappings: str, species: str = "Homo sapiens"):
        """
        Creates a GeneSet object based on a REACTOME mapping specification.
        :param mappings: A string containing the REACTOME mappings
        :param species: The species to fetch the pathway mappings for.
        :return: The generated GeneSet object
        :raise SyntaxError: In case the mapping specification is invalid.
        """
        pathways = dict()
        pathway_names = dict()

        mapping_lines = mappings.split("\n")

        for mapping in mapping_lines:
            mapping = mapping.strip()

            # ignore empty lines
            if len(mapping) == 0:
                continue

            fields = mapping.split("\t")

            if len(fields) < 6:
                raise SyntaxError("Invalid mapping specification passed. Must contain at least 6 fields.")

            molecule_id = fields[0]
            pathway_id = fields[1]
            pathway_name = fields[3]
            pathway_species = fields[5]

            if pathway_species != species:
                continue

            if pathway_id not in pathways:
                pathways[pathway_id] = set()
                pathway_names[pathway_id] = pathway_name

            pathways[pathway_id].add(molecule_id)

        # create the GeneSet object
        return GeneSet(gene_sets=pathways, gene_set_names=pathway_names)

    def add_interactors(self, interactors: dict):
        """
        Add the specified interactors to the already present gene sets.
        :param interactors: A dict in the form of key interacts with all members of the subsequent list.
        """
        for pathway_id in self.gene_sets.keys():
            pathway_genes = self.gene_sets[pathway_id]
            pathway_interactors = set()

            # get all interactors for that pathway
            for gene in pathway_genes:
                if gene in interactors:
                    pathway_interactors.update(interactors[gene])

            # save the added interactors per pathway
            self.interactors[pathway_id] = set([interactor for interactor in pathway_interactors
                                                if interactor not in pathway_genes])

            self.gene_sets[pathway_id].update(pathway_interactors)

        # update the gene set sizes
        self.gene_set_size = dict([(gene_set_id, len(self.gene_sets[gene_set_id]))
                                   for gene_set_id in self.gene_sets])

        # get the number of interactors by subtracting from the curated size
        self.n_interactors = dict([(gene_set_id, self.gene_set_size[gene_set_id] - self.n_curated[gene_set_id])
                                   for gene_set_id in self.gene_sets])

    def copy(self):
        """
        Returns a copy of this object
        :return: A GeneSet object
        """
        copy = GeneSet(dict())

        copy.gene_sets = deepcopy(self.gene_sets)
        copy.gene_set_size = deepcopy(self.gene_set_size)
        copy.interactors = deepcopy(self.interactors)
        copy.n_curated = deepcopy(self.n_curated)
        copy.n_interactors = deepcopy(self.n_interactors)

        return copy

    def save(self, filename: str):
        """
        Saves the gene set to a file
        :param filename: The target filename. If it exists, an exception will be thrown
        :raise FileExistsError: In case the target filename exists
        """
        if os.path.isfile(filename):
            raise FileExistsError("Target file exists")

        with open(filename, "wb") as binary_writer:
            pickle.dump([self.gene_sets, self.gene_set_size, self.gene_set_names,
                         self.interactors, self.n_curated, self.n_interactors],
                        binary_writer)

    @staticmethod
    def create_from_file(filename: str):
        """
        Creates a GeneSet object from a file that was created using the `save` function.
        :param filename: Path to the filename
        :return: The GeneSet object
        :raise FileNotFoundError
        """
        if not os.path.isfile(filename):
            raise FileNotFoundError("Failed to find GeneSet file " + filename)

        with open(filename, "rb") as binary_reader:
            (gene_sets, gene_set_size, gene_set_names, interactors,  n_curated, n_interactors) = \
                pickle.load(binary_reader)

            gene_set = GeneSet(dict())
            gene_set.gene_sets = gene_sets
            gene_set.gene_set_size = gene_set_size
            gene_set.gene_set_names = gene_set_names
            gene_set.interactors = interactors
            gene_set.n_curated = n_curated
            gene_set.n_interactors = n_interactors

            return gene_set
