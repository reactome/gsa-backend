from reactome_analysis_worker.models.gene_set import GeneSet


class GeneSetMapping:
    def __init__(self, gene_set_indices: dict, n_found: dict, n_curated: dict, n_interactors: dict):
        """
        Create a new GeneSetMapping object. All parameters are dicts with the gene set ids as keys.
        :param gene_set_indices: Holds the genes as indexes of the gene in the expression list
        :param n_found: Number of found genes
        :param n_curated: Number of found curated genes
        :param n_interactors: Number of found interactors
        """
        self.gene_set_indices = gene_set_indices
        self.n_found = n_found
        self.n_curated = n_curated
        self.n_interactors = n_interactors

    @staticmethod
    def create_mapping(gene_set: GeneSet, identifiers: list, identifier_mapping: dict):
        """
        Maps the the passed identifiers to the passed gene_set taking the identifier_mapping
        results into consideration.
        :param gene_set: The `GeneSet` to map against
        :param identifiers: The identifiers to map. The resulting mapping will only store the identifiers index.
        :param identifier_mapping: The identifier mapping to use. This must be a dict with the identifier as found
                                   in the `identifiers` list as key and all mappings as values
        :return: A GeneSetMapping object
        """
        # store every identifier's position
        identifier_positions = dict()
        for i in range(len(identifiers)):
            # use 1-based index since R is 1-based
            identifier_positions[identifiers[i]] = i + 1

        # map every found identifier to its position in the identifiers list
        mapped_identifier_positions = dict()
        for org_identifier in identifier_mapping:
            # ignore identifiers that were already removed from the dataset
            if org_identifier not in identifier_positions:
                continue

            position = identifier_positions[org_identifier]
            for mapped_identifier in identifier_mapping[org_identifier]:
                mapped_identifier_positions[mapped_identifier] = position

        # create the mapped gene set object
        gene_set_indices = dict()
        n_found = dict([(gene_set_id, 0) for gene_set_id in gene_set.gene_sets])
        n_curated = dict([(gene_set_id, 0) for gene_set_id in gene_set.gene_sets])
        n_interactors = dict([(gene_set_id, 0) for gene_set_id in gene_set.gene_sets])

        for gene_set_id in gene_set.gene_sets:
            indices = set()

            for gene in gene_set.gene_sets[gene_set_id]:
                if gene in mapped_identifier_positions:
                    gene_index = mapped_identifier_positions[gene]

                    # only count new mappings
                    if gene_index not in indices:
                        n_found[gene_set_id] += 1
                        if gene_set_id in gene_set.interactors and gene in gene_set.interactors[gene_set_id]:
                            n_interactors[gene_set_id] += 1
                        else:
                            n_curated[gene_set_id] += 1

                        indices.add(mapped_identifier_positions[gene])

            # save the updated set
            if len(indices) > 0:
                gene_set_indices[gene_set_id] = indices

        return GeneSetMapping(gene_set_indices, n_found, n_curated, n_interactors)
