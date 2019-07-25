#' Removes unmapped genes
#'
#' Removes all genes/proteins from the expression matrix that do not appear in the mapping result ie. that
#' do not correspond to names in the \code{identifier_mapping} list.
#'
#' @param expression.data data.frame containing genes in rows and samples in columns
#' @param identifier_mapping names list with the genes as names and their mappings as values
remove_unmapped <- function(expression.data, identifier_mapping) {
    expression.data <- expression.data[rownames(expression.data) %in% names(identifier_mapping), ]

    return(expression.data)
}

#' Remove missing / infrequently expressed genes
#'
#' Tests in how many samples of the main comparison group a gene / proteins
#' is found. Removes all entities that are found in less than \code{max.missing.values}
#' of the most frequent group.
#'
#' @param expression.data data.frame containing genes in rows and samples in columns
#' @param sample.data data.frame with sample metadata. Each sample is represented as one row.
#'        The data in the \code{analysis_group} column will be evaluated to estimate the frequency
#'        of a gene.
#' @param max.missing.values The maximum allowed relative proprtion of missing values (0 - 1).
remove_missing <- function(expression.data, sample.data, max.missing.values) {
    rel.missing.values <- apply(expression.data, 1, function(gene.expr) {
        missing.values <- by(gene.expr, sample.data$analysis_group, function(x) sum(x == 0))
        # get the relative numbers per group
        missing.values <- as.numeric(missing.values) / as.numeric(table(sample.data$analysis_group))
        return(min(missing.values))
    })

    # every gene must be observed in at least max.missing.values of the samples of one group
    filtered.expressions <- expression.data[rel.missing.values < max.missing.values,]

    return(filtered.expressions)
}

#' Replicates \code{\link[limma]{ids2indices}} while supporting multiple mappings
#'
#' When mapping identifiers, some identifiers often map to multiple entries in the target namespace.
#' This function replicates \code{\link[limma]{ids2indices}} but adds support for these multiple mappings.
#'
#' If a protein / gene maps to multiple entries in the target namespace, the gene / protein will be referenced
#' once in every pathway where any of these mapped entries occurs.
#'
#' @param identifier_mapping A named list with the original gene/protein identifier as key and all mappings as
#'        a vector
#' @param gene_set A named list with pathway ids as names and their member genes / proteins as vactors
#' @param expression.data \code{data.frame} containing the expression data for genes / proteins (rows) per sample
#'        (column)
ids2indices_multiple <- function(identifier_mapping, gene_set, expression.data) {
    # make sure the identifier_mapping has the same sort order as the rows in the expression.data
    identifier_mapping <- identifier_mapping[rownames(expression.data)]

    # map every identifier to a position
    identifier_position = c()
    for (i in 1:length(identifier_mapping)) {
        identifier_position[identifier_mapping[[i]]] = i
    }

    # get the gene position per gene set
    gene.indices <- lapply(gene_set, function(x) {
        existing.genes <- x[x %in% names(identifier_position)]
        positions <- as.numeric(identifier_position[existing.genes])
        positions <- unique(positions)
        return(positions)
    })

    # only keep gene sets with mappings
    gene.indices <- gene.indices[lengths(gene.indices) != 0L]

    return(gene.indices)
}

create_design <- function(sample.data, group_1) {
    # re-level the comparison group to use group_1 as the reference
    if (!"analysis_group" %in% colnames(sample.data)) {
        stop("Error: analysis_group not found in the sample.data")
    }

    all_levels <- unique(as.character(sample.data$analysis_group))
    all_levels <- all_levels[all_levels != group_1]
    sample.data$analysis_group <- factor(sample.data$analysis_group, levels = c(group_1, all_levels))

    # create the design taking all properties into consideration
    design <- model.matrix(~ 0 + ., data=sample.data)

    return(design)
}

#' Converts a data.frame to a string representation
#'
#' A data.frame is converted into a single string using
#' `\t` (the characters, not tab) as field delimiter and
#' `\n` (the characters, not newline) as line delimiter
#'
#' @param data The data.frame to convert
#' @param keep_rownames If set, the rownames will be stored as this column name
#' @return A string representing the passed data.frame
data_frame_as_string <- function(data, rowname_column = NULL) {
    text.connection <- textConnection(NULL, open="w")
    if (!is.null(rowname_column)) {
        data[, rowname_column] <- rownames(data)
        # use as first column
        data <- data[, c(rowname_column, colnames(data)[1:ncol(data)-1])]
    }

    write.table(data, text.connection, sep = "\t", quote = F, row.names = F, col.names = T)
    text.result <- textConnectionValue(text.connection)
    close(text.connection)

    text.result <- paste0(text.result, collapse = "\n")

    return(text.result)
}

#' Add pathway names to the data.frame
#'
#' @param pathway_data data.frame with pathway values. Must contain a
#'        column 'Pathway' with the pathway identifiers
#' @param pathway_names list with pathway identifiers as key and the name
#'        as value
add_pathway_names <- function(pathway_data, pathway_names) {
    if (!"Pathway" %in% colnames(pathway_data)) {
        stop("Error: Missing required column 'Pathway' in pathway_data")
    }

    # convert the names to a named vector
    pathway_name_vector <- unlist(pathway_names)

    pathway_data$Name <- pathway_name_vector[pathway_data$Pathway]

    # use the name as the second column
    col_order <- colnames(pathway_data)[c(1, ncol(pathway_data), 2:(ncol(pathway_data)-1))]
    pathway_data <- pathway_data[, col_order]

    return(pathway_data)
}

#' Add average fold changes for pathway result
#'
#' @param result Result of the GSA algorithm as a data.frame
#' @param fold_changes Molecule level fold_changes
#' @param gene_index Index of genes mapping to a pathway
#' @param expression_data A data.frame containing the expression matrix
#' @return An updated version of result with an additional column 'av_foldchange'
add_pathway_foldchanges <- function(result, fold_changes, gene_index, expression_data) {
    # make sure the fold_changes contain the required columns
    for (required_col in c("Identifier", "logFC")) {
        if (!required_col %in% colnames(fold_changes)) {
            stop("Error: Missing ", required_col, " in fold_changes: ", paste(colnames(fold_changes)))
        }
    }

    # change the fold_changes to use the identifier as rownames
    rownames(fold_changes) <- fold_changes$Identifier

    # calculate the average fold-change for every pathway
    av_fc <- sapply(result$Pathway, function(pathway_id) {
        # get all genes for that pathway
        if (!pathway_id %in% names(gene_index)) {
            stop("Error: Missing index for ", pathway_id)
        }

        gene_ids <- rownames(expression_data)[gene_index[[pathway_id]]]

        this_av_fc <- mean(fold_changes[gene_ids, "logFC"], na.rm = T)

        # return the average logFC
        return(this_av_fc)
    })

    # add the average fold-changes as a new column
    result$av_foldchange <- as.numeric(av_fc)

    # add the direction if it isn't present
    if (!"Direction" %in% colnames(result)) {
        org_colnames <- colnames(result)
        result$Direction <- ifelse(result$av_foldchange > 0, "Up", "Down")

        result <- result[, c(org_colnames[1:2], "Direction", org_colnames[3:length(org_colnames)])]
    }

    return(result)
}
