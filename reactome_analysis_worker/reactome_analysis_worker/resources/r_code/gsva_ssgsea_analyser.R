#' "Public" function called by the analyzer to load the required libraries.
#' This is mainly put in a separate function to be able to quicky see if
#' required libraries are not available on the system.
#'
load_libraries <- function() {
    suppressPackageStartupMessages(library(GSVA))
}

#' Main function to process the dataset
#'
#' @param expression.data data.frame containing the expression values
#' @param gene.indices a named list with each gene set as entry (and its id as name) and the index of the member
#'        genes based on the expression.data data.frame as values
#' @param data.type type of data submitted (ie. rnaseq, proteomics-sc, proteomics-int)
#' @param min.sz minimum pathway size to consider
#' @param max.sz maximum pathway size to consider
#' @param pathways A comma delimited list of pathways to include in the analysis
process <- function(expression.data, gene.indices, data.type, min.sz, max.sz, pathways) {
    kcdf <- "Gaussian"

    if (data.type %in% c("rnaseq_counts", "proteomics_sc")) {
        kcdf <- "Poisson"
    }

    # if pathways is set, limit the pathways to use for the comparison
    if (!is.null(pathways) && nchar(pathways) > 0) {
        pathways_to_include <- strsplit(pathways, ",")[[1]]
        # remove white-spaces from the pathway names
        pathways_to_include <- gsub(" ", "", pathways_to_include)

        # filter the gene set
        gene.indices <- gene.indices[names(gene.indices) %in% pathways_to_include]
    }

    # change gene names to indexes to match the pathway mapping result
    rownames(expression.data) <- 1:nrow(expression.data)

    # perform the GSVA analysis
    res <- gsva(as.matrix(expression.data), gene.indices, method="ssgsea", parallel.sz = 1,
                kcdf = kcdf, min.sz = min.sz, max.sz = max.sz, verbose = FALSE)

    res <- data.frame(res, stringsAsFactors = FALSE)

    # replace the rownames with a new column
    res$Pathway <- rownames(res)
    rownames(res) <- NULL

    res <- res[, c("Pathway", colnames(res)[1:ncol(res) - 1])]

    return(res)
}

get_gene_fc <- function(expression.data, sample.data, design, data.type, analysis.group.1, analysis.group.2) {
    # this is not supported
    return(data.frame())
}
