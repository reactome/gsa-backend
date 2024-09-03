#' Function to convert (for discrete values) and normalise (if set) the data.
#'
#' This function evaluates two global values as parameters: \code{edger.norm.function}
#' and \code{continuous.norm.function} which control the normalisation methods used
#' for discrete and continuous data respectively.
#'
#' @param expression.data data.frame containing the expression values
#' @param sample.data data.frame with all sample annotations (one sample per row)
#' @param design model.matrix specifying the experimental design
#' @param data.type type of data submitted (ie. rnaseq, proteomics-sc, proteomics-int)
prepareData <- function(expression.data, sample.data, design, data.type) {
    if (data.type %in% c("rnaseq_counts", "proteomics_sc", "metabolomics")) {
        expression.data <- DGEList(counts=expression.data, samples = sample.data, group = sample.data$analysis_group)

        # normalize using edgeR's function
        expression.data <- calcNormFactors(expression.data, method = edger.norm.function)

        expression.data <- voom(expression.data, design=design, plot = FALSE)
    } else {
        if (continuous.norm.function == "none") {
            warning("Not performing normalization for ", data.type, " data")
            expression.data <- new("EList", list(E = expression.data, targets = sample.data))
        } else {
            # create the EList object
            expression.data <- new("EListRaw", list(E = expression.data, targets = sample.data))

            # normalize if set (may be set to "none")
            expression.data <- normalizeBetweenArrays(expression.data, method = continuous.norm.function)
        }
    }

    return(expression.data)
}

#' "Public" function called by the analyzer to load the required libraries.
#' This is mainly put in a separate function to be able to quicky see if
#' required libraries are not available on the system.
#'
load_libraries <- function() {
    suppressPackageStartupMessages(library(edgeR))
    suppressPackageStartupMessages(library(limma))
}

#' Main function to process the dataset
#'
#' @param expression.data data.frame containing the expression values
#' @param sample.data data.frame with all sample annotations (one sample per row)
#' @param design model.matrix specifying the experimental design
#' @param gene.indices a named list with each gene set as entry (and its id as name) and the index of the member
#'        genes based on the expression.data data.frame as values
#' @param data.type type of data submitted (ie. rnaseq, proteomics-sc, proteomics-int)
#' @param analysis.group.1 name of the first coefficient to test based on the experimental design. Must correspond
#'        to a colname in the \code{design}
#' @param analysis.group.2 name of the second coefficient to test based on the experimental design
process <- function(expression.data, sample.data, design, gene.indices, data.type, analysis.group.1, analysis.group.2) {
    # prepare the data
    expression.data <- prepareData(expression.data, sample.data, design, data.type)

    # create the contrast vector
    contrasts <- rep(0, ncol(design))
    contrasts[which(colnames(design) == analysis.group.1)] <- -1
    contrasts[which(colnames(design) == analysis.group.2)] <- 1

    # perform the analysis
    camera.result <- camera(y = expression.data, index = gene.indices, design = design, contrast = contrasts)

    # convert to the required result output
    gsa.result <- camera.result[, c("Direction", "FDR", "PValue", "NGenes")]
    gsa.result$Pathway <- rownames(gsa.result)

    return(gsa.result[, c("Pathway", "Direction", "FDR", "PValue", "NGenes")])
}

get_gene_fc <- function(expression.data, sample.data, design, data.type, analysis.group.1, analysis.group.2) {
    # prepare the data
    expression.data <- prepareData(expression.data, sample.data, design, data.type)

    # create the contrast vector
    contrasts <- rep(0, ncol(design))
    contrasts[which(colnames(design) == analysis.group.1)] <- -1
    contrasts[which(colnames(design) == analysis.group.2)] <- 1

    # create the fit
    fit <- lmFit(expression.data, design)
    cont.fit <- contrasts.fit(fit, contrasts)
    fit1 <- eBayes(cont.fit)

    result <- na.omit(topTable(fit1, adjust="fdr", number="all"))

    result$Identifier <- rownames(result)

    # move "Identifier" as first column
    col_order <- c("Identifier", colnames(result)[1:ncol(result)-1])

    return(result[, col_order])
}
