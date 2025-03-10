#' Function to convert (for discrete values) and normalise (if set) the data.
#'
#' This function evaluates two global values as parameters: \code{edger.norm.function}
#' and \code{continuous.norm.function} which control the normalisation methods used
#' for discrete and continuous data respectively.
#'
#' @param expression.data data.frame containing the expression values
#' @param sample.data data.frame with all sample annotations (one sample per row)
#' @param design model.matrix specifying the experimental design
#' @param data.type type of data submitted only riboseq data is allowed in this case
prepareData <- function(expression.data, sample.data, design, analysis.group.1, analysis.group.2) {
 # Copy sample.data to the exp_de dataframe
    exp_de <- sample.data

    # Adds a column with the comparison groups, as defined in the design
    exp_de$Group <- NA
    exp_de$Group[design[, analysis.group.1] == 1] <- "c"
    exp_de$Group[design[, analysis.group.2] == 1] <- "d"

    # Checks that the number of samples matches in the data matches the one in the metadata
    if (nrow(exp_de) != ncol(expression.data)) {
        stop("Error Number of samples in submitted data (columns) do not match the metadata (rows in sample.data)")
    }
    # Checks that the metadata has the required info
    if (length(setdiff(c("SampleID", "SeqType", "SampleName"), colnames(exp_de))) > 0) {
        stop("Error: column SampleID, SampleName and/or SeqType missing from sample.data")
    }

    # It is important to ensure all columns are converted to factors!
    exp_de$Group <- factor(exp_de$Group, levels = c("c", "d"))
    exp_de$SeqType <- factor(exp_de$SeqType, levels = c("RNA", "RIBO"))

    # sample.groups is a global variable (string) defined in another script (needed for paired studies).
    # Checks if the column corresponding to the string "sample.groups" is among the sample.data columns
    # If so, it renames it to "Block" (used below to define paired-sample structure in analysis).

    if (exists("sample.groups") && sample.groups %in% colnames(sample.data)) {
        colnames(exp_de)[colnames(exp_de) == sample.groups] <- "Block"
        # This ensures the Block column is a factor.
        exp_de$Block <- factor(exp_de$Block)
    }


    # Checks if the expression.data submitted is integers or floats.
    are_columns_integer <- sapply(expression.data, function(x) {all(x == as.integer(x), na.rm = TRUE)})

    # If not all columns are integers, round and convert to integer (required for DeltaTE)
    if (!all(are_columns_integer)) {
        row_names <- rownames(expression.data)
        expression.data <- data.frame(lapply(expression.data, function(x) as.integer(round(x))))
        rownames(expression.data) <- row_names
    }

    # Changes the dataframe into a matrix
    expression.data <- as.matrix(expression.data)

    return(list(expression.data = expression.data, exp_de = exp_de))
}

#' "Public" function called by the analyzer to load the required libraries.
#' This is mainly put in a separate function to be able to quicky see if
#' required libraries are not available on the system.
#'
load_libraries <- function() {
    suppressPackageStartupMessages(library(edgeR))
    suppressPackageStartupMessages(library(limma))
    suppressPackageStartupMessages(library(DESeq2))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(plyr))
    suppressPackageStartupMessages(library(apeglm))
    suppressPackageStartupMessages(library(terapadog))
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

    # Prepare the data for DeltaTE (expression.data and exp_de)
    prepared_data <- prepareData(expression.data, sample.data, design, analysis.group.1, analysis.group.2)
    exp_de <- prepared_data$exp_de
    expression.data <- prepared_data$expression.data

    # Creates a contrast vector where samples are categorised between "c" and "d" for the comparison.
    padog_group <- exp_de$Group

    # Creates a filter for samples that are not in any comparison group: a vector with all SampleIDs were Group == NA.
    na_samples <- exp_de$SampleID[is.na(exp_de$Group)]

    # Remove the NA samples from exp_de
    exp_de <- exp_de[!exp_de$SampleID %in% na_samples, ]

    # Remove the NA samples from sample_data as well.
    sample.data <- sample.data[!sample.data$SampleID  %in% na_samples, ]

    # Remove NA samples from expression.data
    expression.data <- expression.data[, !colnames(expression.data) %in% na_samples, drop = FALSE]

    # Retrieves list of gene IDs from expression.data matching the ones in gene.indices
    gene_identifier_set <- lapply(gene.indices, function(gene_ids) {
        rownames(expression.data)[gene_ids]
    })

    # Sets the paired design off by default
    is_paired <- FALSE

    # The following is needed for PAIRED samples.
    # NB. It calls it "sample.groups" but is a different thing from the comparison groups ("c", "d").
    # This just defines the paired structure in the data submitted.
    if (exists("sample.groups")) {
        if (!sample.groups %in% colnames(sample.data)) {
            stop("Error: Failed to find defined sample.groups '", sample.groups, "' in the sample metadata. ",
                 "In the ReactomeGSA R package, this must also be specified as an 'additional_factor'")
        }
        # Updates to TRUE if a paired design has been given
        is_paired <- TRUE
    }

    # check padog's requirements and create nice error messages
    found_genes <- length(unique(as.numeric(unlist(gene.indices))))
    if (found_genes <= 10) {
        stop("Error: PADOG requires more than 10 proteins/genes to be found in the gene sets. Only ", found_genes,
        " identifiers be mapped to Reactome's pathways")
    }

    # there must be at least 6 samples
    if (any(table(padog_group) < 3)) {
        stop("Error: PADOG requires at least 3 samples per group.")
    }

    # implementation is below
    terapadog_result <- terapadog(
        esetm = as.matrix(expression.data),
        exp_de = exp_de,
        paired = is_paired,
        gslist = gene_identifier_set,
        NI = 1000, # number of iterations
        Nmin = 0) # minimum gene set size

    # convert to the required result output
    colnames(terapadog_result) <- plyr::mapvalues(
        from = c("ID", "Size", "meanAbsT0", "padog0", "PmeanAbsT", "Ppadog"),
        to = c("Pathway", "NGenes", "MeanAbsT0", "MeanWeightT0", "PValue", "FDR"),
        x = colnames(terapadog_result))

    return(terapadog_result[, c("Pathway", "FDR", "PValue", "NGenes", "MeanAbsT0", "MeanWeightT0")])

}


get_gene_fc <- function(expression.data, sample.data, design, data_type, analysis.group.1, analysis.group.2) {

    # prepare the data
    prepared_data <- prepareData(expression.data, sample.data, design, analysis.group.1, analysis.group.2)
    exp_de <- prepared_data$exp_de
    expression.data <- prepared_data$expression.data

    # DeltaTE cannot manage NA among the levels of a factor, so I need to remove them (code same as above)
    na_samples <- exp_de$SampleID[is.na(exp_de$Group)]
    # Remove the NA samples from exp_de
    exp_de <- exp_de[!exp_de$SampleID %in% na_samples, ]
    # Remove NA samples from expression.data
    expression.data <- expression.data[, !colnames(expression.data) %in% na_samples, drop = FALSE]

    # Similarly to above, pired
    paired <- FALSE

    if (exists("sample.groups")) {
        if (!sample.groups %in% colnames(sample.data)) {
            stop("Error: Failed to find defined sample.groups '", sample.groups, "' in the sample metadata. ",
                 "In the ReactomeGSA R package, this must also be specified as an 'additional_factor'")
        }
        paired <- TRUE
    }

    if (paired == TRUE) {
        design_TE <- ~ Block + Group + SeqType+ Group:SeqType
        design_R <- ~ Block + Group
    }
    else {
        design_TE <- ~ Group + SeqType+ Group:SeqType
        design_R <- ~ Group
    }

    # run DeltaTE to get the TE linear model
    ddsMat <- DESeqDataSetFromMatrix(
        countData = expression.data, # Ribo_counts and rna_counts should have been provided as a single dataframe already.
        colData = exp_de, # This is where we give DeltaTe the info
        design = design_TE
    )

    # Run DESeq
    ddsMat <- suppressMessages(DESeq(ddsMat))

    # Extract the fold changes
    res <- results(ddsMat, name = "Groupd.SeqTypeRIBO")

    # FC analysis for RNA counts:
    # Create filter using exp_de to find samples with SeqType "RNA"
    rna_samples <- exp_de %>%
      filter(SeqType == "RNA") %>%
      select(SampleID)

    # Extract Sample IDs as a vector
    rna_sample_ids <- rna_samples$SampleID

    # Select columns from expression.data
    expression.data.rna <- expression.data[, rna_sample_ids]

    # Execute the analysis for RNA counts
    ddsMat_rna <- DESeqDataSetFromMatrix(
        countData = expression.data.rna,
        colData=exp_de[which(exp_de$SeqType == "RNA"),],
        design = design_R)

    ddsMat_rna <- DESeq(ddsMat_rna)

    # Extract results from the ddsMat object
    res_rna <- results(ddsMat_rna, name="Group_d_vs_c")
    res_rna <- lfcShrink(ddsMat_rna,coef="Group_d_vs_c", res=res_rna)


    # FC analsysis for RIBO counts:
    # Create filter using exp_de to find samples with SeqType  "RIBO"
    ribo_data <- exp_de %>%
      filter(SeqType == "RIBO") %>%
      select(SampleID)

    # Extract Sample IDs as a vector
    ribo_sample_ids <- ribo_data$SampleID

    # Select columns from expression.data
    expression.data.ribo <- expression.data[, ribo_sample_ids]

    # Getting FC for Ribo
    ddsMat_ribo <- DESeqDataSetFromMatrix(
        countData = expression.data.ribo,
        colData = exp_de[which(exp_de$SeqType == "RIBO"),],
        design = design_R)

    ddsMat_ribo <- DESeq(ddsMat_ribo)

    # Extract results from the ddsMat object
    res_ribo <- results(ddsMat_ribo,name="Group_d_vs_c")
    res_ribo <- lfcShrink(ddsMat_ribo,coef="Group_d_vs_c", res=res_ribo)

    # Converts results for Translational Efficiency changes (TE) to dataframe
    res$Identifier <- rownames(res)
    res_combined <- as.data.frame(res)

    # Extracts info from RNA results (res_rna)
    res_rna_slim <- as.data.frame(res_rna) %>%
      select(RNA_FC = log2FoldChange, RNA_padj = padj)

    # Extracts info from RIBO results (res_ribo)
    res_ribo_slim <- as.data.frame(res_ribo) %>%
      select(RIBO_FC = log2FoldChange, RIBO_padj = padj)

    res_rna_slim$Identifier <- rownames(res_rna_slim)
    res_ribo_slim$Identifier <- rownames(res_ribo_slim)

    # Combines the results into a single dataframe
    res_combined <- res_combined %>%
      left_join(res_ribo_slim, by = "Identifier") %>%
      left_join(res_rna_slim, by = "Identifier")

    # Assignes the Regulatory Mode to each gene in the results df
    res_combined <- assign_Regmode(res_combined)

    # Rename and reorder columns
    names(res_combined)[names(res_combined) == "log2FoldChange"] <- "logFC"
    names(res_combined)[names(res_combined) == "pvalue"] <- "P.Value"
    names(res_combined)[names(res_combined) == "padj"] <- "adj.P.Val"
    col_order <- c("Identifier", colnames(res_combined)[1:ncol(res_combined)])

    return(res_combined[, col_order])
}
