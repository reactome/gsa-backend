prepareData <- function(expression.data, sample.data, design, analysis.group.1, analysis.group.2) {
 # Copy sample.data to the exp_de dataframe
    exp_de <- sample.data

    # Adds a column with the comparison groups =, as defined in the design
    exp_de$Group <- NA
    exp_de$Group[design[, analysis.group.1] == 1] <- "c"
    exp_de$Group[design[, analysis.group.2] == 1] <- "d"

    # Debug Check: controls that nrow(exp_de) == ncol(expression.data)
    if (nrow(exp_de) != ncol(expression.data)) {
        stop("Error Number of samples in submitted data (columns) do not match the metadata (rows in sample.data)")
    }

    if (length(setdiff(c("SampleID", "SeqType", "SampleName"), colnames(exp_de))) > 0) {
        stop("Error: column SampleID and/or SeqType missing from sample.data")
    }

    # It is important to ensure all columns are converted to factors!
    exp_de$Group <- factor(exp_de$Group, levels = c("c", "d"))
    exp_de$SeqType <- factor(exp_de$SeqType, levels = c("RNA", "RIBO"))

    # sample.groups is a global variable (string) defined in another script (needed for paired studies).
    # Checks if the column corresponding to the string "sample.groups" is among the sample.data columns
    # If so, it renames it to "Block" (used below to define paired structure in analysis).

    if (exists("sample.groups") && sample.groups %in% colnames(sample.data)) {
        colnames(exp_de)[colnames(exp_de) == sample.groups] <- "Block"
        # This ensures the Block column is a factor. levels are not specified though. How to?
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


load_libraries <- function() {
    suppressPackageStartupMessages(library(edgeR))
    suppressPackageStartupMessages(library(limma))   # maybe not needed
    suppressPackageStartupMessages(library(DESeq2))
    suppressPackageStartupMessages(library(dplyr)) # Needed
    suppressPackageStartupMessages(library(plyr)) # Needed
    suppressPackageStartupMessages(library(apeglm)) # Needed
}


process <- function(expression.data, sample.data, design, gene.indices, data.type, analysis.group.1, analysis.group.2) {

    # Prepare the data for DeltaTE (expression.data and exp_de)
    prepared_data <- prepareData(expression.data, sample.data, design, analysis.group.1, analysis.group.2)
    exp_de <- prepared_data$exp_de
    expression.data <- prepared_data$expression.data

    # Creates a contrast vector (padog_group) where samples are categorised between "c" and "d" for the comparison.
    padog_group <- exp_de$Group

    # Creates a filter for samples that are not in any comparison group: a vector with all SampleIDs were Group == NA
    na_samples <- exp_de$SampleID[is.na(exp_de$Group)]

    # Remove the NA samples from exp_de
    exp_de <- exp_de[!exp_de$SampleID %in% na_samples, ]

    # Remove the NA samples from sample_data as well (same logic, here just to make sure to bring changes across both dfs.
    sample.data <- sample.data[!sample.data$SampleID  %in% na_samples, ]

    # Remove NA samples from expression.data
    expression.data <- expression.data[, !colnames(expression.data) %in% na_samples, drop = FALSE]

    # convert the gene.indices back to the identifiers (not sure what THIS does)
    gene_identifier_set <- lapply(gene.indices, function(gene_ids) {
        rownames(expression.data)[gene_ids]
    })

    # Sets the paired design off by default
    is_paired <- FALSE

    # The following is needed for PAIRED samples.
    # NB. It calls it "sample.groups" but is a different thing from the comparison groups ("c", "d").
    # This just defines the paired structure int he data submitted.
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

    # DeltaTE cannot manage NA among the levels of a factor, so I need to remove them! (code same as above)
    na_samples <- exp_de$SampleID[is.na(exp_de$Group)]
    # Remove the NA samples from exp_de
    exp_de <- exp_de[!exp_de$SampleID %in% na_samples, ]
    # Remove NA samples from expression.data
    expression.data <- expression.data[, !colnames(expression.data) %in% na_samples, drop = FALSE]


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
    res_rna <- results(ddsMat_rna, name="Group_d_vs_c") # Check name of result! we encode it as "d" and "c"
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
    res_ribo <- results(ddsMat_ribo,name="Group_d_vs_c")
    res_ribo <- lfcShrink(ddsMat_ribo,coef="Group_d_vs_c", res=res_ribo)

    # convert results to dataframe
    res$Identifier <- rownames(res)
    res_combined <- as.data.frame(res)

    res_rna_slim <- as.data.frame(res_rna) %>%
      select(RNA_FC = log2FoldChange, RNA_padj = padj)

    res_ribo_slim <- as.data.frame(res_ribo) %>%
      select(RIBO_FC = log2FoldChange, RIBO_padj = padj)

    res_rna_slim$Identifier <- rownames(res_rna_slim)
    res_ribo_slim$Identifier <- rownames(res_ribo_slim)

    res_combined <- res_combined %>%
      left_join(res_ribo_slim, by = "Identifier") %>%
      left_join(res_rna_slim, by = "Identifier")

    res_combined <- assign_Regmode(res_combined)

    # Rename and reorder columns
    names(res_combined)[names(res_combined) == "log2FoldChange"] <- "logFC"
    names(res_combined)[names(res_combined) == "pvalue"] <- "P.Value"
    names(res_combined)[names(res_combined) == "padj"] <- "adj.P.Val"
    col_order <- c("Identifier", colnames(res_combined)[1:ncol(res_combined)])

    return(res_combined[, col_order])
}


# Implementation of TerraPadog algorithm


terapadog <- function (esetm = NULL, exp_de = NULL, paired = FALSE,
    gslist = "KEGGRESTpathway", organism = "hsa", annotation = NULL,
    gs.names = NULL, NI = 1000, Nmin = 3, verbose = TRUE, parallel = FALSE, dseed = NULL,
    ncr = NULL) {

    # validity checks of the data
    # Initial checks on the data (as PADOG would do). Some have been modified/removed to fit new kind of data.
    if (length(gslist) == 1 && gslist == "KEGGRESTpathway") {
        stopifnot(nchar(organism) == 3)
        res <- keggLink("pathway", organism)
        a <- data.frame(path = gsub(paste0("path:", organism),
                                    "", res), gns = gsub(paste0(organism, ":"),
                                                         "", names(res)))
        gslist <- tapply(a$gns, a$path, function(x) {
            as.character(x)
        })
        gs.names <- keggList("pathway", organism)[paste0("path:",
                                                        organism, names(gslist))]
        names(gs.names) <- names(gslist)
        stopifnot(length(gslist) >= 3)
        rm(res, a)
    }
    stopifnot(is.matrix(esetm))
    stopifnot(all(dim(esetm) > 4))
    stopifnot(mode(gslist) == "list")
    stopifnot(length(gslist) >= 3)

    if (!is.null(gs.names)) {
        stopifnot(length(gslist) == length(gs.names))
    }
    stopifnot(class(NI) == "numeric")
    stopifnot(NI > 5)
    if (!is.null(annotation)) {
        if (!annotation %in% c("hgu133a.db", "hgu133plus2.db")) {
            stopifnot(require(annotation, character.only = TRUE))
        }
        stopifnot(sum(rownames(esetm) %in% mappedkeys(get(paste0(substr(annotation,
            1, nchar(annotation) - 3), "ENTREZID")))) > 4)
    }
    else {
        stopifnot(sum(rownames(esetm) %in% as.character(unlist(gslist))) >
            10 & !any(duplicated(rownames(esetm))))
    }

    # Initial operations on the gene set lists. Unchanged from PADOG.
    gf <- table(unlist(gslist))
    if (!(var(gf) == 0)) {
        if (quantile(gf, 0.99) > mean(gf) + 3 * sd(gf)) {
            gf[gf > quantile(gf, 0.99)] <- quantile(gf, 0.99)
        }
        gff <- function(x) {
            1 + ((max(x) - x)/(max(x) - min(x)))^0.5
        }
        gf <- gff(gf)
    }
    else {
        fdfd <- unique(unlist(gslist))
        gf <- rep(1, length(fdfd))
        names(gf) <- fdfd
    }
    allGallP <- unique(unlist(gslist))

    restg <- setdiff(rownames(esetm), names(gf))
    appendd <- rep(1, length(restg))
    names(appendd) <- restg
    gf <- c(gf, appendd)
    stopifnot(all(!duplicated(rownames(esetm))))
    stopifnot(sum(rownames(esetm) %in% allGallP) > 10)
    if (verbose) {
        cat(paste0("Starting with ", length(gslist), " gene sets!"))
        cat("\n")
    }
    gslist <- gslist[unlist(lapply(gslist, function(x) {
        length(intersect(rownames(esetm), x)) >= Nmin
    }))]
    gs.names <- gs.names[names(gslist)]
    stopifnot(length(gslist) >= 3)
    if (verbose) {
        cat(paste0("Analyzing ", length(gslist), " gene sets with ", Nmin, " or more genes!"))
        cat("\n")
    }

    if (!is.null(dseed))
        set.seed(dseed)

    # Here starts what was modified more heavily.

    # This section groups matching RNA and RIBO samples under a single index
    # (to keep the pair together during permutations).

    # Orders matching RIBO and RNA samples by the SampleName value (which is provided by the user when submitting data)
    exp_de_ordered <- exp_de[do.call(order, exp_de["SampleName"]),]

    # Retrieves the indexes (from exp_de) for each RNA and matching RIBO Sample.
    ordered_indexes <- rownames(exp_de_ordered)
    grouped_indexes <- as.data.frame(matrix(ordered_indexes, ncol = 2, byrow = TRUE))

    # Adds to the dataframe info on the Group for each couple of indices
    grouped_indexes$Group <- exp_de_ordered[exp_de_ordered$SeqType != "RIBO", ]$Group

    block <- NULL
    if (paired) {
        # Retrieves the info for paired samples from the dataframe
        grouped_indexes$Block <- exp_de_ordered[exp_de_ordered$SeqType != "RIBO", ]$Block
        block <- factor(grouped_indexes$Block)
    }

    # Extracts the group vector from the dataframe, to use it for permutations as PADOG does
    group <- grouped_indexes$Group

    # Global variables

    G <- factor(group) # creates a vector of factors ("c" or "d" in riboNAV)
    Glen <- length(G) # number of samples (since each factor is tied to a sample)
    tab <- table(G) # counts frequency for each factor
    idx <- which.min(tab) # finds factors with least number of samples assigned
    minG <- names(tab)[idx] # retrieve its name
    minGSZ <- tab[idx] # retrieves number of samples assigned to it
    bigG <- rep(setdiff(levels(G), minG), length(G))

    # topSigNum <- dim(esetm)[1]

    # Creates all the possible "group" label permutations for the given data.
    combFun <- function(gi, countn = TRUE) {
        g <- G[gi]
        tab <- table(g)
        if (countn) {
            minsz <- min(tab)
            ifelse(minsz > 10, -1, choose(length(g), minsz))
        }
        else {
            dup <- which(g == minG)
            cms <- combn(length(g), tab[minG])
            del <- apply(cms, 2, setequal, dup)
            if (paired) {
                cms <- cms[, order(del, decreasing = TRUE), drop = FALSE]
                cms[] <- gi[c(cms)]
                cms
            }
            else {
                cms[, !del, drop = FALSE]
            }
        }
    }


    if (paired) {
        bct <- tapply(seq_along(G), block, combFun, simplify = TRUE)
        nperm <- ifelse(any(bct < 0), -1, prod(bct))
        if (nperm < 0 || nperm > NI) {
            btab <- tapply(seq_along(G), block, `[`, simplify = FALSE)
            bSamp <- function(gi) {
                g <- G[gi]
                tab <- table(g)
                bsz <- length(g)
                minsz <- tab[minG]
                cms <- do.call(cbind, replicate(NI, sample.int(bsz,
                  minsz), simplify = FALSE))
                cms[] <- gi[c(cms)]
                cms
            }
            combidx <- do.call(rbind, lapply(btab, bSamp))
        }
        else {
            bcomb <- tapply(seq_along(G), block, combFun, countn = FALSE,
                simplify = FALSE)
            colb <- expand.grid(lapply(bcomb, function(x) 1:ncol(x)))[-1,
                , drop = FALSE]
            combidx <- mapply(function(x, y) x[, y, drop = FALSE],
                bcomb, colb, SIMPLIFY = FALSE)
            combidx <- do.call(rbind, combidx)
        }
    }
    else {
        # Checks number of permutations that would be created. Must be >0 and < 1000
        nperm <- combFun(seq_along(G))
        if (nperm < 0 || nperm > NI) {
            combidx <- do.call(cbind, replicate(NI, sample.int(Glen,
                minGSZ), simplify = FALSE))
        } else {
            # If nperm is acceptable, then permutation matrix combidx is generated (is the product of "cms" of combFun)
            combidx <- combFun(seq_along(G), countn = FALSE)
        }
    }

    NI <- ncol(combidx)
    if (verbose) {
        cat("# of permutations used:", NI, "\n")
    }

    deINgs <- intersect(rownames(esetm), unlist(gslist))
    gslistINesetm <- lapply(gslist, match, table = deINgs, nomatch = 0)
    MSabsT <- MSTop <- matrix(NA, length(gslistINesetm), NI +
        1)

    gsScoreFun <- function(G, block) {
        force(G) # Forces evaluation of G
        force(block) # Forces evaluation of block

        # This bit re-assigns the "group" label to the samples according to the combidx matrix, if past the first iteration

        if (ite > 1) {
            G <- bigG # Need to check the original PADOG for bigG
            G[combidx[, ite - 1]] <- minG
            G <- factor(G)

            # Unpacks grouped_indexes, so to have "group" labels for each sample (RNA count and RIBO count) in a vector
            for (i in seq_along(G)) {

                # Get the new label from G
                value <- as.character(G[i])

                # Updates exp_de according to the retrieved indexes (from grouped_indexes) with the new "group" label
                exp_de[grouped_indexes[i,]$V1, ]$Group <- value
                exp_de[grouped_indexes[i,]$V2, ]$Group <- value
            }
        }

        if (paired) {
            design <- ~ Block + Group + SeqType+ Group:SeqType # Paired designs do not have batch effect correction!
        }
        else { # Add if statement, if there is batch, then do this, else do not
            design <- ~ Group + SeqType+ Group:SeqType
        }

        # Setup the ddsMat object
        ddsMat <- DESeqDataSetFromMatrix(
        countData = esetm,
        colData = exp_de,
        design = design
        )

        # Calculate results (without printing messages on console)
        ddsMat <- suppressMessages(DESeq(ddsMat))

        # Extract specific comparison of interest
        res <- results(ddsMat, name = "Groupd.SeqTypeRIBO")

        # originally, moderated t-values (abs value) were retrived and stored in a df.
        # Now it just extracts adjusted p-values.
        de <- res$padj
        names(de) <- rownames(res)

        # The scaling happens.
        degf <- scale(cbind(de, de * gf[names(de)]))
        rownames(degf) <- names(de)
        degf <- degf[deINgs, , drop = FALSE]

        sapply(gslistINesetm, function(z) {
            X <- na.omit(degf[z, , drop = FALSE])
            colMeans(X, na.rm = TRUE) * sqrt(nrow(X))
        })
    }

    if (parallel && requireNamespace("doParallel", quietly = TRUE) &&
        requireNamespace("parallel", quietly = TRUE)) {
        ncores <- parallel::detectCores()
        if (!is.null(ncr))
            ncores <- min(ncores, ncr)
        if (verbose) {
            clust <- parallel::makeCluster(ncores, outfile = "")
        }
        else {
            clust <- parallel::makeCluster(ncores)
        }
        doParallel::registerDoParallel(clust)
        tryCatch({
            parRes <- foreach(ite = 1:(NI + 1), .combine = "c",
                .packages = "DESeq2") %dorng% { # original: .packages = "limma"
                Sres <- gsScoreFun(G, block)
                tmp <- list(t(Sres))
                names(tmp) <- ite
                if (verbose && (ite%%10 == 0)) {
                  cat(ite, "/", NI, "\n")
                }
                tmp
            }
            parRes <- do.call(cbind, parRes[order(as.numeric(names(parRes)))])
            evenCol <- (1:ncol(parRes))%%2 == 0
            MSabsT[] <- parRes[, !evenCol]
            MSTop[] <- parRes[, evenCol]
            rm(parRes)
        }, finally = parallel::stopCluster(clust))
    }
    else {
        if (parallel)
            message("Execute in serial! Packages 'doParallel' and 'parallel'\n needed for parallelization!")
        for (ite in 1:(NI + 1)) {
            Sres <- gsScoreFun(G, block)
            MSabsT[, ite] <- Sres[1, ]
            MSTop[, ite] <- Sres[2, ]
            if (verbose && (ite%%10 == 0)) {
                cat(ite, "/", NI, "\n")
            }
        }
    }
    meanAbsT0 <- MSabsT[, 1]
    padog0 <- MSTop[, 1]
    MSabsT <- scale(MSabsT)
    MSTop <- scale(MSTop)

    # mff is what checks the real score against the iterative ones.
    # Originally, PADOG compares t-scores (the bigger, the more significant).
    # Since we are using p-value, the smaller the more significant, so the comparison sign was changed.

    mff <- function(x) {
        mean(x[-1] < x[1], na.rm = TRUE) # Original mean(x[-1] > x[1], na.rm = TRUE)
    }
    PSabsT <- apply(MSabsT, 1, mff)
    PSTop <- apply(MSTop, 1, mff)
    PSabsT[PSabsT == 0] <- 1/NI/100
    PSTop[PSTop == 0] <- 1/NI/100

    # Removed part for plots
    if (!is.null(gs.names)) {
        myn <- gs.names
    }
    else {
        myn <- names(gslist)
    }
    SIZE <- unlist(lapply(gslist, function(x) {
        length(intersect(rownames(esetm), x))
    }))
    res <- data.frame(Name = myn, ID = names(gslist), Size = SIZE,
        meanAbsT0, padog0, PmeanAbsT = PSabsT, Ppadog = PSTop,
        stringsAsFactors = FALSE)
    ord <- order(res$Ppadog, -res$padog0)
    res <- res[ord, ]
    return(res)
}


assign_Regmode <- function(res_df) {
    # Function checks padj (the padj value for the TE change), the RIBO_padj (same but for the FC from Ribo-Seq counts)
    # and the RNA_padj (same as above, but for RNA FC).
    res_df$RegMode <- NA
    res_df$RegModeExplicit <- NA


    #Applying the "Forwarded" RegMode (genes regulated by mRNA abundance, with no signficant changes in Translational Efficiency)
    forwarded <- res_df$padj > 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj < 0.05
    res_df$RegMode[forwarded] <- "Forwarded"
    # Checks directionality of FC and assigns a more "explicit" RegMode value
    up_forwarded <- res_df$RegMode == "Forwarded" & res_df$RNA_FC > 0
    res_df$RegModeExplicit[up_forwarded] <- "(Up)regulated, driven by mRNA transcription only"
    down_forwarded <- res_df$RegMode == "Forwarded" & res_df$RNA_FC < 0
    res_df$RegModeExplicit[down_forwarded] <- "(Down)regulated, driven by mRNA transcription only"


    # Applying "Exclusive" RegMode (genes regulated by translation effiency, with no change in mRNA abundance.
    exclusive <- res_df$padj < 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj > 0.05
    res_df$RegMode[exclusive] <- "Exclusive"
    # Checks directionality of FC and assigns a more "explicit" RegMode value
    up_exclusive <- res_df$RegMode == "Exclusive" & res_df$log2FoldChange > 0
    res_df$RegModeExplicit[up_exclusive] <- "(Up)regulated, driven by mRNA translation only"
    down_exclusive <- res_df$RegMode == "Exclusive" & res_df$log2FoldChange < 0
    res_df$RegModeExplicit[down_exclusive] <- "(Down)regulated, driven by mRNA translation only"


    # Applying "Buffered" RegMode. This requires a more complex condition. All padjs must be significant, but
    # The directionality of the Fold change between TE and RNA must be opposite.
    buffered <- (res_df$padj < 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj < 0.05) &
      res_df$log2FoldChange * res_df$RNA_FC < 0
    res_df$RegMode[buffered] <- "Buffered"
    # Buffered (special case, when the effect of TE and RNA cancels out the change in RIBO)
    buffered_special <- res_df$padj < 0.05 & res_df$RIBO_padj > 0.05 & res_df$RNA_padj < 0.05
    res_df$RegMode[buffered_special] <- "Buffered"
    # Checks directionality of FC and assigns a more "explicit" RegMode value
    buffered_mRNA_down <- res_df$RegMode == "Buffered" & res_df$RNA_FC < 0
    res_df$RegModeExplicit[buffered_mRNA_down] <- "Buffered, decrease in mRNA levels counteracted by increase in translation"
    buffered_mRNA_up <- res_df$RegMode == "Buffered" & res_df$RNA_FC > 0
    res_df$RegModeExplicit[buffered_mRNA_up] <- "Buffered, increase in mRNA levels counteracted by decrease in translation"


    # Applying "Intensified" RegMode. This requires a more complex condition. All padjs must be significant, but
    # The directionality of the Fold change between TE and RNA must be the same.
    intensified <- (res_df$padj < 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj < 0.05) &
      res_df$log2FoldChange * res_df$RNA_FC > 0
    res_df$RegMode[intensified] <- "Intensified"
    # Checks directionality of FC and assigns a more "explicit" RegMode value
    intensified_down <- res_df$RegMode == "Intensified" & res_df$RNA_FC < 0
    res_df$RegModeExplicit[intensified_down] <- "Synergic decrease in both transcription and translation"
    insensified_up <- res_df$RegMode == "Intensified" & res_df$RNA_FC > 0
    res_df$RegModeExplicit[insensified_up] <- "Synergic increase in both transcription and translation"


    # No change case
    no_change <- res_df$padj > 0.05 & res_df$RIBO_padj > 0.05 & res_df$RNA_padj > 0.05
    res_df$RegMode[no_change] <- "No Change"
    res_df$RegModeExplicit[no_change] <- "No significant change detected in transcription or translation"


    # When a padj value is NA -> A padj is not calculated if the gene data is below a quality treshold enforced by DESeq2.
    no_data <- is.na(res_df$padj) | is.na(res_df$RIBO_padj) | is.na(res_df$RNA_padj)
    res_df$RegMode[no_data] <- "Undeterminable"
    res_df$RegModeExplicit[no_data] <- "One or more adjusted p-values are missing (NA)"


    # Combination of padjs not categorised by DeltaTE (Defined as "Undetermined" by Chotani et al. 2019)
    undetermined <- (res_df$padj > 0.05 & res_df$RIBO_padj > 0.05 & res_df$RNA_padj < 0.05)  |
        (res_df$padj > 0.05 & res_df$RIBO_padj < 0.05 & res_df$RNA_padj > 0.05) |
        (res_df$padj < 0.05 & res_df$RIBO_padj > 0.05 & res_df$RNA_padj > 0.05)
    res_df$RegMode[undetermined] <- "Undetermined"
    res_df$RegModeExplicit[undetermined] <- "Cannot be assigned to any Regulatory Mode"


    return(res_df)
}