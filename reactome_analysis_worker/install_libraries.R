# ------------------------------
# Script to install all required libraries
# in the docker container
# ------------------------------

# set the CRAN mirror to use
options(repos=structure(c(CRAN = 'http://cran.ma.imperial.ac.uk')))

required_libraries <- c("limma", "edgeR", "GSVA", "PADOG", "DESeq2", "apeglm", "terapadog")

if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

if (!requireNamespace("plyr", quietly = TRUE))
	install.packages("plyr")

if (!requireNamespace("dplyr", quietly = TRUE))
	install.packages("dplyr")

for (library_name in required_libraries) {a
	BiocManager::install(library_name)
}

