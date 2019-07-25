# ------------------------------
# Script to install all required libraries
# in the docker container
# ------------------------------

# set the CRAN mirror to use
options(repos=structure(c(CRAN = 'http://cran.ma.imperial.ac.uk')))

required_libraries <- c("limma", "edgeR", "GSVA", "PADOG")

if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

if (!requireNamespace("plyr", quietly = TRUE))
	install.packages("plyr")

for (library_name in required_libraries) {
	BiocManager::install(library_name)
}

