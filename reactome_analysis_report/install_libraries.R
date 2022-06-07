# ------------------------------
# Script to install all required libraries
# in the docker container
# ------------------------------

# set the CRAN mirror to use
options(repos=structure(c(CRAN = 'http://cran.ma.imperial.ac.uk')))

# set the correct site library
.libPaths("/usr/local/lib/R/site-library")

install.packages("openxlsx")
install.packages("httr")
install.packages("progress")
install.packages("jsonlite")
# necessary as this is otherwise outdated
install.packages("vctrs")

# install the ReactomeGSA package
# install devtools if needed
if (!require(devtools)) {
  install.packages("devtools")
}

# install the ReactomeGSA package - version 1.3.7
devtools::install_github("reactome/ReactomeGSA")

# install the ReactomeGSA report package - version 1.0.7
devtools::install_github("reactome/ReactomeGSA.report")