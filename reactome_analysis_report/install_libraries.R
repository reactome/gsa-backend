# ------------------------------
# Script to install all required libraries
# in the docker container
# ------------------------------

# set the CRAN mirror to use
options(repos=structure(c(CRAN = 'http://cran.ma.imperial.ac.uk')))

install.packages("xlsx")
install.packages("httr")
install.packages("progress")
install.packages("jsonlite")

# install the ReactomeGSA package
# install devtools if needed
if (!require(devtools)) {
  install.packages("devtools")
}

# install the ReactomeGSA package
devtools::install_github("reactome/ReactomeGSA")

# install the ReactomeGSA report package - version 1.0.3
devtools::install_github("reactome/ReactomeGSA.report")