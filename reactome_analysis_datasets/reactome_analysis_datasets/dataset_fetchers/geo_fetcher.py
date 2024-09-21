import logging
import os
from typing import Tuple

import rpy2.robjects as ro
import rpy2.rinterface as ri
from rpy2.robjects import pandas2ri
from reactome_analysis_api.models.external_data_sample_metadata import ExternalDataSampleMetadata
from reactome_analysis_worker.analysers import ReactomeRAnalyser
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, \
    DatasetFetcherException

LOGGER = logging.getLogger(__name__)


class GeoFetcher(DatasetFetcher):
    """
    A DatasetFetcher to retrieve data from GEO
    """

    def __init__(self):
        # constructor of abstract super class
        super().__init__()

    def get_dataset_id(self, parameters: list) -> str:
        """
        Returns the dataset identifier
        :param parameters: A list of DatasetRequestParameter objects.
        :returns: The dataset identifier
        """
        return self._get_parameter(name="dataset_id", parameters=parameters)

    def load_dataset(self, parameters: list, reactome_mq) -> Tuple[str, ExternalData]:
        """
        Load the specified GEO dataset based on GSE id
        :param parameters: GSE id from GEO
        :returns empty string instead of count matrix, metadata of all datasets
        """
        identifier = self._get_parameter(name="dataset_id", parameters=parameters).strip()
        if not identifier:
            raise DatasetFetcherException(
                "Missing required parameter 'dataset_id' to load the example dataset.")
        if not identifier[0:3] == "GSE":
            raise DatasetFetcherException(
                f"{identifier} is not a valid geo_accession")
        
        # download the R GSE object
        self._update_status(progress=0.2, message="Downloading data from GEO...")
        gse_object = self._download_gse(identifier)

        self._update_status(progress=0.7, message="Extracting metadata...")
        metadata_obj = self._get_metadata(identifier, gse_object)

        # load the data
        self._update_status(progress=0.8, message="Extracting count matrix...")
        count_matrix = self._create_count_matrix(gse_object)

        return (count_matrix, metadata_obj)

    def _download_gse(self, gse_identifier: str):
        """Downloads the GSE object using the R package GEOquery

        :param gse_identifier: The identifier of the dataset to download
        :type gse_identifier: str
        """
        # activate R in Python
        pandas2ri.activate()
        # Load R GEOquery library for temporary R code
        ro.r('library(GEOquery)')
        ro.r(f'gse <- getGEO("{gse_identifier}", GSEMatrix = TRUE)')

        return ri.globalenv["gse"][0]
    
    def _get_metadata(self, gse_identifier: str, gse_object) -> ExternalData:
        """Extract the sample metadata from the GSE object

        :param gse_object: _description_
        :type gse_object: _type_
        :return: The metadata of the dataset as an ExternalData object
        :rtype: ExternalData
        """
        ri.globalenv["gse"] = gse_object

        ro.r("""
             # get the pheno data
             pheno_data <- pData(phenoData(gse))

             # remove all contact columns
             pheno_data <- pheno_data[, !grepl("contact_.*", colnames(pheno_data))]

             # remove all supplementary file columns
             pheno_data <- pheno_data[, !grepl("supplementary_file.*", colnames(pheno_data))]

             # remove keyword columns
             pheno_data <- pheno_data[, !grepl("[Kk]eywords.*", colnames(pheno_data))]

             # remove specific columns
             cols_to_remove <- c("geo_accession", "status", "submission_date", "last_update_date")
             pheno_data <- pheno_data[, !colnames(pheno_data) %in% cols_to_remove]

             # make some column names nicer
             colnames(pheno_data) <- gsub("_ch1", "", colnames(pheno_data))

             # add the sample_id as first column and remove the rownames
             org_colnames <- colnames(pheno_data)
             pheno_data$sample_id <- rownames(pheno_data)
             pheno_data <- pheno_data[, c("sample_id", org_colnames)]
             rownames(pheno_data) <- NULL

             # change the description columns
             descr_cols <- colnames(pheno_data)[grepl("description.*", colnames(pheno_data))]

             for (descr_col in descr_cols) {
                # get the real name
                nice_name <- gsub(":.*", "", pheno_data[1, descr_col])
             
                # replace with nice data
                pheno_data[, nice_name] <- gsub(".*:", "", pheno_data[, descr_col])
             
                # delete the old one
                pheno_data[, descr_col] <- NULL
             }

             # get the other metadata fields
             title <- experimentData(gse)@title
             geo_accession <- experimentData(gse)@other$geo_accession
             abstract <- abstract(gse)
             type <- experimentData(gse)@other$type

             # finally, save the colnames
             pheno_cols <- colnames(pheno_data)

             # save the sample ids
             sample_ids <- pheno_data$sample_id
        """)

        # sample-level metadata
        sample_metadata = list()

        for index, sample_field in enumerate(ri.globalenv["pheno_cols"]):
            field_metadata = {
                "name": sample_field,
                "values": [str(value) for value in ri.globalenv["pheno_data"][index]]
            }

            sample_metadata.append(field_metadata)

        # get the experiment type
        stored_type = ri.globalenv["type"][0]

        if stored_type == "Expression profiling by array":
            LOGGER.info("Dataset is microarray")
            experiment_type = "microarray_norm"
        elif stored_type == "Expression profiling by high throughput sequencing":
            LOGGER.info("Dataset is RNA-seq")
            experiment_type = "rnaseq_counts"
        else:
            LOGGER.warning(f"Unknown experiment type {stored_type}")
            # use microarray data as default
            experiment_type = "microarray_norm"

        metadata_obj = ExternalData(id=gse_identifier,
                                    title=ri.globalenv["title"][0],
                                    type=experiment_type,
                                    description=ri.globalenv["abstract"][0],
                                    group=None,
                                    sample_ids=[i for i in ri.globalenv["sample_ids"]],
                                    sample_metadata=sample_metadata,
                                    default_parameters=None)
        
        return metadata_obj


    def _create_count_matrix(self, gse_object) -> str:
        """
        Create count matrix based on the gse_object

        :param gse_object: The R gse object
        :returns: A string containing the count matrix as a tab delimited table.
        """
        ri.globalenv["gse"] = gse_object
        ro.r(f'count_matrix <- data.frame( exprs(gse) )')

        # Convert the R count_matrix to string with seperation
        count_matrix = ri.globalenv["count_matrix"]
        count_matrix_tsv = ReactomeRAnalyser.data_frame_to_string(count_matrix, add_rownames=True)
        return count_matrix_tsv
