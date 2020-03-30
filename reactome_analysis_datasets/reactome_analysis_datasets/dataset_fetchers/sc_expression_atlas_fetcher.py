"""
Class to fetch datasets from the
Single Cell Expression Atlas resource
"""

import urllib3
import logging
import tempfile
import csv
import os
import zipfile
import scipy.io
import numpy
import io
from reactome_analysis_utils import reactome_mq
from reactome_analysis_datasets.dataset_fetchers import expression_atlas_fetcher
from reactome_analysis_datasets.dataset_fetchers.abstract_dataset_fetcher import DatasetFetcher, ExternalData, DatasetFetcherException


logger = logging.getLogger(__name__)


class ScExpressionAtlasFetcher(DatasetFetcher):
    """
    A DatasetFetcher to retrieve experiments from EBI's
    Single Cell Expression Atlas resource
    """
    def get_dataset_id(self, parameters: list) -> str:
        """
        Returns the dataset identifier. This identifier is the
        combination of the dataset id and "k"
        :param parameters: A list of DatasetRequestParameter objects.
        :returns: The dataset identifier
        """
        dataset_id = self._get_parameter("dataset_id", parameters)
        k = self._get_parameter("k", parameters)

        if not dataset_id or not k:
            return None

        return "{}_{}".format(dataset_id, k)

    def load_dataset(self, parameters: list, reactome_mq: reactome_mq.ReactomeMQ) -> (str, ExternalData):
        """
        Load the specified Single Cell ExpressionAtlas experiment
        :param parameters: A list of DatasetRequestParameter objects
        :param reactome_mq: The MQ used to process messages.
        :returns: (data, summary)
        """
        # get the id and k parameter
        identifier = self._get_parameter(name="dataset_id", parameters=parameters)
        k = int(self._get_parameter(name="k", parameters=parameters))

        if not identifier:
            raise DatasetFetcherException("Missing required parameter 'dataset_id' to load the Single Cell Expression Atlas dataset")
        if not k:
            raise DatasetFetcherException("Missing required parameter 'k' to load the Single Cell Expression Atlas dataset")

        # get the clustering data
        cell_clusterings = self._get_cell_clusterings(dataset_id=identifier, k=k)

        # download the matrix files
        logger.debug("Downloading norm counts for {id}".format(id=identifier))
        file_url = "https://www.ebi.ac.uk/gxa/sc/experiment/{id}/download/zip?fileType=normalised&accessKey=".format(id=identifier)
        matrix_file_dir = self._download_zip_file(file_url=file_url)

        # get the average expression per cluster
        logger.debug("Calculating average expression per cluster...")
        (av_exp, rownames, colnames) = self._get_av_cluster_expression(matrix_file_dir=matrix_file_dir, cell_clusterings=cell_clusterings)

        # create the tab-delimited string to represent the expression table
        logger.debug("Creating expression table...")
        expression_table = self._create_expression_table(expression=av_exp, rownames=rownames, colnames=colnames)

        # create a sensible summary object
        logger.debug("Creating summary object...")
        summary = self._create_summary(dataset_id=identifier, k=k, sample_ids=colnames, cell_clusterings=cell_clusterings)

        # return the data
        return (expression_table, summary)

    def _get_cell_clusterings(self, dataset_id: str, k: int) -> dict:
        """
        Retrieves the cell clustering data for the
        specific dataset. Extracts the clustering result
        for the specified `k` value.
        :param dataset_id: The dataset's id
        :param k: K used for the clustering
        :returns: A dict with the cluster id as key and the cell ids as values in a list
        """
        file_url = "https://www.ebi.ac.uk/gxa/sc/experiment/{id}/download?fileType=cluster&accessKey=".format(id=dataset_id)

        # download the clustering tsv
        clustering_tsv_string = self._download_file(file_url)

        # process the tsv data
        reader = csv.DictReader(clustering_tsv_string.decode().split("\n"), delimiter='\t')

        for row in reader:
            if not "K" in row:
                logger.error("Failed to find K in cluster file for {}".format(dataset_id))
                raise DatasetFetcherException("Invalid cell clustering result retrieved for dataset {}".format(dataset_id))

            # process the clusterings if K matches
            if int(row["K"]) == k:
                result_dict = dict()

                for cell_id in row:
                    # ignore K and sel.K
                    if cell_id == "K" or cell_id == "sel.K":
                        continue

                    # change the purely numeric ids to "Cluster #"
                    cluster_id = "Cluster {num_id}".format(num_id=row[cell_id])

                    if not cluster_id in result_dict:
                        result_dict[cluster_id] = [cell_id]
                    else:
                        result_dict[cluster_id].append(cell_id)

                # return the result
                return result_dict

        raise DatasetFetcherException("K = {k} does not exist for dataset {dataset_id}".format(
                                      k=str(k), dataset_id=dataset_id))

    def _download_file(self, file_url: str):
        """
        Downloads the specified file / page and returns
        its content.
        :param file_url: The file's / page's url
        :returns: The content as a binary
        """
        # download the file
        http = urllib3.PoolManager()
        download_request = http.request("GET", file_url)

        if download_request.status == 404:
            logger.info("Failed to find file {}".format(file_url))
            raise DatasetFetcherException("Unknown Single Cell Expression Atlas dataset")

        if download_request.status != 200:
            logger.error("Failed to download ZIP file from scExpressionAtlas({}): {}"
                .format(str(download_request.status), file_url))

        return download_request.data

    def _download_zip_file(self, file_url: str) -> tempfile.TemporaryDirectory:
        """
        Downloads the ZIP file containing The ZIP file is automatically
        extracted and the path to the (temporary) directory
        returned.
        :param dataset_id: The dataset's id
        :returns: The path to the newly created directory containing the files
        """
        zip_file_content = self._download_file(file_url=file_url)

        # create the temporary directory
        try:
            tmp_dir = tempfile.TemporaryDirectory()

            # store the zip file there
            logger.debug("Storing zip file at " + tmp_dir.name)
            zip_file_path = os.path.join(tmp_dir.name, "counts.zip")
            with open(zip_file_path, "wb") as writer:
                writer.write(zip_file_content)

            # free the memory
            del(zip_file_content)
        
            logger.debug("Extracting ZIP file content")
            zip_file = zipfile.ZipFile(file=zip_file_path)
            zip_file.extractall(path=tmp_dir.name)

            # delete the zip file again
            os.remove(zip_file_path)

            return tmp_dir
        except Exception as e:
            if tmp_dir:
                tmp_dir.cleanup()

            logger.error("Failed to process ZIP file for {}: {}".format(file_url, str(e)))
            raise DatasetFetcherException("Failed to retrieve data for {}".format(file_url))

    def _get_av_cluster_expression(self, matrix_file_dir: tempfile.TemporaryDirectory, cell_clusterings: dict) -> tuple:
        """
        Calculates the average gene expression per cluster.

        :param matrix_file_dir: The temporary directory containing the three matrix files.
        :param cell_clusterings: A dict with the cluster ids as key and the cells as values of a list.
        :returns: A tuple with (numpy.ndarray with expression values, rownames, columnnames)
        """
        # get the required filenames
        mtx_file = None
        col_file = None
        row_file = None

        for filename in os.listdir(matrix_file_dir.name):
            if filename[-4:] == ".mtx":
                mtx_file = os.path.join(matrix_file_dir.name, filename)
            if filename[-9:] == ".mtx_cols":
                col_file = os.path.join(matrix_file_dir.name, filename)
            if filename[-9:] == ".mtx_rows":
                row_file = os.path.join(matrix_file_dir.name, filename)

        if not mtx_file or not col_file or not row_file:
            raise DatasetFetcherException("Failed to retrieve required matrix files")

        # load the matrix
        matrix_data = scipy.io.mmread(mtx_file).tocsr()

        # load the colnames
        colnames = list()
        with open(col_file, "r") as reader:
            for line in reader:
                colname = line.strip()

                if len(colname) > 0:
                    colnames.append(colname)

        # load the rownames
        rownames = list()
        with open(row_file, "r") as reader:
            for line in reader:
                rowname = line.strip()

                if len(rowname) > 0:
                    rownames.append(rowname)

        # process every cluster
        cluster_ids = list()
        cluster_av_expressions = list()

        for cluster_id in cell_clusterings:
            # get the cols for this cluster
            cell_ids = cell_clusterings[cluster_id]
            cluster_cols = list()

            for col_index in range(0, len(colnames)):
                if colnames[col_index] in cell_ids:
                    cluster_cols.append(col_index)

            # slice the matrix
            cluster_col_data = matrix_data[:,cluster_cols]

            # get the average expression
            av_cluster_expression = cluster_col_data.mean(axis=1)

            # save the data
            cluster_ids.append(cluster_id)
            cluster_av_expressions.append(av_cluster_expression.flatten().tolist()[0])

        # create the numpy array
        av_exp_array = numpy.column_stack(cluster_av_expressions)
    
        if len(av_exp_array.shape) != 2:
            logger.error("Invalid shape of numpy array: " + str(av_exp_array.shape))
            raise DatasetFetcherException("Failed to extract read count data")

        if av_exp_array.shape[0] != len(rownames):
            logger.error("Different number of rows ({len_row}) and rownames ({len_names})".format(
                len_row=str(av_exp_array[0], len_names=str(len(rownames)))
            ))
            raise DatasetFetcherException("Failed to extract read count data")

        if av_exp_array.shape[1] != len(cluster_ids):
            logger.error("Different number of columns ({len_col}) and colnames ({len_names})".format(
                len_col=str(av_exp_array[1]), len_names=str(len(colnames))
            ))
            raise DatasetFetcherException("Failed to extract read count data")

        return (av_exp_array, rownames, cluster_ids)

    def _create_expression_table(self, expression: numpy.ndarray, rownames: list, colnames: list) -> str:
        """
        Creates a tab-delimited table representing the expression values, as well
        as the matching column and row names.
        :param expression: The expression values as a ndarray
        :param rownames: A list with rownames
        :param colnames: A list with colnames
        """
        # make sure everything fits
        if len(expression.shape) != 2 or expression.shape[0] != len(rownames) or expression.shape[1] != len(colnames):
            logger.error("Could not create expression table. Rownames, colnames and expression values do not match.")
            raise DatasetFetcherException("Failed to convert average cell counts")

        # first line it he colnames
        lines = ["\t" + "\t".join(colnames)]

        # add the rows
        for row_index in range(0, expression.shape[0]):
            string_list = [str(value) for value in expression[row_index, :].tolist()]
            lines.append(rownames[row_index] + "\t" + "\t".join(string_list))

        # return the final string
        return "\n".join(lines)

    def _create_summary(self, dataset_id: str, k: int, sample_ids: list, cell_clusterings: dict) -> ExternalData:
        """
        Creates an ExternalData object describing the datasets
        metadata.
        :param dataset_id: The dataset's id
        :param k: As defined in the k-means clustering
        :param sample_ids: A list containing all sample ids in the dataset.
        :param cell_clusterings: A dict with the cluster ids as key and the cells as values of a list.
        :returns: An ExternalData object
        """
        # get the metadata file
        file_url = "https://www.ebi.ac.uk/gxa/sc/experiment/{id}/download/zip?fileType=experiment-metadata&accessKey=".format(id=dataset_id)
        file_directory = self._download_zip_file(file_url=file_url)

        # get the idf and sdrf file
        idf_file = None
        sdrf_file = None

        for filename in os.listdir(file_directory.name):
            if filename[-7:] == "idf.txt":
                idf_file = os.path.join(file_directory.name, filename)
            if filename[-8:] == "sdrf.txt":
                sdrf_file = os.path.join(file_directory.name, filename)

        # Create the dict to create the ExternalData object from
        summary = {"type": "rnaseq_counts", 
                   "id": "{id}_k{k}".format(id=dataset_id, k=str(k)),
                   "title": "Single Cell Expression Atlas dataset {id} clustered at k = {k}".format(id=dataset_id, k=str(k)),
                   "description": "Enternal dataset loaded from Single Cell Expression Atlas",
                   "sample_ids": sample_ids,
                   "default_parameters": [{"name": "discrete_norm_function", "value": "TMM"}]}
        
        # add the publication information
        if idf_file:
            with open(idf_file, "r") as reader:
                for line in reader:
                    if line.startswith("Investigation Title"):
                        summary["title"] = line[line.find("\t"):].strip() + " clustered at k = " + str(k)
                    if line.startswith("Publication DOI"):
                        summary["description"] = "From publication " + line[line.find("\t"):].strip()

        # add the sample metadata
        cell_factors = None

        if sdrf_file:
            # extract all factors
            cell_factors = self._extract_factor_values_from_sdrf(sdrf_file=sdrf_file)

            # check if the file is valid
            first_cell_id = cell_clusterings[sample_ids[0]][0]

            if first_cell_id not in cell_factors:
                cell_factors = None

        # load the factors from experiment design file instead
        if not cell_factors:
            cell_factors = self._load_experiment_design_factors(dataset_id)

        # merge cell factors on the cluster level
        cluster_factors = dict()
        unique_factors = set()

        for cluster_id in cell_clusterings:
            this_cluster_factors = dict()

            for cell_id in cell_clusterings[cluster_id]:
                # simply ignore missing annotations - not all SDRF files are formatted the same way
                if cell_id not in cell_factors:
                    continue

                for factor_name in cell_factors[cell_id]:
                    cell_factor_value = cell_factors[cell_id][factor_name]

                    if factor_name not in this_cluster_factors:
                        this_cluster_factors[factor_name] = set()

                    this_cluster_factors[factor_name].add(cell_factor_value)
                    unique_factors.add(factor_name)

            cluster_factors[cluster_id] = this_cluster_factors

        # add the factors on the cluster level
        summary["sample_metadata"] = list()

        for factor_name in unique_factors:
            factor_values = list()

            for cluster_id in sample_ids:
                cluster_values = list(cluster_factors[cluster_id].get(factor_name, set()))

                # the property is the sorted list of all available properties
                cluster_values.sort()
                factor_values.append(",".join(cluster_values))

            summary["sample_metadata"].append({"name": factor_name, "values": factor_values})

        # return the object
        return ExternalData.from_dict(summary)

    def _load_experiment_design_factors(self, dataset_id: str) -> dict:
        """
        Downloads the experiment design file and returns all characteristics
        in a dict with the cell id as key and a second dict as value that
        contains all characteristics and values.
        :param dataset_id: The dataset's identifier
        :return: dict {cell_id => {param_name => param_value, ...}}
        """
        # Download the file
        logger.debug("Loading experimental design for {}...".format(dataset_id))

        exp_design_data = self._download_file(
            "https://www.ebi.ac.uk/gxa/sc/experiment/{}/download?fileType=experiment-design&accessKey="
            .format(dataset_id)).decode()

        # Clean the data
        clean_exp_design_data = expression_atlas_fetcher.ExpressionAtlasFetcher._filter_metadata(exp_design_data)

        # Convert to a numpy array
        exp_design = numpy.genfromtxt(io.StringIO(clean_exp_design_data), delimiter="\t", names = True, dtype = None)

        # process the data for every cell
        cell_factors = dict()

        for row in exp_design:
            # first column is the id
            cell_id = row[0].decode().replace("\"", "").strip()
            cell_factors[cell_id] = dict()

            # process every column
            for col_index in range(1, len(row)):
                col_name = exp_design.dtype.names[col_index]
                
                # save the value
                value = row[col_index].decode()
                value = value.replace("\"", "").strip()

                cell_factors[cell_id][col_name] = value

        return cell_factors

    def _extract_factor_values_from_sdrf(self, sdrf_file: str) -> dict:
        """
        Extracts all factor values from an SDRF file.
        :param sdrf_file: Path to the file
        :returns: A dict with the sample ids as keys and a second dict as value with the factors as key and the values
        """
        with open(sdrf_file, "r") as file_reader:
            dict_reader = csv.DictReader(file_reader, delimiter="\t")

            # get all factors
            factors = [factor for factor in dict_reader.fieldnames 
                        if "Factor Value[" in factor and "identifier" not in factor]
            factors_clean = [factor[13:-1] for factor in dict_reader.fieldnames 
                        if "Factor Value[" in factor and "identifier" not in factor]

            cell_factors = dict()

            for line in dict_reader:
                this_cell_dict = dict()

                for i in range(0, len(factors)):
                    this_cell_dict[factors_clean[i]] = line[factors[i]]

                cell_factors[line["Scan Name"]] = this_cell_dict

        return cell_factors
