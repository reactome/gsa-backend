"""
R-based reactome analyser
"""

import logging
import gc
from collections.abc import Iterable

import rpy2.rinterface as ri
import rpy2.robjects as ro
from pkg_resources import resource_string, resource_listdir
from reactome_analysis_api.models.analysis_result_results import AnalysisResultResults
from rpy2.robjects.packages import importr

from reactome_analysis_worker.analysers.analyser import ReactomeAnalyser, AnalysisException
from reactome_analysis_worker.result_converter import ReactomeResultTypes

# initialize R
ri.initr()
# create a standard logger object
LOGGER = logging.getLogger(__name__)


def load_r_code_file(filename: str) -> ro.packages.SignatureTranslatedAnonymousPackage:
    """
    Load the code from the specified R code resource file and convert into
    an anonymously loaded package.

    :param filename: The filename to load. This file must recide in r_code
    :return: A SignatureTranslatedAnonymousPackage representing the package
    """
    # make sure the file is a resource file
    resource_files = resource_listdir('reactome_analysis_worker.resources.r_code', '')

    if filename not in resource_files:
        LOGGER.error("Failed to load R analyser code: {} does not exist.".format(filename))
        raise AnalysisException("{} is not a valid resource file".format(filename))

    # load the code from file
    r_code = resource_string('reactome_analysis_worker.resources.r_code', filename).decode("UTF-8")

    r_package = ro.packages.SignatureTranslatedAnonymousPackage(r_code, "reactome_" + filename[:-2])

    return r_package


class ReactomeRAnalyser(ReactomeAnalyser):
    methods = {"camera": load_r_code_file("camera_analyser.R"),
               "padog": load_r_code_file("padog_analyser.R")}
    preprocess = load_r_code_file("preprocessing_functions.R")

    def __init__(self):
        """
        Only used to initialize member variables
        """
        super().__init__()
        self.reactome_result_types = [ReactomeResultTypes.gsa]
        self.r_messages = list()

        # Note: This changed in rp2 version 3.x
        # see https://rpy2.github.io/doc/latest/html/callbacks.html?highlight=warn
        ri.set_writeconsole_warnerror(self._catch_message)
        ri.set_writeconsole_regular(self._catch_message)

    def _catch_message(self, message):
        """
        Function used to append messages received from R to the buffer
        :param message: The received message
        """
        # triggers a "heartbeat"
        self._heartbeat()

        message = message.strip()
        if len(message) == 0:
            return

        LOGGER.debug("R Message: " + message)
        self.r_messages.append(message)

    """
    Performs pathway analysis using different R methods
    """
    def analyse_request(self, request, gene_set_mappings, identifier_mappings, gene_set):
        # clean the environment
        ro.reval("rm(list=ls())")

        # get the pathway id to name mapping
        pathway_names = self.dict_of_list_to_r(gene_set.gene_set_names)

        # process every dataset separately
        analysis_results = list()
        previous_progress = 0.3

        for dataset in request.datasets:
            # make sure the dataset has a design
            if dataset.design is None:
                raise AnalysisException("Dataset '" + dataset.name + "' does not contain an experimental design.")

            # get the gene index
            gene_index = self.dict_of_list_to_r(gene_set_mappings[dataset.name].gene_set_indices, value_type=int)

            # prepare the dataset for the analysis - including pre-processing
            (expression_data, sample_data, design) = \
                self._prepare_dataset(dataset=dataset)

            self._update_status("Analysing dataset '{}' using {}".format(dataset.name, request.method_name),
                                complete=previous_progress + (0.3 / len(request.datasets)))

            result = self._perform_gsa(method=request.method_name,
                                       parameters=getattr(dataset, "parameter_dict", dict()),
                                       expression_data=expression_data, sample_data=sample_data, design=design,
                                       gene_index=gene_index, data_type=dataset.type,
                                       pathway_names=pathway_names,
                                       comparison_group_1=dataset.design.comparison.group1,
                                       comparison_group_2=dataset.design.comparison.group2)

            self._update_status("Analysing dataset '{}' using {}".format(dataset.name, request.method_name),
                                complete=previous_progress + (0.5 / len(request.datasets)))

            fold_changes = self._estimate_gene_fc(method=request.method_name,
                                                  parameters=getattr(dataset, "parameter_dict", dict()),
                                                  expression_data=expression_data, sample_data=sample_data,
                                                  design=design, data_type=dataset.type,
                                                  comparison_group_1=dataset.design.comparison.group1,
                                                  comparison_group_2=dataset.design.comparison.group2)

            self._update_status("Analysing dataset '{}' using {}".format(dataset.name, request.method_name),
                                complete=previous_progress + (0.7 / len(request.datasets)))

            # add average fold-changes to the analysis result
            result = ReactomeRAnalyser.preprocess.add_pathway_foldchanges(result, fold_changes, gene_index,
                                                                          expression_data)

            analysis_result = AnalysisResultResults(name=dataset.name,
                                                    pathways=str(ReactomeRAnalyser.preprocess.data_frame_as_string(result)[0]),
                                                    fold_changes=str(ReactomeRAnalyser.preprocess.data_frame_as_string(fold_changes)[0]))

            analysis_results.append(analysis_result)

            previous_progress += 0.7 / len(request.datasets)

        return analysis_results

    def _prepare_dataset(self, dataset):
        # create save sample names
        sample_names = self._create_save_names(dataset.design.samples)

        # convert the expression data to an R matrix
        expression_data = self._convert_dataset(dataset, sample_names)

        # create the sample annotation matrix
        sample_data = self._create_sample_data(dataset.design, sample_names)

        # create the model matrix
        LOGGER.debug("Creating design matrix...")
        design = ReactomeRAnalyser.preprocess.create_design(sample_data, group_1=dataset.design.comparison.group1)

        return expression_data, sample_data, design

    def _perform_gsa(self, method, parameters, expression_data, sample_data, design, gene_index, data_type,
                     pathway_names, comparison_group_1, comparison_group_2):
        """
        Perform the GSA on the specified dataset
        :param dataset:
        :param method:
        :param parameters:
        :return:
        """
        # parameters are currently stored in globalenv
        ri.globalenv["edger.norm.function"] = ri.StrSexpVector([parameters.get("discrete_norm_function", "TMM")])
        ri.globalenv["continuous.norm.function"] = ri.StrSexpVector([parameters.get("continuous_norm_function", "none")])
        ri.globalenv["sample.groups"] = ri.StrSexpVector([parameters.get("sample_groups", "")])

        # load the analyser R code
        LOGGER.debug("Processing data using R analysis code for {}".format(method.lower()))
        analysis_package = ReactomeRAnalyser.methods[method.lower()]

        # load the required limma package
        try:
            analysis_package.load_libraries()
        except ri.RRuntimeError as e:
            LOGGER.critical("Failed to load required package: " + str(e))
            raise AnalysisException("Failed to load required R package")

        # run the analysis
        try:
            gsa_result = analysis_package.process(expression_data, sample_data, design, gene_index,
                                                  ri.StrSexpVector([data_type]),
                                                  ri.StrSexpVector(["analysis_group" + comparison_group_1]),
                                                  ri.StrSexpVector(["analysis_group" + comparison_group_2]))
        except ri.RRuntimeError as e:
            LOGGER.debug("Failed to perform analysis: " + str(e))

            # return a more verbose message
            raise AnalysisException(". ".join(self.r_messages))

        # add the pathway names
        gsa_result = ReactomeRAnalyser.preprocess.add_pathway_names(gsa_result, pathway_names)

        return gsa_result

    def _estimate_gene_fc(self, method, parameters, expression_data, sample_data, design, data_type,
                          comparison_group_1, comparison_group_2):
        """
        Estimate the gene / protein-wise fold-change.
        :param dataset:
        :param method:
        :param parameters:
        :return:
        """
        # parameters are currently stored in globalenv
        ri.globalenv["edger.norm.function"] = ri.StrSexpVector([parameters.get("discrete_norm_function", "TMM")])
        ri.globalenv["continuous.norm.function"] = ri.StrSexpVector([parameters.get("continuous_norm_function", "none")])

        # load the analyser R code
        LOGGER.debug("Processing data using R analysis code for {}".format(method.lower()))
        analysis_package = ReactomeRAnalyser.methods[method.lower()]

        # load the required limma package
        try:
            analysis_package.load_libraries()
        except ri.RRuntimeError as e:
            LOGGER.critical("Failed to load required package: " + str(e))
            raise AnalysisException("Failed to load required R package")

        # run the analysis
        fc_result = analysis_package.get_gene_fc(expression_data, sample_data, design,
                                              ri.StrSexpVector([data_type]),
                                              ri.StrSexpVector(["analysis_group" + comparison_group_1]),
                                              ri.StrSexpVector(["analysis_group" + comparison_group_2]))

        return fc_result

    @staticmethod
    def dict_of_list_to_r(dict_to_convert, value_type=str):
        """
        Converts a dict with lists of strings as values to
        a named R list()
        :param dict_to_convert: The dict to convert
        :param value_type: The type of the values in the dict
        :return A named R list() object.
        """
        # convert to a dict of string vectors
        r_lists = dict()
        for item_name in dict_to_convert:
            value = dict_to_convert[item_name]

            # encapsulate single value items in lists
            if not isinstance(value, Iterable) or isinstance(value, str):
                value = [value]

            # convert non-list iterables to list
            if not isinstance(value, list):
                value = list(value)

            if value_type == str:
                r_lists[item_name] = ri.StrSexpVector(value)
            elif value_type == int:
                r_lists[item_name] = ri.IntSexpVector(value)
            elif value_type == float:
                r_lists[item_name] = ri.FloatSexpVector(value)
            else:
                raise AnalysisException("Failed to convert R vector of type " + str(value_type))

        rlist = ri.baseenv["list"]
        r_object = rlist(**r_lists)

        return r_object

    @staticmethod
    def _create_save_names(samples: list) -> ri.StrSexpVector:
        """
        Convert a list of strings to R save names.

        :param samples: A list of strings to convert
        :return: The save R names in an StrSexpVector
        """
        sample_names = ri.baseenv["make.names"](ri.StrSexpVector(samples))

        return sample_names

    @staticmethod
    def _convert_dataset(dataset, sample_names: ri.StrSexpVector):
        """
        Convert the dataset (a Dataset object) into an R data.frame with the
        passed sample_names used as colnames and the gene names as rownames.
        :param dataset: The dataset to convert.
        :param sample_names: A StrSexpVector containing the (converted) sample names to use
        :return: An R data.frame object
        """
        # first column contains the gene names
        gene_names = dataset.df[dataset.df.dtype.names[0]].tolist()

        # column names in the array contain the sample names
        colnames = dataset.df.dtype.names[1:]

        # turn all expression values into a single vector
        expression_values = list()
        for sample_name in colnames:
            expression_values += dataset.df[sample_name].tolist()

        # create the R vector
        r_vector = ri.FloatSexpVector(expression_values)

        rmatrix = ri.baseenv["matrix"]
        rlist = ri.baseenv["list"]

        if len(sample_names) != len(colnames):
            # this should never happen since it was checked before
            raise AnalysisException("Different number of samples in the experimental design and the expression matrix.")

        expression_data = rmatrix(r_vector, ncol=ri.IntSexpVector([len(colnames)]),
                                  dimnames=rlist(ri.StrSexpVector(gene_names), sample_names))

        # convert to data.frame
        rdataframe = ri.baseenv["data.frame"]

        return rdataframe(expression_data)

    @staticmethod
    def _create_sample_data(design, sample_names):
        """
        Convert the Design object into an R data.frame with one sample per row
        and every property stored in a corresponding column.
        :param design: The Design object to convert
        :param sample_names: The sample names to use
        :return: An R data.frame with the samples as rows and their properties as columns
        """
        # create a dict to hold all sample properties
        sample_properties = dict()

        # add the analysis group
        sample_properties["analysis_group"] = ri.StrSexpVector(design.analysis_group)

        # add additional properties
        for additional_property in getattr(design, "additional_properties", dict()):
            sample_properties[additional_property] = ri.StrSexpVector(
                design.additional_properties[additional_property])

        # convert the dict to a data.frame
        rdataframe = ri.baseenv["data.frame"]
        sample_data = rdataframe(**sample_properties)

        # use the passed sample_names as row names
        sample_data.do_slot_assign("row.names", sample_names)

        return sample_data


ReactomeAnalyser.register_analyser("camera", ReactomeRAnalyser)
ReactomeAnalyser.register_analyser("padog", ReactomeRAnalyser)
