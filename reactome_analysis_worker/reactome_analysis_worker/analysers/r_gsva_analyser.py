"""
R-based reactome analyser
"""

import logging

import rpy2.rinterface as ri
import rpy2.robjects as ro
from reactome_analysis_api.models.analysis_input import AnalysisInput
from reactome_analysis_api.models.analysis_result_results import AnalysisResultResults

from reactome_analysis_worker.analysers.analyser import ReactomeAnalyser, AnalysisException
from reactome_analysis_worker.analysers.r_analyser import load_r_code_file, ReactomeRAnalyser
from reactome_analysis_worker.result_converter import ReactomeResultTypes

# initialize R
ri.initr()
# create a standard logger object
LOGGER = logging.getLogger(__name__)


class ReactomeGSVARAnalyser(ReactomeRAnalyser):
    # maximum number of samples allowed before the visualization is disabled
    MAX_SAMPLES = 15

    methods = {"ssgsea": load_r_code_file("gsva_ssgsea_analyser.R")}

    def __init__(self):
        super().__init__()
        # only set the correct supported result types
        self.reactome_result_types = [ReactomeResultTypes.gsva]

    """
    Performs GSVA analysis using different R methods
    """
    def analyse_request(self, request: AnalysisInput, gene_set_mappings, identifier_mappings, gene_set):
        # clean the environment
        ro.reval("rm(list=ls())")

        # get the pathway id to name mapping
        pathway_names = self.dict_of_list_to_r(gene_set.gene_set_names)

        # load the analyser R code
        LOGGER.debug("Processing data using R analysis code for {}".format(request.method_name.lower()))
        analysis_package = self.methods[request.method_name.lower()]

        LOGGER.debug("Retrieved analysis package")

        # load the libraries
        try:
            analysis_package.load_libraries()
        except Exception as e:
            LOGGER.critical("Failed to load required package: " + str(e))
            raise AnalysisException("Failed to load required R package")

        LOGGER.debug("R libraries loaded")

        # get the analysis-level parameters
        analysis_parameters = getattr(request, "parameter_dict", dict())

        # indicates whether the visualization should be disabled
        disable_visualization = False

        # if pathways are filtered using a list of pathways disable visualization
        if len(analysis_parameters.get("pathways", "")) > 0:
            disable_visualization = True

        # process every dataset separately
        analysis_results = list()
        previous_progress = 0.3

        for dataset in request.datasets:
            # create save sample names
            org_names = dataset.design.samples if dataset.design else dataset.df.dtype.names[1:]
            sample_names = self._create_save_names(org_names)

            # if there are more then MAX_SAMPLES, disable the visualization
            if len(sample_names) > ReactomeGSVARAnalyser.MAX_SAMPLES:
                disable_visualization = True

            LOGGER.debug("Converting expression data")

            # convert the expression data to an R matrix
            expression_data = self._convert_dataset(dataset, sample_names)

            # convert the fold_changes to the text representation
            # pylint: disable=no-member
            expression_data_id = ReactomeRAnalyser.preprocess.change_first_column(expression_data, 
                                                                               rowname_column=ri.StrSexpVector(["Identifier"]))
            r_fold_change_text = ReactomeRAnalyser.data_frame_to_string(expression_data_id)

            LOGGER.debug("Converting gene_index")

            # get the gene index
            gene_index = self.dict_of_list_to_r(gene_set_mappings[dataset.name].gene_set_indices, value_type=int)

            self._update_status("Analysing dataset '{}' using {}".format(dataset.name, request.method_name),
                                complete=previous_progress + (0.3 / len(request.datasets)))

            # perform the analysis
            LOGGER.debug("Starting GSVA analysis for {}".format(dataset.name))

            # use float before int to support scientific notation for max_size. This happens for large
            # numbers in the R package
            max_size = int(float(analysis_parameters.get("max_size", 1_000_000)))

            result = analysis_package.process(expression_data, gene_index, ri.StrSexpVector([dataset.type]),
                                              # These parameters are currently not visible to the user as
                                              # it might cause inconsistencies in the reactome result conversion
                                              ri.IntSexpVector([int(analysis_parameters.get("min_size", 0))]),
                                              ri.IntSexpVector([max_size]),
                                              ri.StrSexpVector([analysis_parameters.get("pathways", "")]))

            # add the pathway's name
            # pylint: disable=no-member
            result = ReactomeRAnalyser.preprocess.add_pathway_names(result, pathway_names)

            LOGGER.debug("GSVA analysis completed for {}".format(dataset.name))

            self._update_status("Analysing dataset '{}' using {}".format(dataset.name, request.method_name),
                                complete=previous_progress + (0.6 / len(request.datasets)))

            # convert the data.frame to a string
            r_text_result = ReactomeRAnalyser.data_frame_to_string(result)

            # add the result
            analysis_results.append(AnalysisResultResults(name=dataset.name,
                                                          pathways=r_text_result,
                                                          fold_changes=r_fold_change_text))

            previous_progress += 0.7 / len(request.datasets)

        LOGGER.debug("Returning combined analysis result")

        # disable the visualization if set
        if disable_visualization:
            if not hasattr(request, "parameter_dict"):
                request.parameter_dict = dict()

            request.parameter_dict["create_reactome_visualization"] = "False"

        return analysis_results


ReactomeAnalyser.register_analyser("ssgsea", ReactomeGSVARAnalyser)
