import copy

from reactome_analysis_api.models.method import Method
from reactome_analysis_api.models.method_parameters import MethodParameters

global_parameters = [
        MethodParameters(name="use_interactors",
                         display_name="Use interactors",
                         type="bool",
                         default="False",
                         scope="analysis",
                         description="Indicates whether interactors from IntAct should be used to extent REACTOME's "
                                     "pathways in the analysis."),
        MethodParameters(name="include_disease_pathways",
                         display_name="Include disease pathways",
                         type="bool",
                         default="True",
                         scope="analysis",
                         description="Disease pathways in Reactome may lead to a skewed analysis result. Note: Excluding disease pathways "
                                     "currently prevents the visualization of results in the Reactome pathway browser."),
        MethodParameters(name="max_missing_values",
                         display_name="Max. missing values",
                         type="float",
                         default="0.5",
                         scope="dataset",
                         description="The maximum (relative) number of missing values within one comparison "
                                     "group before a gene / protein is removed from analysis. If no comparison groups "
                                     "are defined (for example for 'ssGSEA', the number of missing values accross all "
                                     "samples is used. Must be between 0-1."),
        MethodParameters(name="create_reactome_visualization",
                         display_name="Create REACTOME visualizations",
                         type="bool",
                         scope="common",
                         default="True",
                         description="If set to 'False', no REACTOME visualization is created for the performed "
                                     "analysis."),
        MethodParameters(name="create_reports",
                         display_name="Create reports",
                         type="bool",
                         scope="common",
                         default="False",
                         description="If set to 'True', additional Microsoft Excel and PDF-based reports of the "
                                     "analysis result will be created."),
        MethodParameters(name="email",
                         display_name="E-Mail",
                         type="string",
                         scope="common",
                         default="",
                         description="If set to a valid e-mail address, links to the analysis result (and report) will "
                                     "be sent once the analysis is complete."),
        MethodParameters(name="reactome_server",
                         display_name="Reactome server",
                         type="string",
                         scope="common",
                         default="production",
                         description="This parameter allows the usage of other reactome servers. Available options are "
                                     "'production', 'dev', 'release'")
    ]

available_methods = [
    Method(name="PADOG", description="Weighted gene set analysis method that down-weighs genes that are present in many"
                                     " pathways. Supports multiple Omics data sources including Ribo-Seq data", parameters=[
        MethodParameters(name="sample_groups",
                         display_name="Sample Groups",
                         type="string",
                         scope="dataset",
                         description="Specifies the sample property name that holds the sample group information. "
                                     "This parameter should be used for matched-pair analyses (f.e., the same patients "
                                     "before and after therapy). If used, every sample must occur exactly twice, once "
                                     "in each of the analysis groups.",
                         default=""),
        MethodParameters(name="discrete_norm_function",
                         display_name="Discrete normalisation function",
                         values=["TMM", "RLE", "upperquartile", "none"],
                         type="string",
                         scope="dataset",
                         default="TMM",
                         description="The normalisation function to use for raw RNA-seq read counts and "
                                     "raw Proteomics spectral counts. By default, the TMM normalisation "
                                     "is used."),
        MethodParameters(name="continuous_norm_function",
                        display_name="Continuous normalisation function",
                        values=["none", "scale", "quantile", "cyclicloess"],
                        type="string",
                        scope="dataset",
                        default="none",
                        description="The normalisation function to use for proteomics intensity data. Note "
                                    "that it is generally advised that normalisation is performed on the "
                                    "PSM or peptide level and not on the protein level.")
    ]),

    Method(name="Camera", description="A gene set analysis algorithm similar to the classical GSEA algorithm "
                                      "as implemented in the limma package.",
           parameters=[
               MethodParameters(name="discrete_norm_function",
                                display_name="Discrete normalisation function",
                                values=["TMM", "RLE", "upperquartile", "none"],
                                type="string",
                                scope="dataset",
                                default="TMM",
                                description="The normalisation function to use for raw RNA-seq read counts and "
                                            "raw Proteomics spectral counts. By default, the TMM normalisation "
                                            "is used."),
               MethodParameters(name="continuous_norm_function",
                                display_name="Continuous normalisation function",
                                values=["none", "scale", "quantile", "cyclicloess"],
                                type="string",
                                scope="dataset",
                                default="none",
                                description="The normalisation function to use for proteomics intensity data. Note "
                                            "that it is generally advised that normalisation is performed on the "
                                            "PSM or peptide level and not on the protein level.")
           ]),
    Method(name="ssGSEA", description="The ssGSEA approach to derive pathway expression values for every sample. " \
                                      "Note: The Reactome visualization is only available for up to 15 samples.",
           parameters=[
               MethodParameters(name="pathways",
                                display_name="Pathways",
                                type="string",
                                scope="analysis",
                                description="A comma delimited list of pathways to include in the analysis. All other "
                                            "pathways will be ignored.",
                                default=""),
               MethodParameters(name="min_size",
                                display_name="Minimum pathway size",
                                type="int",
                                scope="analysis",
                                description="The minimum pathway size (determined as the number of submitted genes "
                                            "mapped to that pathway) to include a pathway in the analysis.",
                                default="1"),
               MethodParameters(name="max_size",
                                display_name="Maximum pathway size",
                                type="int",
                                scope="analysis",
                                description="The maximum pathway size (determined as the number of submitted genes "
                                            "mapped to that pathway) to include a pathway in the analysis.",
                                default="1000")
           ])
]


def get_available_methods():
    """
    Returns all available methods with the global parameters added
    to each method's parameters.
    :return: A list of Methods
    """
    methods = copy.deepcopy(available_methods)

    # add the global parameters to all methods (first)
    for method in methods:
        method.parameters = global_parameters + method.parameters

    return methods


def get_parameters_for_method(method_name: str) -> list():
    """
    Returns a list of parameters defined for the specified method. Returns None if
    no method with that name exists.
    :param method_name: The method's name (case-insensitive)
    :return: A list of MethodParameters
    """
    method_name = method_name.lower().strip()

    all_methods = get_available_methods()

    for method in all_methods:
        if method.name.lower() == method_name:
            return method.parameters

    return None
