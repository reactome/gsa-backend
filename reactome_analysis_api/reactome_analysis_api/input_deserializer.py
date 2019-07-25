import logging
from copy import deepcopy

from reactome_analysis_api import methods
from reactome_analysis_api.models.analysis_input import AnalysisInput

LOGGER = logging.getLogger(__name__)


def create_analysis_input_object(input_dict: dict) -> AnalysisInput:
    """
    Convert the dictionary retrieved from the analysis JSON object into
    a AnalysisInput object.

    This method circumvents the current bug in swagger that additional properties are
    not correctly supported. These are now stored in the "additional_properties" dict.

    Additionally, request parameters are also available in the parameter_dict object (instead
    of only the list).

    :param input_object: The AnalysisInput object
    :param input_dict: The dict representing the analysis input
    :return: The AnalysisInput object with additional properties.
    """
    if not isinstance(input_dict, dict):
        raise Exception("create_analysis_input_object must be called with a dict object.")

    input_object = AnalysisInput.from_dict(input_dict)

    # convert the request parameters to a dict
    LOGGER.debug("Existing user parameters:")
    params = dict()
    object_params = getattr(input_object, "parameters", list())
    if object_params:
        for parameter in object_params:
            LOGGER.debug(parameter.name + " = " + parameter.value)
            params[parameter.name] = parameter.value

    # set all non-existing parameters to their default value
    method_parameters = methods.get_parameters_for_method(method_name=input_object.method_name)

    if not method_parameters:
        raise Exception("Unknown analysis method '{}' selected".format(input_object.method_name))

    for parameter in method_parameters:
        if parameter.name not in params:
            params[parameter.name] = parameter.default

    input_object.parameter_dict = params

    # get the names of all dataset-level parameters
    dataset_parameter_names = set([parameter.name for parameter in method_parameters if parameter.scope == "dataset"])

    # get the current (default) values for all dataset-level parameters
    dataset_default_parameters = dict([(name, value) for (name, value) in input_object.parameter_dict.items()
                                       if name in dataset_parameter_names])

    # adapt every dataset
    for i in range(0, len(input_object.datasets)):
        existing_dataset_parameters = getattr(input_object.datasets[i], "parameters", list())
        # start with a copy of the default parameters
        dataset_parameter_dict = deepcopy(dataset_default_parameters)

        # if dataset level parameters are present, they overwrite the analysis wide ones
        if existing_dataset_parameters:
            for existing_parameter in existing_dataset_parameters:
                dataset_parameter_dict[existing_parameter.name] = existing_parameter.value

        input_object.datasets[i].parameter_dict = dataset_parameter_dict

    # fix the additional property issue in swagger
    for i in range(0, len(input_object.datasets)):
        # ignore datasets without design
        if "design" not in input_dict["datasets"][i]:
            continue

        for property in input_dict["datasets"][i]["design"]:
            if property not in ["samples", "comparison", "analysisGroup"]:
                if not hasattr(input_object.datasets[i].design, "additional_properties"):
                    input_object.datasets[i].design.additional_properties = dict()
                input_object.datasets[i].design.additional_properties[property] = \
                    input_dict["datasets"][i]["design"][property]

    return input_object