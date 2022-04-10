"""
This class tests the validity of experimental designs
"""

import enum
import logging
from reactome_analysis_api.models.design import Design
import numpy


# default logger
logger = logging.getLogger(__name__)


class ExperimentalDesignExceptionType(enum.Enum):
    FORMAT_ERROR = 1
    SAME_GROUP = 2
    EMPTY_GROUP = 3
    FEW_SAMPLES = 4
    FEW_SAMPLES_COVAR = 5
    INVALID_DESIGN = 6


class ExperimentalDesignException(Exception):
    def __init__(self, msg: str, type: ExperimentalDesignExceptionType, dataset_name: str):
        """Initialise an ExperimentalDesignException.

        :param msg: The error message
        :type msg: str
        :param type: The type of error
        :type type: ExperimentalDesignExceptionType
        :param dataset_name: The dataset's name that is affected.
        :type dataset_name: str
        """
        super().__init__(msg)
        self.type = type
        self.dataset_name = dataset_name


class ExperimentalDesignTester:
    @staticmethod
    def test_experimental_designs(datasets: list) -> None:
        """
        Checks whether the structure of the experimental design matches the structure of
        the expression matrices in every dataset.

        Raises an ExperimentalDesignException in case the design is not valid.

        :param datasets: The datasets to test
        """
        # run the test for every dataset
        try:
            for dataset in datasets:
                # ignore requests that do not contain a design
                if not dataset.design:
                    continue

                # test for equal samples
                if not ExperimentalDesignTester._has_equal_samples(dataset.design, dataset.df):
                    raise ExperimentalDesignException(
                            msg="Different number of samples in the experimental design",
                            type=ExperimentalDesignExceptionType.FORMAT_ERROR,
                            dataset_name=dataset.name)

                # test for different comparison groups
                if not ExperimentalDesignTester._has_different_comparison_groups(dataset.design):
                    raise ExperimentalDesignException(
                            msg="Comparison group 1 and 2 must be different.",
                            type=ExperimentalDesignExceptionType.SAME_GROUP,
                            dataset_name=dataset.name)
                        
                # make sure the sample groups are present
                missing_group = ExperimentalDesignTester._find_missing_comparison_group(dataset.design)

                if missing_group:
                    raise ExperimentalDesignException(
                            msg="No samples annotated as comparison group '{group}'".format(group=missing_group),
                            type=ExperimentalDesignExceptionType.EMPTY_GROUP,
                            dataset_name=dataset.name)

                # make sure a suffient number of samples are annotated for every group
                too_few_samples_group = ExperimentalDesignTester._test_sufficient_samples_per_group(dataset.design)

                if too_few_samples_group:
                    raise ExperimentalDesignException(
                            msg="Insufficient samples in group '{group}'. Each group must at least contain 3 samples for accurate results.".format(group=too_few_samples_group),
                            type=ExperimentalDesignExceptionType.FEW_SAMPLES,
                            dataset_name=dataset.name)

                # make sure that additional groups have sufficient samples
                insufficient_samples_in_covariate = ExperimentalDesignTester._test_sufficient_samples_per_covariate(dataset.design)

                if insufficient_samples_in_covariate:
                    raise ExperimentalDesignException(
                            msg="Insufficient samples in covariate '{covariate}'. Each covariate must at least have two samples per group.".format(covariate=insufficient_samples_in_covariate),
                            type=ExperimentalDesignExceptionType.FEW_SAMPLES_COVAR,
                            dataset_name=dataset.name)
        except ExperimentalDesignException as err:
            # simply re-raise any already handled issue
            raise err
        except Exception as err:
            # this is unexpected, but also most likely caused by an invlid design
            logger.error("Failed to check experimental design: %s", str(err))
            logger.exception(err)

            raise ExperimentalDesignException(
                msg="Invalid experimental design", 
                type=ExperimentalDesignExceptionType.INVALID_DESIGN, 
                dataset_name="") from err

    @staticmethod
    def _has_equal_samples(design: Design, data: numpy.ndarray) -> bool:
        """Tests whether the samples mentioned in the design are present in the
        data matrix.

        :param design: The design to test
        :type design: Design
        :param data: The data as a numpy ndarray
        :type data: numpy.ndarray
        :return: Indicates whether the samples are equal
        :rtype: bool
        """
        design_samples = design.samples
        # first column contains the gene / protein id
        data_samples = data.dtype.names[1:]

        if not len(design_samples) == len(data_samples):
            return False

        return True

    @staticmethod
    def _has_different_comparison_groups(design: Design) -> bool:
        """Ensures that the comparison groups have different names

        :param design: The design to test
        :type design: Design
        :return: Indicates whether the test passed
        :rtype: bool
        """
        return design.comparison.group1 != design.comparison.group2

    @staticmethod
    def _find_missing_comparison_group(design: Design) -> str:
        """Checks whether the comparison groups are also present in the
        analysis_groups.

        :param design: The design to test
        :type design: Design
        :return: The missing comparison group's name or None if everything is fine.
        :rtype: str
        """
        analysis_groups = set(design.analysis_group)

        if design.comparison.group1 not in analysis_groups:
            return design.comparison.group1

        if design.comparison.group2 not in analysis_groups:
            return design.comparison.group2

        return None

    @staticmethod
    def _test_sufficient_samples_per_group(design: Design, min_samples:int = 3) -> str:
        """Checks whether at least min_samples are present in each comparison group.

        :param design: The design to test
        :type design: Design
        :param min_samples: Minimum number of samples required per group
        :type min_samples: int
        :return: The missing comparison group's name or None if everything is fine.
        :rtype: str
        """
        # get the number of samples per group
        n_group_1 = 0
        n_group_2 = 0

        for sample_group in design.analysis_group:
            if sample_group == design.comparison.group1:
                n_group_1 += 1
            if sample_group == design.comparison.group2:
                n_group_2 += 1

        if n_group_1 < min_samples:
            return design.comparison.group1

        if n_group_2 < min_samples:
            return design.comparison.group2

        return None

    @staticmethod
    def _test_sufficient_samples_per_covariate(design: Design) -> str:
        """Tests whether all defined covariate groups contain at least
        two samples.

        :param design: The design to test
        :type design: Design
        :return: Name of the group that does not contain enough samples (group name = instance name) or None if everything is fine.
        :rtype: str
        """
        # don't run the test if there are no additional properties
        if not hasattr(design, "additional_properties"):
            return None

        for property_name in design.additional_properties:
            # get the counts per sample
            instance_counts = dict()

            for sample_label in design.additional_properties[property_name]:
                if sample_label is None:
                    raise ExperimentalDesignException(
                        f"Missing sample label in {property_name}", 
                        ExperimentalDesignExceptionType.INVALID_DESIGN, 
                        "")

                # ignore empty lables
                sample_label = sample_label.strip()

                if len(sample_label) < 1:
                    continue

                # count the occurrence
                if sample_label not in instance_counts:
                    instance_counts[sample_label] = 1
                else:
                    instance_counts[sample_label] += 1

            # make sure all instances are > 1
            for instance_name in instance_counts:
                if instance_counts[instance_name] < 2:
                    return property_name + " = " + instance_name

        return None