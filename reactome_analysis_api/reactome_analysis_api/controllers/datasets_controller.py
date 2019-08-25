import connexion
import six

from reactome_analysis_api.models.external_data import ExternalData  # noqa: E501
from reactome_analysis_api import util


def get_examples():  # noqa: E501
    """Lists the available example datasets

     # noqa: E501


    :rtype: ExternalData
    """
    # this method is currently only intended for local testing

    return [ExternalData(id="griss_melanoma_rnaseq", title="Griss Melanoma RNA-seq", type="rnaseq_counts", 
                         description="RNA-seq analysis of melanoma induced B cells.")]
