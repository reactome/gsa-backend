"""
Contains all functions + abstract class for the Analyse implementations
"""

from reactome_analysis_worker.models import gene_set


class AnalysisException(Exception):
    """
    Exceptions related to the analysis.
    """
    pass


class ReactomeAnalyser:
    """
    Abstract class specifying how a ReactomeAnalyser should work.

    :ivar status_callback: The callback function to use to update the status. This variable should not be used
                           directly but set using `set_status_callback`. Analysers should then update their
                           status using the `_update_status` function.
    :ivar reactome_result_types: Array of `result_converter.ReactomeResultTypes` that are supported by the
                                 respective analyser. Empty list if the conversion is generally not supported.
    """
    _method_to_analyser = dict()

    def __init__(self):
        self.status_callback = None
        self.heartbeat_callback = None
        self.reactome_result_types = list()
        self.completion = 0

    def set_status_callback(self, callback):
        """
        Set the callback to use to allow status updates.
        :param callback: A callback following the schema (message: str, complete: float)
        """
        self.status_callback = callback

    def set_heartbeat_callback(self, callback):
        """
        A message that is called to simply indicate that the
        analysis is running
        """
        self.heartbeat_callback = callback


    def _heartbeat(self):
        """
        Heartbeat function to show that the analysis is still running.
        """
        if self.heartbeat_callback:
            self.heartbeat_callback()

    def _update_status(self, message: str, complete: float = None):
        """
        Uses the defined callback to update the status of the current process.
        :param message: The message to set
        :param complete: Percent completion (0 - 1). If not set, the previous completion value is used.
        """        
        # also triggers a heartbeat
        self._heartbeat()

        if self.status_callback:
            # use the previous completion if not set
            if not complete:
                complete = self.completion

            # update the previous completion
            self.completion = complete

            self.status_callback(message, complete)

    def analyse_request(self, request, gene_set_mappings: dict, identifier_mappings: dict, gene_set: gene_set.GeneSet):
        """
        Analyse the passed request. The result MUST contain the following columns in the pathways and
        fold_changes tables:

        ## Pathway table:

          * Pathway
          * Direction ["Up" / "Down"]
          * FDR

        ## Fold-change table:

          * Identifier
          * logFC
          * adj.P.Val

        :param request: The request sent by the user (as a `AnalysisInput` object)
        :param gene_set_mappings: A dict with the dataset names as keys and the gene set mapping result as value.
        :param identifier_mappings: The identifier mappings created by the REACTOME mapping service
        :param gene_set: The GeneSet used to perform the mapping.
        :return: The result as an `AnalysisResult` object.
        """
        raise NotImplementedError("This class must be implemented by any child class")

    @staticmethod
    def get_analyser_for_method(method: str):
        if method in ReactomeAnalyser._method_to_analyser:
            return ReactomeAnalyser._method_to_analyser[method]()

        return None

    @staticmethod
    def register_analyser(method: str, analyser):
        ReactomeAnalyser._method_to_analyser[method] = analyser

    @staticmethod
    def register_analyser_multiple(methods: list, analyser):
        for method in methods:
            ReactomeAnalyser.register_analyser(method, analyser)
