import multiprocessing
import logging
import time
import queue

from reactome_analysis_worker import result_converter


LOGGER = logging.getLogger(__name__)


class ReactomeOraAnalyser:
    def __init__(self, identifiers: list, use_interactors: bool, reactome_server: str, include_disease: bool):
        """Fetch the reactome blueprint in a parallel process

        :param identifiers: Identifiers found
        :type identifiers: list
        :param use_interactors: Indicates whether interactors are used
        :type use_interactors: bool
        :param reactome_server: Server to use for the analysis
        :type reactome_server: str
        :param include_disease: Indicates whether disease pathways should be included
        :type include_disease: bool
        """
        self.on_complete = multiprocessing.Event()
        self.result_queue = multiprocessing.Queue()

        # create the process to fetch the result
        self.process = ReactomeAnalysisProcess(result_queue=self.result_queue,
                                               on_complete=self.on_complete, 
                                               identifiers=identifiers, 
                                               use_interactors=use_interactors, 
                                               reactome_server=reactome_server, 
                                               include_disease=include_disease)

    def start(self):
        """Start the process to fetch the reactome blueprint
        """
        self.process.start()

    def get_result(self, sleep_function=time.sleep) -> dict:
        """Retrieves the result. This function blocks until the result is retrieved.
           While waiting for the completion of the process, the sleep function is called.

        :param sleep_function: The function called to cause the thread to sleep. Must take one argument (seconds as float)
        :raises Exception: Any error results in a raised exception
        :return: The blueprint as a dict or None if nothing was returned
        :rtype: dict
        """
        while self.process.is_alive() and not self.on_complete.is_set():
            LOGGER.debug("Sleeping (Reactome blueprint)...")
            sleep_function(1)
                
        # for potential cleanup
        self.process.join(1)

        # retrieve the result from the queue
        try:
            result = self.result_queue.get(block=True, timeout=0.5)
        except queue.Empty:
            result = None

        # check if an exception was raised
        if result and isinstance(result, Exception):
            raise Exception(result)

        return result

class ReactomeAnalysisProcess(multiprocessing.Process):
    def __init__(self, result_queue: multiprocessing.Queue, on_complete: multiprocessing.Event,
                identifiers: list, use_interactors: bool, reactome_server: str, include_disease: bool):
        """Create a new process to fetch the reactome blueprint

        :param result_queue: Queue to use to post the result
        :type result_queue: multiprocessing.Queue
        :param on_complete: Set as soon as the process is done
        :type on_complete: multiprocessing.Event
        :param identifiers: Identifiers found
        :type identifiers: list
        :param use_interactors: Indicates whether interactors are used
        :type use_interactors: bool
        :param reactome_server: Server to use for the analysis
        :type reactome_server: str
        :param include_disease: Indicates whether disease pathways should be included
        :type include_disease: bool
        """
        super().__init__()
        
        self.result_queue = result_queue
        self.on_complete = on_complete
        self.identifiers = identifiers
        self.use_interactors = use_interactors
        self.reactome_server = reactome_server
        self.include_disease = include_disease

    def run(self) -> None:
        LOGGER.debug("Fetching blueprint for Reactome result conversion")
        try:
            reactome_blueprint = result_converter.perform_reactome_gsa(identifiers=self.identifiers,
                                                                        use_interactors=self.use_interactors,
                                                                        reactome_server=self.reactome_server,
                                                                        include_disease=self.include_disease)

            self.result_queue.put(reactome_blueprint)
        except Exception as e:
            self.result_queue.put(e)
        finally:
            self.on_complete.set()
