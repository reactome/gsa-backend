import logging
import time
import os
import multiprocessing
from contextlib import contextmanager


_LOGGER = logging.getLogger(__name__)


@contextmanager
def timeout_killer(alive_file: str, timeout: int) -> None:
    """Kills a running pod by deleting the keep alive
       file after a defined timeout. The process runs in a separate
       thread to work even if the main thread dies.

        :param alive_file: Path to the file to delete if the timeout is reached.
        :type alive_file: str
        :param timeout: The timeout in seconds after which the file is deleted.
        :type timeout: int
    """
    cancel_event = multiprocessing.Event()
    kill_process = TimeoutKillerProcess(alive_file=alive_file, timeout=timeout, cancel_event=cancel_event)
    kill_process.start()

    try:
        yield kill_process
    finally:
        cancel_event.set()


class TimeoutKillerProcess(multiprocessing.Process):
    """The process running in the background
    """
    def __init__(self, alive_file: str, timeout: int, cancel_event: multiprocessing.Event) -> None:
        """Initialize a new TimoutKilerProcess

        :param alive_file: Path to the file to delete
        :type alive_file: str
        :param timeout: The timeout in seconds after which to
                        delete the file
        :type timeout: int
        :param cancel_event: The event triggered if the file should not be deleted
        :type cancel_event: multiprocessing.Event
        """
        super().__init__()

        self._alive_file = alive_file
        self._timeout = timeout
        self._cancel_event = cancel_event
        

    def run(self) -> None:
        _LOGGER.debug("TimeoutKiller started. Waiting for %d seconds before deleting '%s'", 
                      self._timeout, self._alive_file)
        start_time = time.time()

        # check only every 10 seconds
        while start_time + self._timeout < time.time() and self.is_alive() and not self._cancel_event.is_set():
            time.sleep(10)

        # time is up, therefore delete the file
        if not self._cancel_event.is_set():
            _LOGGER.info("Timeout reached. Deleting keep alive file at '%s'", self._alive_file)
            os.unlink(self._alive_file)
        else:
            _LOGGER.debug("TimeoutKiller cancelled.")
