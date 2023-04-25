import unittest
import tempfile
import logging
import os
import time
from reactome_analysis_utils.timeout_killer import timeout_killer


class TimeoutKillerTest(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)

    def test_timeout_reached(self):
        # create the test file
        filename = tempfile.mktemp(prefix="timeout_killer_test")
        with open(filename, "w") as writer:
            writer.write("A")
        
        # create the killer
        start_time = time.time()
        with timeout_killer(alive_file=filename, timeout=2):
            while start_time + 1.5 > time.time():
                self.assertTrue(os.path.isfile(filename))
                time.sleep(0.1)

            # ensure that the timeout is reached
            time.sleep(0.7)

        # make sure the file was removed
        self.assertFalse(os.path.isfile(filename))

    def test_timeout_cancelled(self):
        # create the test file
        filename = tempfile.mktemp(prefix="timeout_killer_test")
        with open(filename, "w") as writer:
            writer.write("A")
        
        # create the killer
        with timeout_killer(alive_file=filename, timeout=2):
            # leave after 1 sec
            time.sleep(1)

        # ensure that the process ended
        time.sleep(1.1)

        # make sure the file was removed
        self.assertTrue(os.path.isfile(filename))
    