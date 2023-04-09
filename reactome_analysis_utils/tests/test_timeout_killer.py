import unittest
import tempfile
import logging
import os
import time
from reactome_analysis_utils.timeout_killer import TimeoutKiller, timeout_killer


class TimeoutKillerTest(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)

    def test_timeout_reached(self):
        # create the test file
        filename = tempfile.mktemp(prefix="timeout_killer_test")
        with open(filename, "w") as writer:
            writer.write("A")
        
        # create the killer
        with timeout_killer(alive_file=filename, timeout=1):
            time.sleep(2)

        # make sure the file was removed
        self.assertFalse(os.path.isfile(filename))

    def test_timeout_cancelled(self):
        # create the test file
        filename = tempfile.mktemp(prefix="timeout_killer_test")
        with open(filename, "w") as writer:
            writer.write("A")
        
        # create the killer
        with timeout_killer(alive_file=filename, timeout=2):
            pass

        time.sleep(2)

        # make sure the file was removed
        self.assertTrue(os.path.isfile(filename))
    