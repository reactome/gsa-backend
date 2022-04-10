from shutil import move
import unittest
import logging
from reactome_analysis_utils.reactome_logging import ReactomeSMTPHandler
from rediscluster.exceptions import MovedError


class TestReactomeSMTPHandler(unittest.TestCase):
    def test_ignore_log(self):
        test_logger = logging.getLogger(__name__)
        handler = ReactomeSMTPHandler(capacity=1)

        # even though the mail support is enabled, this should not cause a problem
        handler.has_mail_support = True
        test_logger.addHandler(handler)

        test_logger.exception("MovedError")
