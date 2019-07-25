import logging
import os
import unittest

from reactome_analysis_report.reactome_analysis_report_generator import ReactomeAnalysisReportGenerator


class ReactomeAnalysisReportGeneratorTest(unittest.TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)

        test_dir = os.path.dirname(os.path.abspath(__file__))
        json_filename = os.path.join(test_dir, "result.json")

        os.environ["REDIS_HOST"] = "192.168.99.100"
        os.environ["REDIS_PORT"] = "31578"
        os.environ["REDIS_PASSWORD"] = "test"
        os.environ["RABBIT_HOST"] = "192.168.99.100"
        os.environ["RABBIT_PORT"] = "32428"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

        with open(json_filename, "r") as reader:
            self.test_json = reader.read()

    def test_process_all_messages(self):
        generator = ReactomeAnalysisReportGenerator()
        generator.start_listening()

    def test_xslx_creation(self):
        result_file = "/tmp/report_result.xlsx"

        ReactomeAnalysisReportGenerator.create_excel_file(analysis_result=self.test_json.encode(),
                                                          filename=result_file)

        self.assertTrue(os.path.isfile(result_file))
