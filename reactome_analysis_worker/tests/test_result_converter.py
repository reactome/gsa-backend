import json
import unittest

from reactome_analysis_api.models.analysis_result import AnalysisResult, AnalysisResultResults

from reactome_analysis_worker import result_converter
from reactome_analysis_worker import util


class ResultConverterTest(unittest.TestCase):
    def setUp(self) -> None:
        self.gsva_result = AnalysisResult(results=[
            AnalysisResultResults(name="proteomics",
                                  pathways="""Pathway\tName\tSample_1\tSample_2\tSample_3
R-HSA-3232118\tName\t1\t2\t3
R-HSA-196807\tName\t0.1\t0.001\t0.02
R-HSA-2219530\tName\t5\t0.5\t0.01
""", fold_changes="""Identifier\tSample_1\tSample_2\tSample_3
CD19\t0.1\t0.2\t2
CD20\t1\t2\t2
MITF\t2\t3\t4
CD38\t1\t2\t3
"""),
            AnalysisResultResults(name="rnaseq",
                                  pathways="""Pathway\tName\tSample_5\tSample_6\tSample_7
R-HSA-3232118\tName\t2\t0.01\t0.00231
R-HSA-196807\tName\t3\t0.1\t0.1231
R-HSA-2219530\tName\t4\t0.5\t0.0001231
R-HSA-983695\tName\t5\t0.01\t0.1231
R-HSA-6811558\tName\t6\t0.02\t0.00001231
""", fold_changes="""Identifier\tSample_1\tSample_2\tSample_3
CD19\t0.1\t0.2\t2
CD20\t1\t2\t2
MITF\t2\t3\t4
CD38\t1\t2\t3
""")
        ])

    def test_reactome_gsa(self):
        identifier = ["MITF", "MS4A1", "CD19"]

        result = result_converter.perform_reactome_gsa(identifiers=identifier)

        self.assertIsNotNone(result)
        self.assertIsNotNone(result["summary"])
        self.assertTrue(result["summary"]["projection"])
        self.assertFalse(result["summary"]["interactors"])

        self.assertEqual(23, len(result["pathways"]))

    def test_get_pathway_changes(self):
        pathway_fc_string_1 = "Pathway\tDirection\tFDR\nP1\tUp\t0.02\nP2\tDown\t0.01\nP3\tUp\t0.07\nP4\tDown\t0.09\n"
        pathway_fc_string_2 = "Pathway\tDirection\tFDR\nP1\tDown\t0.02\nP2\tUp\t0.01\nP5\tDown\t0.07\nP4\tUp\t0.09\n"

        pathway_fcs = [util.string_to_array(pathway_fc_string_1), util.string_to_array(pathway_fc_string_2)]
        all_pathways = {"P1", "P2", "P3", "P4", "P5"}

        pathway_changes = result_converter._get_pathway_changes(pathway_fcs, all_pathways, 0.05)

        self.assertIsNotNone(pathway_changes)
        self.assertEqual(5, len(pathway_changes))

        for p in all_pathways:
            self.assertTrue(p in pathway_changes)
            self.assertEqual(2, len(pathway_changes[p]))

    def test_get_pathway_changes_p(self):
        pathway_fc_string_1 = "Pathway\tDirection\tFDR\nP1\tUp\t0.02\nP2\tDown\t0.01\nP3\tUp\t0.07\nP4\tDown\t0.09\n"
        pathway_fc_string_2 = "Pathway\tDirection\tFDR\nP1\tDown\t0.02\nP2\tUp\t0.01\nP5\tDown\t0.07\nP4\tUp\t0.09\n"

        pathway_fcs = [util.string_to_array(pathway_fc_string_1), util.string_to_array(pathway_fc_string_2)]
        all_pathways = {"P1", "P2", "P3", "P4", "P5"}

        pathway_changes = result_converter._get_pathway_changes(pathway_fcs, all_pathways, 0.05, return_p=True)

        self.assertIsNotNone(pathway_changes)
        self.assertEqual(5, len(pathway_changes))

        for p in all_pathways:
            self.assertTrue(p in pathway_changes)
            self.assertEqual(2, len(pathway_changes[p]))

        self.assertEqual(0.02, pathway_changes["P1"][0])
        self.assertEqual(0.02, pathway_changes["P1"][1])
        self.assertEqual(0.01, pathway_changes["P2"][0])
        self.assertEqual(0.01, pathway_changes["P2"][1])
        self.assertEqual(0.07, pathway_changes["P3"][0])
        self.assertEqual(1, pathway_changes["P3"][1])

    def test_get_identifier_changes(self):
        identifier_fc_string_1 = "Identifier\tlogFC\tadj.P.Val\nP1\t0.01\t1\nP2\t0.02\t1\nP3\t0.03\t1\n"
        identifier_fc_string_2 = "Identifier\tlogFC\tadj.P.Val\nP1\t0.01\t1\nP2\t0.02\t1\nP4\t0.03\t1\n"

        all_identifiers = {"P1", "P2", "P3", "P4"}
        identifier_fcs = [util.string_to_array(identifier_fc_string_1), util.string_to_array(identifier_fc_string_2)]

        identifier_changes = result_converter._get_identifier_changes(identifier_fcs, all_identifiers)

        self.assertIsNotNone(identifier_changes)
        self.assertEqual(4, len(identifier_changes))

        for i in all_identifiers:
            self.assertTrue(i in identifier_changes)
            self.assertEqual(2, len(identifier_changes[i]))

    def test_get_identifier_p_values(self):
        identifier_fc_string_1 = "Identifier\tlogFC\tadj.P.Val\nP1\t0.01\t1\nP2\t0.02\t1\nP3\t0.03\t0.01\n"
        identifier_fc_string_2 = "Identifier\tlogFC\tadj.P.Val\nP1\t0.01\t1\nP2\t0.02\t1\nP4\t0.03\t1\n"

        all_identifiers = {"P1", "P2", "P3", "P4"}
        identifier_fcs = [util.string_to_array(identifier_fc_string_1), util.string_to_array(identifier_fc_string_2)]

        identifier_changes = result_converter._get_identifier_changes(identifier_fcs, all_identifiers, return_p=True)

        self.assertIsNotNone(identifier_changes)
        self.assertEqual(4, len(identifier_changes))

        for i in all_identifiers:
            self.assertTrue(i in identifier_changes)
            self.assertEqual(2, len(identifier_changes[i]))

        self.assertEqual(1, identifier_changes["P1"][0])
        self.assertEqual(1, identifier_changes["P1"][1])
        self.assertEqual(1, identifier_changes["P2"][0])
        self.assertEqual(1, identifier_changes["P2"][1])
        self.assertEqual(0.01, identifier_changes["P3"][0])
        self.assertEqual(1, identifier_changes["P3"][1])

    def test_get_pathway_p_values(self):
        pathway_fc_string_1 = "Pathway\tDirection\tPValue\nP1\tUp\t0.02\nP2\tDown\t0.01\nP3\tUp\t0.07\nP4\tDown\t0.09\n"
        pathway_fc_string_2 = "Pathway\tDirection\tPValue\nP1\tUp\t0.02\nP2\tUp\t0.01\nP5\tDown\t0.07\nP4\tUp\t0.09\n"

        pathway_fcs = [util.string_to_array(pathway_fc_string_1), util.string_to_array(pathway_fc_string_2)]

        pathway_p = result_converter._get_pathway_p_values(pathway_fcs)

        self.assertIsNotNone(pathway_p)
        self.assertEqual(5, len(pathway_p))

        for pathway in {"P1", "P2", "P3", "P4", "P5"}:
            self.assertTrue(pathway in pathway_p)

        self.assertGreater(0.02, pathway_p["P1"]["p"])
        self.assertLess(0.01, pathway_p["P2"]["p"])
        self.assertEqual(0.07, pathway_p["P3"]["p"])

    def test_convert_gsa_result(self):
        identifier_fc_string_1 = "Identifier\tlogFC\tadj.P.Val\nMITF\t0.01\t1\nCD20\t0.02\t1\nCD19\t0.03\t1\n"
        identifier_fc_string_2 = "Identifier\tlogFC\tadj.P.Val\nMITF\t0.01\t1\nCD38\t-0.02\t0.001\nCD20\t0.03\t1\n"

        pathway_fc_string_1 = """Pathway\tDirection\tFDR\tPValue
R-HSA-3232118\tUp\t0.01\t0.02
R-HSA-196807\tDown\t0.001\t0.02
R-HSA-2219530\tUp\t0.5\t0.01
R-HSA-983695\tUp\t0.01\t0.001
R-HSA-6811558\tDown\t0.02\t0.002
R-HSA-199418\tUp\t0.09\t0.004
R-HSA-3108232\tUp\t0.03\t0.1
R-HSA-983705\tDown\t0.01\t0.003
R-HSA-2990846\tDown\t0.01\t0.002
R-HSA-196849\tDown\t0.01\t0.004
R-HSA-198933\tDown\t0.01\t0.123
R-HSA-1257604\tDown\t0.01\t0.00123
R-HSA-9006925\tDown\t0.01\t0.002321
R-HSA-196854\tDown\t0.01\t0.231
R-HSA-5663202\tDown\t0.01\t0.12312
R-HSA-1280218\tDown\t0.01\t0.0002131
R-HSA-168249\tDown\t0.01\t0.0123
R-HSA-1643685\tDown\t0.01\t0.001231
R-HSA-597592\tDown\t0.01\t0.00231
R-HSA-392499\tDown\t0.01\t0.0123
R-HSA-168256\tDown\t0.01\t0.1231
R-HSA-162582\tDown\t0.01\t0.1231
R-HSA-1430728\tDown\t0.01\t0.123
"""
        pathway_fc_string_2 = """Pathway\tDirection\tFDR\tPValue
R-HSA-3232118\tUp\t0.01\t0.00231
R-HSA-196807\tDown\t0.1\t0.1231
R-HSA-2219530\tUp\t0.5\t0.0001231
R-HSA-983695\tUp\t0.01\t0.1231
R-HSA-6811558\tDown\t0.02\t0.00001231
R-HSA-199418\tUp\t0.09\t0.1231
R-HSA-2219528\tDown\t0.02\t0.0001231
R-HSA-977606\tDown\t0.01\t0.321
R-HSA-166658\tUp\t0.01\t0.00012311
R-HSA-3108232\tUp\t0.03\t0.01231
R-HSA-983705\tDown\t0.01\t0.01231
R-HSA-2990846\tDown\t0.01\t0.1231
R-HSA-196849\tDown\t0.01\t0.0001231
R-HSA-198933\tDown\t0.01\t0.1231
R-HSA-1257604\tUp\t0.01\t0.1231
R-HSA-9006925\tUp\t0.01\t0.0001231
R-HSA-196854\tUp\t0.01\t0.0001231
R-HSA-5663202\tUp\t0.01\t0.01231
R-HSA-1280218\tDown\t0.01\t0.0002131
R-HSA-168249\tDown\t0.01\t0.23232
R-HSA-1643685\tUp\t0.01\t0.00001231
R-HSA-168256\tUp\t0.01\t0.01231
R-HSA-162582\tDown\t0.01\t0.001231
R-HSA-1430728\tDown\t0.01\t0.012312
"""
        result = result_converter.perform_reactome_gsa(["MITF", "CD20", "CD19", "CD38"])

        self.assertIsNotNone(result)

        # create the analysis result object
        analysis_result = AnalysisResult(results=[
            AnalysisResultResults(name="proteomics", pathways=pathway_fc_string_1, fold_changes=identifier_fc_string_1),
            AnalysisResultResults(name="rnaseq", pathways=pathway_fc_string_2, fold_changes=identifier_fc_string_2)
        ])

        converted_result = result_converter._convert_gsa_result(result=analysis_result, reactome_blueprint=result)

        # make sure the original object is not changed
        self.assertEqual("OVERREPRESENTATION", result["summary"]["type"])
        self.assertEqual("GSA_REGULATION", converted_result["summary"]["type"])

        with open("/tmp/test.json", "w") as writer:
            writer.write(json.dumps(converted_result))

        self.assertIsNotNone(converted_result)

        self.assertEqual("R-HSA-3232118", converted_result["pathways"][0]["stId"])
        self.assertEqual("MITF", converted_result["pathways"][0]["data"]["entities"][0]["id"])
        self.assertEqual(1, converted_result["pathways"][0]["data"]["entities"][0]["exp"][0])
        self.assertEqual(1, converted_result["pathways"][0]["data"]["entities"][0]["exp"][1])
        self.assertEqual(2, converted_result["pathways"][0]["data"]["statistics"][0]["exp"][0])
        self.assertEqual(2, converted_result["pathways"][0]["data"]["statistics"][0]["exp"][1])

        self.assertEqual("R-HSA-196807", converted_result["pathways"][1]["stId"])
        self.assertEqual("CD38", converted_result["pathways"][1]["data"]["entities"][0]["id"])
        self.assertEqual(0, converted_result["pathways"][1]["data"]["entities"][0]["exp"][0])
        self.assertEqual(-2, converted_result["pathways"][1]["data"]["entities"][0]["exp"][1])
        self.assertEqual(-2, converted_result["pathways"][1]["data"]["statistics"][0]["exp"][0])
        self.assertEqual(-1, converted_result["pathways"][1]["data"]["statistics"][0]["exp"][1])

    def test_get_gsva_pathway_expression(self):
        pathway_expression = result_converter._get_gsva_pathway_expression(self.gsva_result)

        self.assertIsNotNone(pathway_expression)
        for p in pathway_expression:
            self.assertEqual(6, len(pathway_expression[p]))

        self.assertEqual(1, pathway_expression["R-HSA-3232118"][0])
        self.assertEqual(2, pathway_expression["R-HSA-3232118"][1])
        self.assertEqual(3, pathway_expression["R-HSA-3232118"][2])
        self.assertEqual(2, pathway_expression["R-HSA-3232118"][3])
        self.assertEqual(0.01, pathway_expression["R-HSA-3232118"][4])
        self.assertEqual(0.00231, pathway_expression["R-HSA-3232118"][5])

        self.assertEqual(0, pathway_expression["R-HSA-6811558"][0])
        self.assertEqual(0, pathway_expression["R-HSA-6811558"][1])
        self.assertEqual(0, pathway_expression["R-HSA-6811558"][2])

    def test_get_identifier_zscores(self):
        identifier_expression = result_converter._get_identifier_zscores(self.gsva_result)

        self.assertIsNotNone(identifier_expression)
        self.assertEqual(4, len(identifier_expression))
        for i in identifier_expression:
            self.assertEqual(6, len(identifier_expression[i]))

    def test_submit_result_to_reactome(self):
        """
        This tests the complete upload of a result object to REACTOME
        :return:
        """
        identifier_fc_string_1 = "Identifier\tlogFC\tadj.P.Val\nMITF\t0.01\t1\nCD20\t0.02\t1\nCD19\t0.03\t1\n"
        identifier_fc_string_2 = "Identifier\tlogFC\tadj.P.Val\nMITF\t0.01\t1\nCD38\t-0.02\t0.001\nCD20\t0.03\t1\n"

        pathway_fc_string_1 = """Pathway\tDirection\tFDR\tPValue
        R-HSA-3232118\tUp\t0.01\t0.02
        R-HSA-196807\tDown\t0.001\t0.02
        R-HSA-2219530\tUp\t0.5\t0.01
        R-HSA-983695\tUp\t0.01\t0.001
        R-HSA-6811558\tDown\t0.02\t0.002
        R-HSA-199418\tUp\t0.09\t0.004
        R-HSA-3108232\tUp\t0.03\t0.1
        R-HSA-983705\tDown\t0.01\t0.003
        R-HSA-2990846\tDown\t0.01\t0.002
        R-HSA-196849\tDown\t0.01\t0.004
        R-HSA-198933\tDown\t0.01\t0.123
        R-HSA-1257604\tDown\t0.01\t0.00123
        R-HSA-9006925\tDown\t0.01\t0.002321
        R-HSA-196854\tDown\t0.01\t0.231
        R-HSA-5663202\tDown\t0.01\t0.12312
        R-HSA-1280218\tDown\t0.01\t0.0002131
        R-HSA-168249\tDown\t0.01\t0.0123
        R-HSA-1643685\tDown\t0.01\t0.001231
        R-HSA-597592\tDown\t0.01\t0.00231
        R-HSA-392499\tDown\t0.01\t0.0123
        R-HSA-168256\tDown\t0.01\t0.1231
        R-HSA-162582\tDown\t0.01\t0.1231
        R-HSA-1430728\tDown\t0.01\t0.123
        """
        pathway_fc_string_2 = """Pathway\tDirection\tFDR\tPValue
        R-HSA-3232118\tUp\t0.01\t0.00231
        R-HSA-196807\tDown\t0.1\t0.1231
        R-HSA-2219530\tUp\t0.5\t0.0001231
        R-HSA-983695\tUp\t0.01\t0.1231
        R-HSA-6811558\tDown\t0.02\t0.00001231
        R-HSA-199418\tUp\t0.09\t0.1231
        R-HSA-2219528\tDown\t0.02\t0.0001231
        R-HSA-977606\tDown\t0.01\t0.321
        R-HSA-166658\tUp\t0.01\t0.00012311
        R-HSA-3108232\tUp\t0.03\t0.01231
        R-HSA-983705\tDown\t0.01\t0.01231
        R-HSA-2990846\tDown\t0.01\t0.1231
        R-HSA-196849\tDown\t0.01\t0.0001231
        R-HSA-198933\tDown\t0.01\t0.1231
        R-HSA-1257604\tUp\t0.01\t0.1231
        R-HSA-9006925\tUp\t0.01\t0.0001231
        R-HSA-196854\tUp\t0.01\t0.0001231
        R-HSA-5663202\tUp\t0.01\t0.01231
        R-HSA-1280218\tDown\t0.01\t0.0002131
        R-HSA-168249\tDown\t0.01\t0.23232
        R-HSA-1643685\tUp\t0.01\t0.00001231
        R-HSA-168256\tUp\t0.01\t0.01231
        R-HSA-162582\tDown\t0.01\t0.001231
        R-HSA-1430728\tDown\t0.01\t0.012312
        """
        result = result_converter.perform_reactome_gsa(["MITF", "CD20", "CD19", "CD38"])

        self.assertIsNotNone(result)

        # create the analysis result object
        analysis_result = AnalysisResult(results=[
            AnalysisResultResults(name="proteomics", pathways=pathway_fc_string_1, fold_changes=identifier_fc_string_1),
            AnalysisResultResults(name="rnaseq", pathways=pathway_fc_string_2, fold_changes=identifier_fc_string_2)
        ])

        token = result_converter.submit_result_to_reactome(result=analysis_result,
                                                           result_type=result_converter.ReactomeResultTypes.gsa,
                                                           reactome_blueprint=result)

        self.assertIsNotNone(token)

        print("Reactome result '{}' ({}) available here: {}".format(token.name, token.description, token.url))
