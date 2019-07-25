import csv
import json
import logging
import os
import unittest
from tempfile import gettempdir

from reactome_analysis_api.input_deserializer import create_analysis_input_object

from reactome_analysis_worker import util, geneset_builder
from reactome_analysis_worker.analysers.r_gsva_analyser import ReactomeGSVARAnalyser
from reactome_analysis_worker.models.gene_set_mapping import GeneSetMapping


class TestReactomeRAnalyzer(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        url_logger = logging.getLogger("urllib3")
        url_logger.setLevel(logging.ERROR)

        self.test_json = """
                {
          "analysisId": "test_01",
          "datasets": [
            {
              "data": "\\tSample 1\\tSample2\\tSample 3\\nCD19\\t10\\t20\\t2\\nMS4A1\\t10\\t20\\t2\\n\
              MITF\\t10\\t0\\t0\\n",
              "design": {
                "analysisGroup": [
                  "Treatment",
                  "Control",
                  "Treatment"
                ],
                "comparison": {
                  "group1": "Control",
                  "group2": "Treatment"
                },
                "samples": [
                  "Sample 1",
                  "Sample 2",
                  "Sample 3"
                ],
                "patient": [
                  "Patient 1",
                  "Patient 2",
                  "Patient 3"
               ]
              },
              "name": "First experiment",
              "type": "rnaseq_counts"
            }
          ],
          "methodName": "ssgsea"
        }
        """

    def _get_gene_set(self):
        # load the gene set
        gene_set_file = os.path.join(gettempdir(), "test_gene_set.pkl")

        if os.path.isfile(gene_set_file):
            gene_set = geneset_builder.GeneSet.create_from_file(gene_set_file)
        else:
            gene_set = geneset_builder.fetch_reactome_geneset(source="https://reactome.org/download/current/",
                                                              species="Homo sapiens")
            gene_set.save(gene_set_file)

        return gene_set

    def test_no_design(self):
        test_json = """
                        {
                  "analysisId": "test_01",
                  "datasets": [
                    {
                      "data": "\\tSample 1\\tSample2\\tSample 3\\nCD19\\t10\\t20\\t2\\nMS4A1\\t10\\t20\\t2\\n\
                      MITF\\t10\\t0\\t0\\n",
                      "name": "First experiment",
                      "type": "rnaseq_counts"
                    }
                  ],
                  "methodName": "ssgsea"
                }
                """
        json_obj = json.loads(test_json)
        request = create_analysis_input_object(json_obj)
        request.datasets[0].df = util.string_to_array(request.datasets[0].data)

        self.assertIsNotNone(request)

        # get the mappings
        mappings = util.map_identifiers({"MITF", "CD19", "MS4A1"})

        gene_set = self._get_gene_set()
        gene_id_colname = request.datasets[0].df.dtype.names[0]
        gene_set_mapping = GeneSetMapping.create_mapping(gene_set, identifier_mapping=mappings,
                                                         identifiers=request.datasets[0].df[:][
                                                             gene_id_colname].tolist())

        analyser = ReactomeGSVARAnalyser()
        result = analyser.analyse_request(request=request,
                                          gene_set_mappings={request.datasets[0].name: gene_set_mapping},
                                          identifier_mappings=mappings,
                                          gene_set=gene_set)

        # test the result
        self.assertEqual(1, len(result))
        self.assertIsNotNone(result[0].pathways)
        self.assertIsNotNone(result[0].fold_changes)

        # test the actual result
        reader = csv.DictReader(result[0].pathways.split("\n"), delimiter="\t")
        self.assertEqual(5, len(reader.fieldnames))

        required_fields = ["Pathway", "Sample_1", "Sample2", "Sample_3"]
        for required_field in required_fields:
            self.assertTrue(required_field in reader.fieldnames, "Missing required field " + required_field)

        # test the pathways
        found_pathways = 0

        for pathway in reader:
            found_pathways += 1

            if pathway["Pathway"] == "R-HSA-1280218":
                self.assertEqual("0.5", pathway["Sample_1"])
                self.assertEqual("0.5", pathway["Sample2"])
                self.assertEqual("0.5", pathway["Sample_3"])

            if pathway["Pathway"] == "R-HSA-392499":
                self.assertEqual("-0.5", pathway["Sample_1"])
                self.assertEqual("-0.5", pathway["Sample2"])
                self.assertEqual("-0.5", pathway["Sample_3"])

        self.assertEqual(23, found_pathways)

    def test_ssgsea(self):
        json_obj = json.loads(self.test_json)
        request = create_analysis_input_object(json_obj)
        request.datasets[0].df = util.string_to_array(request.datasets[0].data)

        # get the mappings
        mappings = util.map_identifiers({"MITF", "CD19", "MS4A1"})

        gene_set = self._get_gene_set()
        gene_id_colname = request.datasets[0].df.dtype.names[0]
        gene_set_mapping = GeneSetMapping.create_mapping(gene_set, identifier_mapping=mappings,
                                                         identifiers=request.datasets[0].df[:][
                                                             gene_id_colname].tolist())

        analyser = ReactomeGSVARAnalyser()
        result = analyser.analyse_request(request=request,
                                          gene_set_mappings={request.datasets[0].name: gene_set_mapping},
                                          identifier_mappings=mappings,
                                          gene_set=gene_set)

        # test the result
        self.assertEqual(1, len(result))
        self.assertIsNotNone(result[0].pathways)
        self.assertIsNotNone(result[0].fold_changes)

        # test the actual result
        reader = csv.DictReader(result[0].pathways.split("\n"), delimiter="\t")
        self.assertEqual(5, len(reader.fieldnames))

        required_fields = ["Pathway", "Name", "Sample.1", "Sample.2", "Sample.3"]
        for required_field in required_fields:
            self.assertTrue(required_field in reader.fieldnames)

        # test the pathways
        found_pathways = 0
        found_p1 = False
        found_p2 = False

        for pathway in reader:
            found_pathways += 1

            if pathway["Pathway"] == "R-HSA-1280218":
                self.assertEqual("0.5", pathway["Sample.1"])
                self.assertEqual("0.5", pathway["Sample.2"])
                self.assertEqual("0.5", pathway["Sample.3"])
                found_p1 = True

            if pathway["Pathway"] == "R-HSA-392499":
                self.assertEqual("-0.5", pathway["Sample.1"])
                self.assertEqual("-0.5", pathway["Sample.2"])
                self.assertEqual("-0.5", pathway["Sample.3"])
                found_p2 = True

        self.assertTrue(found_p1)
        self.assertTrue(found_p2)
        self.assertEqual(23, found_pathways)

    def test_pathway_string(self):
        json_obj = json.loads(self.test_json)

        # add the parameters
        json_obj["parameters"] = [{"name": "pathways", "value": "R-HSA-1280218,R-HSA-392499"},
                                  {"name": "create_reactome_visualization", "value": "False"}]

        request = create_analysis_input_object(json_obj)
        request.datasets[0].df = util.string_to_array(request.datasets[0].data)

        # get the mappings
        mappings = util.map_identifiers({"MITF", "CD19", "MS4A1"})

        gene_set = self._get_gene_set()
        gene_id_colname = request.datasets[0].df.dtype.names[0]
        gene_set_mapping = GeneSetMapping.create_mapping(gene_set, identifier_mapping=mappings,
                                                         identifiers=request.datasets[0].df[:][
                                                             gene_id_colname].tolist())

        analyser = ReactomeGSVARAnalyser()
        result = analyser.analyse_request(request=request,
                                          gene_set_mappings={request.datasets[0].name: gene_set_mapping},
                                          identifier_mappings=mappings,
                                          gene_set=gene_set)

        # test the result
        self.assertEqual(1, len(result))
        self.assertIsNotNone(result[0].pathways)
        self.assertIsNotNone(result[0].fold_changes)

        # test the actual result
        reader = csv.DictReader(result[0].pathways.split("\n"), delimiter="\t")
        self.assertEqual(5, len(reader.fieldnames))

        # there should only be two entries
        n_entries = 0

        for line in reader:
            n_entries += 1

        self.assertEqual(2, n_entries)
