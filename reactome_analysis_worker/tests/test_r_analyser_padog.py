import csv
import json
import logging
import os
import unittest
from tempfile import gettempdir

from reactome_analysis_api.input_deserializer import create_analysis_input_object

from reactome_analysis_worker import util, geneset_builder
from reactome_analysis_worker.analysers.r_analyser import ReactomeRAnalyser
from reactome_analysis_worker.models.gene_set_mapping import GeneSetMapping
from reactome_analysis_worker.reactome_analysis_worker import ReactomeAnalysisWorker


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
              "data": "\\tSample 1\\tSample2\\tSample 3\\tSample 4\\tSample 5\\tSample 6\\tSample 7\\n\
              CD19\\t10\\t20\\t2\\t1\\t7\\t20\\t1\\n\
              MS4A1\\t10\\t20\\t2\\t50\\t9\\t24\\t1\\n\
              MITF\\t10\\t0\\t0\\t20\\t1\\t10\\t1\\n\
              SDC1\\t5\\t10\\t2\\t8\\t5\\t18\\t1\\n\
              CD38\\t10\\t2\\t7\\t0\\t10\\t9\\t1\\n\
              EGFR\\t10\\t2\\t7\\t0\\t10\\t9\\t1\\n\
              IL10\\t10\\t2\\t7\\t0\\t10\\t9\\t1\\n\
              IL6\\t10\\t2\\t7\\t0\\t10\\t9\\t1\\n\
              GRB2\\t10\\t2\\t7\\t0\\t10\\t9\\t1\\n\
              GAB1\\t10\\t2\\t7\\t0\\t10\\t9\\t1\\n\
              SHC1\\t10\\t2\\t7\\t0\\t10\\t9\\t1\\n",
              "design": {
                "analysisGroup": [
                  "Treatment",
                  "Control",
                  "Treatment",
                  "Control",
                  "Control",
                  "Treatment",
                  "Other"
                ],
                "comparison": {
                  "group1": "Control",
                  "group2": "Treatment"
                },
                "samples": [
                  "Sample 1",
                  "Sample 2",
                  "Sample 3",
                  "Sample 4",
                  "Sample 5",
                  "Sample 6",
                  "Sample 7"
                ]
              },
              "name": "First experiment",
              "type": "rnaseq_counts"
            }
          ],
          "methodName": "padog",
          "parameters": [
            {
              "name": "permutations",
              "value": "10"
            },
            {
              "name": "permutations",
              "value": "10"
            }
          ]
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

    def test_analysis(self):
        json_obj = json.loads(self.test_json)
        request = create_analysis_input_object(json_obj)
        request.datasets[0].df = util.string_to_array(request.datasets[0].data)

        # get the mappings
        mappings = util.map_identifiers({"MITF", "CD19", "MS4A1", "SDC1", "CD38", "EGFR", "IL10", "IL6", "GRB2", "GAB1", "SHC1"})

        # filter the dataset
        request.datasets[0].df = ReactomeAnalysisWorker._filter_dataset(request.datasets[0].df, mappings,
                                                                        request.datasets[0].design, 1)

        gene_set = self._get_gene_set()
        gene_id_colname = request.datasets[0].df.dtype.names[0]
        gene_set_mapping = GeneSetMapping.create_mapping(gene_set, identifier_mapping=mappings,
                                                         identifiers=request.datasets[0].df[:][gene_id_colname].tolist())

        analyser = ReactomeRAnalyser()
        result = analyser.analyse_request(request=request,
                                          gene_set_mappings={request.datasets[0].name: gene_set_mapping},
                                          identifier_mappings=mappings,
                                          gene_set=gene_set)

        # test the result
        self.assertEqual(1, len(result))
        self.assertIsNotNone(result[0].pathways)

        result_lines = result[0].pathways.split("\n")
        self.assertEqual(233, len(result_lines))

        reader = csv.DictReader(result_lines, delimiter="\t")
        required_fields = ("Pathway", "Name", "Direction", "FDR", "PValue", "NGenes")
        for field in required_fields:
            self.assertTrue(field in reader.fieldnames, "Missing field " + field)
