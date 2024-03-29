import csv
import json
import logging
import os
import unittest
import time
from tempfile import gettempdir
import rpy2.rinterface as ri
import rpy2.robjects as ro
import numpy
from io import StringIO

from reactome_analysis_api.input_deserializer import create_analysis_input_object

from reactome_analysis_worker import util, geneset_builder
from reactome_analysis_worker.analysers.r_analyser import ReactomeRAnalyser
from reactome_analysis_worker.models.gene_set_mapping import GeneSetMapping
from reactome_analysis_worker.reactome_analysis_worker import ReactomeAnalysisWorker


ri.initr()


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
          "methodName": "camera",
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

    def test_parameter_passing(self):
        json_obj = json.loads(self.test_json)
        json_obj["parameters"].append({"name": "max_missing_values", "value": "1"})

        # remove the patient since this coefficient cannot be estimated
        json_obj["datasets"][0]["design"].pop("patient")

        request = create_analysis_input_object(json_obj)
        request.datasets[0].df = util.string_to_array(request.datasets[0].data)

        self.assertEqual(3, len(request.parameters))
        # default values inserted automatically
        self.assertEqual(6, len(request.parameter_dict))
        self.assertTrue("max_missing_values" in request.parameter_dict)

        # get the mappings
        mappings = util.map_identifiers({"MITF", "CD19", "MS4A1"})

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
        self.assertEqual(24, len(result_lines))

        reader = csv.DictReader(result_lines, delimiter="\t")
        required_fields = ("Pathway", "Name", "Direction", "FDR", "PValue", "NGenes")
        for field in required_fields:
            self.assertTrue(field in reader.fieldnames, "Missing field " + field)

        pathways_up = ("R-HSA-392499", "R-HSA-597592", "R-HSA-2990846", "R-HSA-3108232", "R-HSA-3232118")
        for row in reader:
            if reader.line_num == 2:
                self.assertTrue(row["Pathway"] == "R-HSA-392499")
            if reader.line_num == 6:
                self.assertTrue(row["Pathway"] == "R-HSA-3232118")
            if reader.line_num == 15:
                self.assertTrue(row["Pathway"] == "R-HSA-162582")
            if reader.line_num == 24:
                self.assertTrue(row["Pathway"] == "R-HSA-6811558")

            if row["Pathway"] in pathways_up:
                self.assertTrue(row["Direction"] == "Down")
                self.assertTrue(float(row["av_foldchange"]) < 0, "Incorrect regulation for " + row["Pathway"])
            else:
                self.assertTrue(row["Direction"] == "Up")
                self.assertTrue(float(row["av_foldchange"]) > 0)

        # test the FC result
        self.assertIsNotNone(result[0].fold_changes)
        fc_lines = result[0].fold_changes.split("\n")
        self.assertEqual(4, len(fc_lines))

        fc_reader = csv.DictReader(fc_lines, delimiter="\t")
        fc_fields = ("logFC", "Identifier")

        for field in fc_fields:
            self.assertTrue(field in fc_reader.fieldnames, "Missing FC field " + field)

        mitf_found = False

        for row in fc_reader:
            if row["Identifier"] == "MITF":
                self.assertAlmostEqual(4.53, float(row["logFC"]), delta=0.01)
                mitf_found = True

        self.assertTrue(mitf_found, "Failed to find MITF in FC data")

    def update_heartbeat(self):
      self.last_heartbeat = int(time.time())

    def test_heartbeat(self):
      json_obj = json.loads(self.test_json)
      json_obj["parameters"].append({"name": "max_missing_values", "value": "1"})

      # remove the patient since this coefficient cannot be estimated
      json_obj["datasets"][0]["design"].pop("patient")

      request = create_analysis_input_object(json_obj)
      request.datasets[0].df = util.string_to_array(request.datasets[0].data)

      # get the mappings
      mappings = util.map_identifiers({"MITF", "CD19", "MS4A1"})

      # filter the dataset
      request.datasets[0].df = ReactomeAnalysisWorker._filter_dataset(request.datasets[0].df, mappings,
                                                                      request.datasets[0].design, 1)

      gene_set = self._get_gene_set()
      gene_id_colname = request.datasets[0].df.dtype.names[0]
      gene_set_mapping = GeneSetMapping.create_mapping(gene_set, identifier_mapping=mappings,
                                                        identifiers=request.datasets[0].df[:][gene_id_colname].tolist())

      analyser = ReactomeRAnalyser()
      analyser.set_heartbeat_callback(self.update_heartbeat)
      start_time = int(time.time()) - 1

      result = analyser.analyse_request(request=request,
                                        gene_set_mappings={request.datasets[0].name: gene_set_mapping},
                                        identifier_mappings=mappings,
                                        gene_set=gene_set)

      # make sure the heartbeat was updated
      self.assertGreater(self.last_heartbeat, start_time)

    def test_df_to_str(self):
      # create the test data.frame
      ro.reval("""
          test_frame = data.frame(
              name = c("John", "Doe"),
              age = c(1, 2),
              row.names = c("Id1", "Id2")
          )

          test_frame_2 = data.frame(
              name = c("John", "Doe"),
              age = c(1.12345, 2.12345),
              row.names = c("Id1", "Id2")
          )
      """)

      r_data_frame = ri.globalenv["test_frame"]

      string_df = ReactomeRAnalyser.data_frame_to_string(r_data_frame)

      self.assertIsNotNone(string_df)
      self.assertEqual("\\tname\\tage\\nId1\\tJohn\\t1.0\\nId2\\tDoe\\t2.0", string_df)

      # check the precision
      r_data_frame2 = ri.globalenv["test_frame_2"]

      string_df2 = ReactomeRAnalyser.data_frame_to_string(r_data_frame2)

      self.assertIsNotNone(string_df2)
      self.assertEqual("\\tname\\tage\\nId1\\tJohn\\t1.12345\\nId2\\tDoe\\t2.12345", string_df2)

    def test_convert_dataset(self):
      # create the sample data
      df = numpy.genfromtxt(StringIO("Gene\tS1\tS2\tS3\nHLA-DQB1\t1\t2\t3\nMITF\t4\t5\t6\n"), dtype=None, delimiter="\t", names=True, 
                            encoding = None)

      class MockDS:
        pass

      dataset = MockDS()
      dataset.df = df

      sample_names = ri.StrSexpVector(["S1", "S2", "S3"])

      # test the function
      converted = ReactomeRAnalyser._convert_dataset(dataset, sample_names)

      self.assertIsNotNone(converted)

      rownames = ri.baseenv["rownames"]

      converted_genes = list(rownames(converted))

      self.assertIsNotNone(converted_genes)
      self.assertEquals(["HLA-DQB1", "MITF"], converted_genes)
