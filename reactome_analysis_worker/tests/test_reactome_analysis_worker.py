#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `reactome_analysis_worker` package."""


import json
import logging
import os
import unittest
from tempfile import gettempdir

from reactome_analysis_api.input_deserializer import create_analysis_input_object
from reactome_analysis_api.models.analysis_result import AnalysisResult
from reactome_analysis_utils import reactome_mq
from reactome_analysis_utils import reactome_storage

from reactome_analysis_worker import geneset_builder
from reactome_analysis_worker import reactome_analysis_worker
from reactome_analysis_worker import util


class TestReactome_analysis_worker(unittest.TestCase):
    """Tests for `reactome_analysis_worker` package."""

    def setUp(self):
        os.environ["REDIS_HOST"] = "192.168.99.100"
        os.environ["REDIS_PORT"] = "31297"
        os.environ["REDIS_PASSWORD"] = "test"
        os.environ["RABBIT_HOST"] = "192.168.99.100"
        os.environ["RABBIT_PORT"] = "31715"
        os.environ["RABBIT_USER"] = "test"
        os.environ["RABBIT_PASSWORD"] = "test"

        logging.basicConfig(level=logging.DEBUG)
        pika_logger = logging.getLogger("pika")
        pika_logger.setLevel(logging.ERROR)

        util_logger = logging.getLogger("reactome_analysis_utils")
        util_logger.setLevel(logging.DEBUG)

    def test_process_all_messages(self):
        if (True):
            return

        if not os.path.isfile("/tmp/reactome_homo_sapiens.pkl"):
            gene_set = self._get_gene_set()
            gene_set.save("/tmp/reactome_homo_sapiens.pkl")
            gene_set.save("/tmp/reactome_homo_sapiens_interactors.pkl")
        worker = reactome_analysis_worker.ReactomeAnalysisWorker()
        worker.start_analyses()

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

    def test_additional_parameter_analysis(self):
        """
        Tests an analysis request that contains multiple blocking factors
        """
        request_json = """
                {
          "analysisId": "test_01",
          "datasets": [
            {
              "data": "\\tSample 1\\tSample2\\tSample 3\\nCD19\\t10\\t2\\t20\\nCD20\\t10\\t20\\t2\\nMITF\\t40\\t20\\t10\\n",
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
                ]
              },
              "name": "First experiment",
              "type": "rnaseq_counts"
            }
          ],
          "methodName": "camera",
          "parametes": [
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

        # make sure the JSON is valid
        obj = json.loads(request_json)
        self.assertIsNotNone(obj)

        # submit the request
        mq = reactome_mq.ReactomeMQ()
        mq.post_analysis(request_json, "camera")

        # download the gene sets
        if not os.path.isfile("/tmp/reactome_homo_sapiens.pkl"):
            geneset = self._get_gene_set()
            geneset.save("/tmp/reactome_homo_sapiens.pkl")

        # enable debug mode
        os.environ["REACTOME_WORKER_DEBUG"] = "True"

        # start to listen to analyses
        worker = reactome_analysis_worker.ReactomeAnalysisWorker()
        worker.process_single_message()

        # fetch the result
        storage = reactome_storage.ReactomeStorage()
        result_text = storage.get_result("test_01")

        self.assertIsNotNone(result_text, "Result was not saved in redis")
        json_obj = json.loads(result_text)
        result = AnalysisResult.from_dict(json_obj)

        self.assertIsNotNone(result)
        self.assertIsNotNone(result.mappings)
        self.assertIsNotNone(result.results)
        self.assertEqual("68", result.release)

        self.assertEqual(1, len(result.results))
        self.assertIsNotNone(result.results[0].pathways)
        self.assertIsNotNone(result.results[0].fold_changes)

        pathway_lines = result.results[0].pathways.split("\n")
        self.assertEqual(23, len(pathway_lines))

        gene_lines = result.results[0].fold_changes.split("\n")
        self.assertEqual(4, len(gene_lines))

    def test_incorrect_genes(self):
        json_request = """
        {
  "analysisId": "test_02",
  "datasets": [
    {
      "data": "\\tSample 1\\tSample2\\tSample3\\nCD19\\t10\\t20\\t5\\nCD20\\t10\\t49\\t2\\nMITF\\t56\\t14\\t24\\n",
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
        ]
      },
      "name": "First experiment",
      "type": "rnaseq_counts"
    },
    {
      "data": "\\tSample 1\\tSample2\\tGene 1\\t10\\t20\\n",
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
        ]
      },
      "name": "First experiment",
      "type": "rnaseq_counts"
    }
  ],
      "methodName": "Camera",
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

        # make sure the JSON is valid
        obj = json.loads(json_request)
        self.assertIsNotNone(obj)

        # submit the request
        mq = reactome_mq.ReactomeMQ()
        mq.post_analysis(json_request, "camera")

        # download the gene sets
        if not os.path.isfile("/tmp/reactome_homo_sapiens.pkl"):
            geneset = self._get_gene_set()
            geneset.save("/tmp/reactome_homo_sapiens.pkl")

        # enable debug mode
        os.environ["REACTOME_WORKER_DEBUG"] = "True"

        # start to listen to analyses
        worker = reactome_analysis_worker.ReactomeAnalysisWorker()
        worker.process_single_message()

        # fetch the result
        storage = reactome_storage.ReactomeStorage()

        status = storage.get_status("test_02")
        status_obj = json.loads(status)
        self.assertEqual("failed", status_obj["status"])

        result_text = storage.get_result("test_02")
        self.assertIsNone(result_text)

    def test_no_design_filtering(self):
        test_json = """
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

        worker = reactome_analysis_worker.ReactomeAnalysisWorker()

        json_obj = json.loads(test_json)
        request_obj = create_analysis_input_object(json_obj)
        worker._convert_datasets(request_obj)
        mappings = util.map_identifiers({"MITF", "CD19", "MS4A1"})

        self.assertEqual(3, len(request_obj.datasets[0].df))

        filtered_df = reactome_analysis_worker.ReactomeAnalysisWorker._filter_dataset(request_obj.datasets[0].df,
                                                                                      mappings,
                                                                                      None,
                                                                                      0.5)

        self.assertIsNotNone(filtered_df)
        self.assertEqual(2, len(filtered_df))
