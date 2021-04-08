import os
import unittest

from reactome_analysis_worker import util


class TestUtil(unittest.TestCase):
    def testConversion(self):
        """
        Test basic data structure conversion
        :return:
        """
        text = "\\tSample 1\\tSample2\\nGene 1\\t10\\t20\\nGene 2\\t10\\t30\\nGene 3\\t10\\t30\\n"

        array = util.string_to_array(text)
        self.assertEqual(1, array.ndim)
        self.assertEqual("Gene 1", array[0][array.dtype.names[0]])
        self.assertEqual("Gene 2", array[1][array.dtype.names[0]])

    def testMissingValuesConversion(self):
        """
        Test basic data structure conversion
        :return:
        """
        text = "\\tSample 1\\tSample2\\nGene 1\\t10\\t20\\nGene 2\\tNA\\t30\\nGene 3\\t10\\t30\\n"

        array = util.string_to_array(text)
        self.assertEqual(1, array.ndim)
        
        if any("U" in dt_type for (_, dt_type) in array.dtype.descr[1:]):
            self.fail("Non-numeric columns")

        self.assertEqual(0, array[1][1])

    def testOneLineConversion(self):
        """
        Test basic data structure conversion
        :return:
        """
        text = "\\tSample 1\\tSample2\\nGene 1\\t10\\t20\\n"

        array = util.string_to_array(text)
        self.assertEqual(0, array.ndim)
        self.assertEqual("Gene 1", array[array.dtype.names[0]].item())

    def testMapping(self):
        # use all genes from one pathway
        genes = set()
        with open(os.path.join(os.path.dirname(__file__), "testfiles", "R-HSA-1980143.uniprot.txt")) as reader:
            for line in reader:
                genes.add(line.strip())

        mapped_identifiers = util.map_identifiers(genes, return_all=True)

        self.assertEqual(len(genes), len(mapped_identifiers))

        # all identifiers should only map to a single one
        for mapped_identifier in mapped_identifiers.values():
            # There are multiple mappings when referring to isoforms (one case)
            self.assertTrue(len(mapped_identifier) < 3, msg="Multiple mappings for {}".format(
                ",".join(mapped_identifier)))

    def test_interactor_mapping(self):
        mapped_identifier = util.map_identifiers(["MS4A1"], return_all=True)

        self.assertEqual(1, len(mapped_identifier))
        self.assertEqual(1, len(mapped_identifier["MS4A1"]))
        self.assertEqual("P11836", mapped_identifier["MS4A1"][0])

    def test_failed_mapping(self):
        genes = {"Gene1", "Gene2"}
        mapped_identifiers = util.map_identifiers(genes, return_all=True)

        self.assertEqual(0, len(mapped_identifiers))

        # test with a single gene
        genes = {"Gene1"}
        mapped_identifiers = util.map_identifiers(genes, return_all=True)

        self.assertEqual(0, len(mapped_identifiers))
