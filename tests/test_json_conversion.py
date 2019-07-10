
import unittest
from collections import OrderedDict

from converter.datamodel import Sample, Attribute, Unit


class TestSampleParsing(unittest.TestCase):

    def test_generate_simple_sample_object(self):
        source = {"alias": "Control 1",
                  "taxon": "Homo sapiens",
                  "taxonId": 9606,
                  "attributes": {
                      "organism": [
                          {"value": "Homo sapiens"}
                      ],
                      "disease": [
                          {"value": "acute myeloid leukemia"}
                      ],
                      "cell line": [
                          {"value": "KG1"}
                      ],
                      "material type": [
                          {"value": "cell"}]
                  }}
        target = Sample("Control 1", None, "Homo sapiens", 9606,
                        OrderedDict([
                            ("organism", Attribute("Homo sapiens", None, None, None)),
                            ("disease", Attribute("acute myeloid leukemia", None, None, None)),
                            ("cell line", Attribute("KG1", None, None, None)),
                            ]), "cell", None)
        converted_sample = Sample.from_json(source)
        self.assertEqual(str(converted_sample), str(target))

    def test_generate_sample_object_with_term_and_unit(self):
        source = {"alias": "Control 1",
                  "taxon": "Homo sapiens",
                  "taxonId": 9606,
                  "attributes": {
                      "organism": [
                          {"value": "Homo sapiens"}
                      ],
                      "age": [
                          {"value": "12",
                           "units": "year"}
                      ],
                      "cell line": [
                          {"value": "KG1",
                           "terms": [{"url": "http://www.ebi.ac.uk/efo/EFO_0002218"}]}
                      ]
                  }}
        target = Sample("Control 1", None, "Homo sapiens", 9606,
                        OrderedDict([
                            ("organism", Attribute("Homo sapiens", None, None, None)),
                            ("age", Attribute("12", Unit("year", None, None, None), None, None)),
                            ("cell line", Attribute("KG1", None, "EFO_0002218", "EFO"))]),
                        None, None)
        converted_sample = Sample.from_json(source)
        self.assertEqual(str(converted_sample), str(target))

    def test_generate_sample_object_with_non_EFO_term(self):
        source = {"alias": "Control 1",
                  "taxon": "Homo sapiens",
                  "taxonId": 9606,
                  "attributes": {
                      "compound": [
                          {"value": "sodium arsanilate",
                           "terms": [{"url": "http://purl.obolibrary.org/obo/CHEBI_36049"}]}
                      ]
                  }}
        target = Sample("Control 1", None, "Homo sapiens", 9606,
                        OrderedDict([
                            ("compound", Attribute("sodium arsanilate", None, "CHEBI_36049", "CHEBI"))]),
                        None, None)
        converted_sample = Sample.from_json(source)
        self.assertEqual(str(converted_sample), str(target))


class TestStudyParsing(unittest.TestCase):

    def test_generate_simple_study(self):
        pass
