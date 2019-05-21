import json
import unittest
from collections import OrderedDict

from converter.datamodel import Sample, Attribute, Unit


class TestJSONparsing(unittest.TestCase):

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
                      "genotype": [
                          {"value": "wild type genotype"}]
                  }}
        target = Sample("Control 1", None, "Homo sapiens", 9606,
                        OrderedDict([
                            ("organism", Attribute("Homo sapiens", None, None, None)),
                            ("disease", Attribute("acute myeloid leukemia", None, None, None)),
                            ("cell line", Attribute("KG1", None, None, None)),
                            ("genotype", Attribute("wild type genotype", None, None, None))]),
                        None, None)
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
                            ("cell line", Attribute("KG1", None, "EFO_0002218", None))]),
                        None, None)
        converted_sample = Sample.from_json(source)
        self.assertEqual(str(converted_sample), str(target))

