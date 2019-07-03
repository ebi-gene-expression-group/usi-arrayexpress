"""Tests for the mapping config guided USI-JSON to data model conversion."""

import unittest
from collections import OrderedDict

from converter.datamodel import Sample, Attribute, Unit
from converter.json2dm import JSONConverter


class TestSampleParsing(unittest.TestCase):

    def setUp(self):
        # Creating mapping instructions. Using Python dictionary directly instead of JSON file
        self.mapping = {
            "sample": {
                "alias": {"type": "string", "import": {"test": {"path": ["alias"], "method": "import_string"}}},
                "taxon": {"type": "string", "import": {"test": {"path": ["taxon"], "method": "import_string"}}},
                "taxonId": {"type": "string", "import": {"test": {"path": ["taxonId"], "method": "import_string"}}},
                "material_type": {"type": "string", "import": {"test": {"path": ["attributes", "material_type"],
                                                                        "method": "generate_attribute_from_json"}}},
                "description": {"type": "string", "import": {"test": {"path": ["attributes", "description"],
                                                                      "method": "generate_attribute_from_json"}}},
                "attributes": {"type": "array", "import": {"test": {"path": ["attributes"],
                                                                    "method": "generate_sample_attribute_dict"}}}}}
        # Create the converter
        self.converter = JSONConverter(self.mapping, import_key="test")

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
                      "material_type": [
                          {"value": "cell"}]
                  }}
        target = Sample("Control 1", None, "Homo sapiens", 9606,
                        OrderedDict([
                            ("organism", Attribute("Homo sapiens", None, None, None)),
                            ("disease", Attribute("acute myeloid leukemia", None, None, None)),
                            ("cell line", Attribute("KG1", None, None, None)),
                            ]), "cell", None)
        sample_dict = self.converter.convert_submittable(source, "sample")
        converted_sample = Sample.from_dict(sample_dict)
        self.assertEqual(str(target), str(converted_sample))

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
                            ("age", Attribute("12", Unit("year", "time unit", None, None), None, None)),
                            ("cell line", Attribute("KG1", None, "EFO_0002218", "EFO"))]),
                        None, None)
        sample_dict = self.converter.convert_submittable(source, "sample")
        converted_sample = Sample.from_dict(sample_dict)
        self.assertEqual(str(target), str(converted_sample))

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
        sample_dict = self.converter.convert_submittable(source, "sample")
        converted_sample = Sample.from_dict(sample_dict)
        self.assertEqual(str(target), str(converted_sample))


if __name__ == '__main__':
    unittest.main()
