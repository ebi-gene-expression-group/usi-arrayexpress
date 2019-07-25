
import unittest
from collections import OrderedDict

from converter.datamodel.sample import Sample
from converter.datamodel.components import Attribute, Unit


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
                      "material_type": [
                          {"value": "cell"}]
                  }}
        target = Sample(alias="Control 1", accession=None, taxon="Homo sapiens", taxonId=9606,
                        attributes=OrderedDict([
                            ("organism",
                             Attribute(value="Homo sapiens", unit=None, term_accession=None, term_source=None)),
                            ("disease",
                             Attribute(value="acute myeloid leukemia", unit=None, term_accession=None, term_source=None)),
                            ("cell line",
                             Attribute(value="KG1", unit=None, term_accession=None, term_source=None))]),
                        material_type="cell", description=None)
        converted_sample = Sample(**source)
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
        target = Sample(alias="Control 1", accession=None, taxon="Homo sapiens", taxonId=9606,
                        attributes=OrderedDict([
                            ("organism", Attribute(value="Homo sapiens",unit=None,
                                                   term_accession=None, term_source=None)),
                            ("age", Attribute(value="12",
                                              unit=Unit(value="year", unit_type=None,
                                                        term_source=None, term_accession=None),
                                              term_source=None, term_accession=None)),
                            ("cell line", Attribute(value="KG1", unit=None,
                                                    term_accession="EFO_0002218", term_source="EFO"))]),
                        material_type=None, description=None)
        converted_sample = Sample(**source)
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
        target = Sample(alias="Control 1", accession=None, taxon="Homo sapiens", taxonId=9606,
                        attributes=OrderedDict([
                            ("compound", Attribute(value="sodium arsanilate", unit=None,
                                                   term_accession="CHEBI_36049", term_source="CHEBI"))]),
                        material_type=None, description=None)
        converted_sample = Sample(**source)
        self.assertEqual(str(converted_sample), str(target))


class TestStudyParsing(unittest.TestCase):

    def test_generate_simple_study(self):
        pass
