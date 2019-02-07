import unittest
import os

from utils.converter_utils import get_taxon
from utils.common_utils import get_term_descendants, create_logger


class TestTaxonRetrieval(unittest.TestCase):

    def test_human(self):
        organism = "Homo sapiens"
        taxon_id = get_taxon(organism)
        self.assertEqual(taxon_id, 9606)

    def test_mixed_species(self):
        mixed_orgs = "Human and mouse"
        taxon_id = get_taxon(mixed_orgs)
        self.assertEqual(taxon_id, 1427524)

    def test_unknown(self):
        print(get_taxon("xxxx"))


class TestOntologyTermRetrieval(unittest.TestCase):
    def setUp(self):
        wd = os.path.dirname(os.path.realpath(__file__))
        # Add logger to hold error messages
        self.logger = create_logger(wd, "testing", "validation", 10)

    def test_get_unit_children(self):
        uri = "http://purl.obolibrary.org/obo/UO_0000000"
        ontology = "EFO"
        terms = get_term_descendants(ontology, uri, self.logger)
        self.assertIn("colony forming unit", terms)


if __name__ == '__main__':
    unittest.main()