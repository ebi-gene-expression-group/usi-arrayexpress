import unittest
from utils.converter_utils import get_taxon


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


if __name__ == '__main__':
    unittest.main()
