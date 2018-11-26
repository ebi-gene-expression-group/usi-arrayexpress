import sys
import unittest


from converter.converting import get_taxon


class TestTaxonRetrieval(unittest.TestCase):

    def test_human_and_mixed(self):
        organism = "Homo sapiens"
        taxon_id = get_taxon(organism)
        self.assertEqual(taxon_id, 9606)
        mixed_orgs = "Human and mouse"
        taxon_id = get_taxon(mixed_orgs)
        self.assertEqual(taxon_id, 1427524)
        print(get_taxon("xxxx"))

if __name__ == '__main__':
    unittest.main()