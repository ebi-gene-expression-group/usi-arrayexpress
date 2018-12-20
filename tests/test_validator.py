import os
import unittest

from converter.converting import data_objects_from_magetab
from validator import metadata_validation
from utils.common_utils import create_logger


class TestMetaDataValidation(unittest.TestCase):

    def setUp(self):
        # Loading the MAGE-TAB and tranforming it into the common datamodel
        wd = os.path.dirname(os.path.realpath(__file__))
        idf = os.path.join(wd, 'test_data', 'E-MTAB-4250.idf.txt')
        sdrf = os.path.join(wd, 'test_data', 'E-MTAB-4250.sdrf.txt')
        self.sub = data_objects_from_magetab(idf, sdrf)
        # Add logger to hold error messages
        self.logger = create_logger(wd, "testing", "validation", 10)

    def test_protocol_validation(self):
        # Everything should be fine, no errors expected
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertEqual(error_codes, [])

        # Creating duplicated protocol name, should give error PROT-E04
        self.sub.protocol[0].alias = "P-MTAB-48204"
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertIn("PROT-E04", error_codes)

        # Empty alias should give error PROT-E02
        self.sub.protocol[0].alias = ""
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertIn("PROT-E02", error_codes)

        # Deleting all protocols should give error PROT-E01
        self.sub.protocol = []
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertIn("PROT-E01", error_codes)


if __name__ == '__main__':
    unittest.main()

