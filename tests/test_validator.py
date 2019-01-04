import os
import unittest
import logging

from converter import datamodel
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
        self.logger = logging.getLogger()
        self.logger.setLevel(10)

    def test_protocol_validation(self):
        # Everything should be fine, no errors expected
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertEqual(error_codes, [])

        # Creating duplicated protocol name, should give error PROT-E04
        self.sub.protocol[0].alias = "P-MTAB-48204"
        # Protocol type not from controlled vocabulary, should give error PROT-E05
        self.sub.protocol[1].protocol_type.value = "nonexiting protocol"
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertIn("PROT-E04", error_codes)
        self.assertIn("PROT-E05", error_codes)

        # Empty alias should give error PROT-E02
        self.sub.protocol[0].alias = ""
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertIn("PROT-E02", error_codes)

        # Deleting all protocols should only give error PROT-E01
        self.sub.protocol = []
        error_codes = metadata_validation.run_protocol_checks(self.sub, self.logger)
        self.assertEqual(["PROT-E01"], error_codes)

    def test_sample_validation(self):
        # Everything should be fine, no errors expected
        error_codes = metadata_validation.run_sample_checks(self.sub, self.logger)
        self.assertEqual(error_codes, [])

        # Sample with no name should give error SAMP-E02
        self.sub.sample[1].alias = ""

        # Empty organism, should give error SAMP-E03
        self.sub.sample[2].taxon = ""

        # Organism not from taxonomy, should give error SAMP-E08
        self.sub.sample[0].taxon = "unknown species"

        # Non-annotated factor should give error SAMP-E05
        self.sub.study.experimental_factor.append(datamodel.Attribute("forgotten_factor", None, None))

        # Add nonsense unit, should give error SAMP-E04
        self.sub.sample[0].attributes["organism"].unit = datamodel.Unit("xxx", "silly unit", None, None)
        error_codes = metadata_validation.run_sample_checks(self.sub, self.logger)
        self.assertIn("SAMP-E02", error_codes)
        self.assertIn("SAMP-E03", error_codes)
        self.assertIn("SAMP-E04", error_codes)
        self.assertIn("SAMP-E05", error_codes)
        self.assertIn("SAMP-E08", error_codes)

        # Removing all samples, should only give error "SAMP-E01"
        self.sub.sample = []
        error_codes = metadata_validation.run_sample_checks(self.sub, self.logger)
        self.assertEqual(["SAMP-E01"], error_codes)

    def test_study_validation(self):
        error_codes = metadata_validation.run_study_checks(self.sub, self.logger)
        self.assertEqual(error_codes, [])
        # No title should give error STUD-E01
        self.sub.study.title = ""
        # No description should give error STUD-E02
        self.sub.study.description = ""
        # Non-existent experiment type should give error STUD-E04
        self.sub.study.experiment_type[0] = "new experiment type"
        # Wrong date format should give error STUD-E05
        self.sub.study.date_of_experiment = "26/08/2011"
        error_codes = metadata_validation.run_study_checks(self.sub, self.logger)
        self.assertIn("STUD-E01", error_codes)
        self.assertIn("STUD-E02", error_codes)
        self.assertIn("STUD-E04", error_codes)
        self.assertIn("STUD-E05", error_codes)

        # Remove experiment types should give STUD-E03
        self.sub.study.experiment_type = []
        error_codes = metadata_validation.run_study_checks(self.sub, self.logger)
        self.assertIn("STUD-E03", error_codes)


if __name__ == '__main__':
    unittest.main()

