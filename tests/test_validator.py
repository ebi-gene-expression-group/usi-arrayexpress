import os
import unittest
import logging

import datamodel
from converter.magetab2dm import data_objects_from_magetab
from validator import metadata_validation
from utils.converter_utils import guess_submission_type_from_sdrf, guess_submission_type_from_idf, \
    read_sdrf_file, read_idf_file


class TestMetaDataValidation(unittest.TestCase):

    def setUp(self):
        # Loading the MAGE-TAB and tranforming it into the common datamodel
        wd = os.path.dirname(os.path.realpath(__file__))
        idf = os.path.join(wd, 'test_data', 'E-MTAB-4250.idf.txt')
        sdrf = os.path.join(wd, 'test_data', 'E-MTAB-4250.sdrf.txt')
        idf_dict = read_idf_file(idf)
        sdrf_data, header, header_dict = read_sdrf_file(sdrf)
        submission_type = guess_submission_type_from_sdrf(sdrf_data, header, header_dict)
        if not submission_type:
            submission_type = guess_submission_type_from_idf(idf_dict)
        self.sub = data_objects_from_magetab(idf, sdrf, submission_type)
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
        self.sub.study.experimental_factor.append(datamodel.components.Attribute(value="forgotten_factor"))
        # Add nonsense unit, should give error SAMP-E04
        self.sub.sample[0].attributes["organism"].unit = datamodel.components.Unit(value="xxx", unit_type="silly unit")
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

    def test_project_validation(self):
        # All should be fine
        error_codes = metadata_validation.run_project_checks(self.sub, self.logger)
        self.assertEqual(error_codes, [])

        # Contact with no last name should give error PROJ-E02
        self.sub.project.contacts[0].lastName = ""
        # False role (and no submitter) should give error PROJ-E03 and PROJ-05
        self.sub.project.contacts[0].roles = ["fake"]
        # No email should give PROJ-E-04 (no submitter details found)
        self.sub.project.contacts[0].email = ""
        # Wrong PubMed ID and DOI should give PROJ-E06 and PROJ-E07
        self.sub.project.publications = [datamodel.components.Publication(articleTitle="Test Article Title",
                                                                          authors="everyone", pubmedId="haha",
                                                                          doi="000")]
        # Wrong release date format should give error PROJ-E09
        self.sub.project.releaseDate = "17.03.2012"
        error_codes = metadata_validation.run_project_checks(self.sub, self.logger)
        self.assertIn("PROJ-E02", error_codes)
        self.assertIn("PROJ-E03", error_codes)
        self.assertIn("PROJ-E04", error_codes)
        self.assertIn("PROJ-E05", error_codes)
        self.assertIn("PROJ-E06", error_codes)
        self.assertIn("PROJ-E07", error_codes)
        self.assertIn("PROJ-E09", error_codes)

        # Deleting all contacts, should give error PROJ-E01
        self.sub.project.contacts = []
        # Deleting release date should give error PROJ-E08
        self.sub.project.releaseDate = ""
        error_codes = metadata_validation.run_project_checks(self.sub, self.logger)
        self.assertIn("PROJ-E01", error_codes)
        self.assertIn("PROJ-E08", error_codes)

    def test_assay_validation(self):
        # All should be fine
        error_codes = metadata_validation.run_assay_checks(self.sub, self.logger)
        self.assertEqual(error_codes, [])

        # Deleting alias from first assay should give ASSA-E02
        self.sub.assay[0].alias = ""
        # Deleting technology type should give ASSA-E03
        self.sub.assay[1].technology_type = ""
        # Changing technology type to unknown type should give ASSA-E04
        self.sub.assay[2].technology_type = "whatever"
        # Changing array design to non AE accession should give ASSA-E08
        self.sub.assay[1].array_design = "this is not in AE"
        error_codes = metadata_validation.run_assay_checks(self.sub, self.logger)
        self.assertIn("ASSA-E02", error_codes)
        self.assertIn("ASSA-E03", error_codes)
        self.assertIn("ASSA-E04", error_codes)
        self.assertIn("ASSA-E08", error_codes)

        # Deleting label should give ASSA-E06
        self.sub.assay[1].label = ""
        # Deleting array design should give ASSA-E07
        self.sub.assay[1].array_design = ""
        error_codes = metadata_validation.run_assay_checks(self.sub, self.logger)
        self.assertIn("ASSA-E06", error_codes)
        self.assertIn("ASSA-E07", error_codes)

        # Deleting all assays should give error ASSA-E01
        self.sub.assay = []
        error_codes = metadata_validation.run_assay_checks(self.sub, self.logger)
        self.assertEqual(["ASSA-E01"], error_codes)


if __name__ == '__main__':
    unittest.main()

