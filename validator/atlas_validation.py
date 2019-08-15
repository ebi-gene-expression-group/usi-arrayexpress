"""Module with validation checks to test for eligibility for Expression Atlas and Single Cell Expression Atlas"""


import pandas

from converter.datamodel.submission import Submission
from utils.converter_utils import simple_idf_parser


def run_singlecell_atlas_checks(sub: Submission):
    pass

    # Check that library_construction is "Smart-seq2", "Drop-Seq




class AtlasMAGETABChecker:
    def __init__(self, idf_file, sdrf_file, submission_type):
        self.idf = idf_file
        self.sdrf = sdrf_file
        self.submission_type = submission_type

        try:
            self.sdrf_dict = pandas.read_csv(sdrf_file, sep='\t', encoding='utf-8', comment='#')
            self.idf_dict = simple_idf_parser(idf_file)
        except Exception as e:
            raise Exception("Failed to open MAGE-TAB files: {}".format(e))

    def run_singlecell_checks(self):
        pass