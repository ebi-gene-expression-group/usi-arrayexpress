"""Module with validation checks to test for eligibility for Expression Atlas and Single Cell Expression Atlas"""
import logging
import re

import pandas

from utils.converter_utils import simple_idf_parser, get_controlled_vocabulary, get_name, get_value, read_sdrf_file
from utils.common_utils import create_logger, is_valid_url


def run_singlecell_atlas_checks(idf, sdrf, header_dict, logger):
    """Check requirements for loading an experiment into Single Cell Expression Atlas"""

    # Single cell IDF checks
    required_comments = get_controlled_vocabulary("required_singlecell_idf_comments", "atlas")
    for comment in required_comments:
        if get_value(comment.lower()) not in idf:
            logger.error("Comment \"{}\" not found in IDF. Required for Single Cell Atlas.".format(comment))

    # Required SDRF fields
    required_sdrf_names = get_controlled_vocabulary("required_singlecell_sdrf_fields", "atlas")
    for field in required_sdrf_names:
        if get_value(field.lower()) not in header_dict:
            logger.error("Required SDRF field \"{}\" not found.".format(field))

    # Check that library_construction is "Smart-seq2", "Drop-Seq


def run_general_atlas_checks(idf, sdrf, sdrf_header, header_dict, submission_type, logger):

    # Extract values from SDRF headers in square brackets e.g. [age]
    sdrf_comment_values = [get_value(sdrf_header[i]).lower() for i in header_dict.get('comment', [])]
    sdrf_charactistics_values = [get_value(sdrf_header[i]).lower()
                                 for i in header_dict.get('characteristics', [])]
    sdrf_factor_values = [get_value(sdrf_header[i]).lower() for i in header_dict.get('factorvalue', [])]
    sdrf_values = sdrf_comment_values + sdrf_charactistics_values + sdrf_factor_values

    # Required IDF fields
    required_idf_fields = get_controlled_vocabulary("required_idf_fields", "atlas")
    for field in required_idf_fields:
        if get_name(field.lower()) not in idf:
            logger.error("No {} found in IDF.".format(field))

    # No duplications between comments and characteristics
    duplicates = set(sdrf_comment_values).intersection(sdrf_charactistics_values)
    duplicates.update(set(sdrf_comment_values).intersection(sdrf_factor_values))
    for c in duplicates:
        logger.error("Column name \"{}\" appears as Comment and Characteristics/Factor Value.".format(c))

    # Sequencing experiments must have RUN or ENA_RUN comment column
    if submission_type in ("sequencing", "singlecell"):
        if not ("run" in sdrf_comment_values or "ena_run" in sdrf_comment_values):
            logger.error("No ENA_RUN or RUN column found in SDRF.")

    # FASTQ_URIs must be valid
    uri_index = [i for i, c in enumerate(sdrf_header) if re.search("fastq_uri", c, flags=re.IGNORECASE)]
    logger.info("Checking FASTQ URIs. This may take a while... (Skip this check with -x option)")
    for row in sdrf:
        for i in uri_index:
            url = row[i]
            # only run the check on internet URLs, for internal files we use just the filename
            if re.match(re.compile("^ftp|^http", re.IGNORECASE), str(url)):
                if not is_valid_url(url):
                    logger.error("FASTQ_URI {} is not valid.".format(url))


    # Warn about technical replicates
    # sdrf_fields = map(get_name, self.sdrf.columns)
    # if 'sourcename' in sdrf_fields:
    #    sample_number = self.sdrf['Source Name'].value_count()
    #    if 'arraydatafile' in sdrf_fields:
    #        pass
    #    elif 'scanname' in sdrf_fields:
    #        pass

    # if files per sample >1 for se

    # if files_per_sample >2 for pe

    # if files_per_sample >3 for 10x


def run_all_atlas_checks(idf, sdrf, sdrf_header, header_dict, submission_type, logger):
    """Trigger all applicable checks"""

    run_general_atlas_checks(idf, sdrf, sdrf_header, header_dict, submission_type, logger)

    if submission_type == "singlecell":
        run_singlecell_atlas_checks(idf, sdrf, header_dict, logger)






class AtlasMAGETABChecker:
    def __init__(self, idf_file, sdrf_file, submission_type, skip_file_checks=False):
        self.idf_file = idf_file
        self.sdrf_file = sdrf_file
        self.submission_type = submission_type
        self.skip_file_checks = skip_file_checks

        try:
            self.sdrf, self.sdrf_header, self.header_dict = read_sdrf_file(sdrf_file)
            self.idf = simple_idf_parser(idf_file)
        except Exception as e:
            raise Exception("Failed to open MAGE-TAB files: {}".format(e))

        self.sdrf_comment_values = [self.normalise_header(self.sdrf_header[i])
                                    for i in self.header_dict.get('comment', [])]
        self.sdrf_charactistics_values = [self.normalise_header(self.sdrf_header[i])
                                          for i in self.header_dict.get('characteristics', [])]
        self.sdrf_factor_values = [self.normalise_header(self.sdrf_header[i])
                                   for i in self.header_dict.get('factorvalue', [])]
        self.sdrf_values = self.sdrf_comment_values + self.sdrf_charactistics_values + self.sdrf_factor_values
        self.idf_values = [self.normalise_header(field) for field in self.idf]
        print(self.idf_values)

    def run_general_checks(self, logger):

        # Warn about technical replicates
        #sdrf_fields = map(get_name, self.sdrf.columns)
        #if 'sourcename' in sdrf_fields:
        #    sample_number = self.sdrf['Source Name'].value_count()
        #    if 'arraydatafile' in sdrf_fields:
        #        pass
        #    elif 'scanname' in sdrf_fields:
        #        pass

        #if files per sample >1 for se

        #if files_per_sample >2 for pe

        #if files_per_sample >3 for 10x


        # Required IDF fields
        required_idf_fields = get_controlled_vocabulary("required_idf_fields", "atlas")
        for field in required_idf_fields:
            if field.lower() not in self.idf_values:
                logger.error("No {} found in IDF.".format(field))

        # No duplications between comments and characteristics
        duplicates = set(self.sdrf_comment_values).intersection(self.sdrf_charactistics_values)
        duplicates.update(set(self.sdrf_comment_values).intersection(self.sdrf_factor_values))
        for c in duplicates:
            logger.error("Column name \"{}\" appears as Comment and Characteristics/Factor Value.".format(c))

        # Sequencing experiments must have RUN or ENA_RUN comment column
        if self.submission_type in ("sequencing", "singlecell"):
            if not ("run" in self.sdrf_values or "ena_run" in self.sdrf_values):
                logger.error("No ENA_RUN or RUN column found in SDRF.")

        # FASTQ_URIs must be valid
        if not self.skip_file_checks:
            uri_index = [i for i, c in enumerate(self.sdrf_header)
                         if re.search("fastq_uri", c, flags=re.IGNORECASE)]
            logger.info("Checking FASTQ URIs. This may take a while... (Skip this check with -x option)")
            for row in self.sdrf:
                for i in uri_index:
                    url = row[i]
                    # only run the check on internet URLs, for internal files we use just the filename
                    if re.match(re.compile("^ftp|^http", re.IGNORECASE), str(url)):
                        if not is_valid_url(url):
                            logger.error("FASTQ_URI {} is not valid.".format(url))




    def run_singlecell_checks(self, logger):
        """Check requirements for loading an experiment into Single Cell Expression Atlas"""

        # Single cell IDF checks
        required_comments = get_controlled_vocabulary("required_singlecell_idf_comments", "atlas")
        for comment in required_comments:
            if comment.lower() not in self.idf_values:
                logger.error("Comment \"{}\" not found in IDF. Required for Single Cell Atlas.".format(comment))

        # Required SDRF fields
        required_sdrf_names = get_controlled_vocabulary("required_singlecell_sdrf_fields", "atlas")
        for field in required_sdrf_names:
            if field.lower() not in self.sdrf_values:
                logger.error("Required SDRF field \"{}\" not found.".format(field))

        # Valid library construction terms
        library_construction_terms = get_controlled_vocabulary("singlecell_library_contruction", "atlas")
        sc_protocol_index = None
        for i, c in enumerate(self.sdrf_header):
            if re.search("library construction", self.normalise_header(c), flags=re.IGNORECASE):
                sc_protocol_index = i
                break
        sc_protocol_values = {row[sc_protocol_index] for row in self.sdrf}
        print(sc_protocol_values)
        if len(sc_protocol_values) > 1:
            logger.warn("Experiment contains more than 1 single cell library construction protocol.")
        for protocol in sc_protocol_values:
            if protocol not in library_construction_terms.get("all", []):
                logger.error("Library construction protocol is not supported for Expression Atlas.".format(protocol))

    def check_all(self, logger):
        """Trigger all applicable checks"""

        if not logger:
            logger = logging.getLogger()

        self.run_general_checks(logger)
        if self.submission_type == "singlecell":
            self.run_singlecell_checks(logger)

    @staticmethod
    def normalise_header(field_name):
        """Strip field names such as Comment and make everything lowercase without spaces"""
        return get_value(field_name).lower()


