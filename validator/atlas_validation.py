"""Module with validation checks to test for eligibility for Expression Atlas and Single Cell Expression Atlas"""

import re

import pandas

from utils.converter_utils import simple_idf_parser, get_controlled_vocabulary, get_name, get_value
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
    sdrf_values = sdrf_comment_values + sdrf_charactistics_values

    # Required IDF fields
    required_idf_fields = get_controlled_vocabulary("required_idf_fields", "atlas")
    for field in required_idf_fields:
        if get_name(field.lower()) not in idf:
            logger.error("No {} found in IDF.".format(field))

    # No duplications between comments and characteristics
    duplicates = set(sdrf_comment_values).intersection(sdrf_charactistics_values)
    for c in duplicates:
        logger.error("Column name \"{}\" appears as comment and characteristics.".format(c))

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
    def __init__(self, idf_file, sdrf_file, submission_type):
        self.idf_file = idf_file
        self.sdrf_file = sdrf_file
        self.submission_type = submission_type

        try:
            self.sdrf = pandas.read_csv(sdrf_file, sep='\t', encoding='utf-8', comment='#')
            self.idf = simple_idf_parser(idf_file)
            # Make list of headers with field names removed
            # (includes stuff like "Derived ArrayExpress FTP file].1")
            self.sdrf_values = [get_value(h).lower() for h in self.sdrf.columns]
            self.idf_values = [get_value(h).lower() for h in self.idf]
            for h in self.sdrf_values:
                print(h)
        except Exception as e:
            raise Exception("Failed to open MAGE-TAB files: {}".format(e))


    def run_general_checks(self, logger):

        # Warn about technical replicates
        sdrf_fields = map(get_name, self.sdrf.columns)
        #if 'sourcename' in sdrf_fields:
        #    sample_number = self.sdrf['Source Name'].value_count()
        #    if 'arraydatafile' in sdrf_fields:
        #        pass
        #    elif 'scanname' in sdrf_fields:
        #        pass

        #if files per sample >1 for se

        #if files_per_sample >2 for pe

        #if files_per_sample >3 for 10x

        print(self.idf_values)
        # Required IDF fields
        required_idf_fields = get_controlled_vocabulary("required_idf_fields", "atlas")
        for field in required_idf_fields:
            if field.lower() not in self.idf_values:
                logger.error("No {} found in IDF.".format(field))

        # No duplications between comments and characteristics
        unique = set()
        for i, c in enumerate(self.sdrf_values):
            if c not in unique:
                unique.add(c)
            elif get_name(self.sdrf.columns[i]) != "factorvalue":
                logger.error("Column name \"{}\" appears more than once.".format(c, self.sdrf.columns[i]))

        # Sequencing experiments must have RUN or ENA_RUN comment column
        if self.submission_type in ("sequencing", "singlecell"):
            if not ("run" in self.sdrf_values or "ena_run" in self.sdrf_values):
                logger.error("No ENA_RUN or RUN column found in SDRF.")

        # FASTQ_URIs must be valid
        file_uris = self.sdrf.filter(regex=re.compile("fastq_uri", re.IGNORECASE), axis=1)
        logger.info("Checking FASTQ URIs. This may take a while... (Skip this next time with -x option)")
        for row in file_uris.itertuples():
            for url in row:
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

    def check_all(self, logger):
        """Trigger all applicable checks"""

        if not logger:
            logger = create_logger()

        self.run_general_checks(logger)
        if self.submission_type == "singlecell":
            self.run_singlecell_checks(logger)

    def _normalise_header(self, field_name):
        """Strip field names such as Comment and make everything lowercase without spaces"""
        return get_value(field_name).lower()


