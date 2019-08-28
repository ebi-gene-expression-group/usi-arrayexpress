"""Module with validation checks to test for eligibility for Expression Atlas and Single Cell Expression Atlas"""
import logging
import re
from collections import defaultdict

from utils.common_utils import is_valid_url
from utils.converter_utils import simple_idf_parser, get_controlled_vocabulary, \
                                  get_name, get_value, read_sdrf_file, get_taxon



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

        # Read in IDF/SDRF files
        try:
            self.sdrf, self.sdrf_header, self.header_dict = read_sdrf_file(sdrf_file)
            self.idf = simple_idf_parser(idf_file)
        except Exception as e:
            raise Exception("Failed to open MAGE-TAB files: {}".format(e))

        # Get metadata fields
        self.sdrf_comment_values = [self.normalise_header(self.sdrf_header[i])
                                    for i in self.header_dict.get('comment', [])]
        self.sdrf_charactistics_values = [self.normalise_header(self.sdrf_header[i])
                                          for i in self.header_dict.get('characteristics', [])]
        self.sdrf_factor_values = [self.normalise_header(self.sdrf_header[i])
                                   for i in self.header_dict.get('factorvalue', [])]
        self.sdrf_values = self.sdrf_comment_values + self.sdrf_charactistics_values + self.sdrf_factor_values
        self.idf_values = [self.normalise_header(field) for field in self.idf]

        self.sample2datafile = self.map_samples_to_files()

    def run_general_checks(self, logger):

        # Warn about technical replicates
        if "technical replicate group" not in self.sdrf_comment_values:
            for s in self.sample2datafile:
                if s.endswith("PAIRED") and len(self.sample2datafile[s]) > 2:
                    logger.warn("Experiment is likely to contain technical replicates"
                                "Please add \"Comment[technical replicate group]\".")
                    break
                elif len(self.sample2datafile[s]) > 1:
                        logger.warn("Experiment is likely to contain technical replicates. "
                                    "Please add \"Comment[technical replicate group]\".")
                        break

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

        # Species checks
        # Only 1 species allowed
        for i, field in enumerate(self.sdrf_header):
            if re.search("organism", get_value(field), re.IGNORECASE):
                organisms = {row[i] for row in self.sdrf}
                if len(organisms) > 1:
                    logger.error("Experiment contains more than 1 organisms.")
                # Organism must be found in NCBI taxonomy
                for o in organisms:
                    taxon_id = get_taxon(o)
                    if not taxon_id:
                        logger.error("Organism \"{}\" is not found in NCBI taxonomy.".format(o))
                break

        # Sequencing experiments must have RUN or ENA_RUN comment column
        if self.submission_type in ("sequencing", "singlecell"):
            if not ("run" in self.sdrf_values or "ena_run" in self.sdrf_values):
                logger.error("No ENA_RUN or RUN column found in SDRF.")

            # Paired-end libraries need two files

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

        # Atlas IDF comment value checks
        for k, attribs in self.idf.items():
            if re.search("EAAdditionalAttributes", k, flags=re.IGNORECASE):
                for attrib in attribs:
                    if attrib and attrib not in self.sdrf_values:
                        logger.error("Additional attribute \"{}\" not found in SDRF.".format(attrib))
            elif re.search("EAExpectedClusters", k, flags=re.IGNORECASE):
                for number in attribs:
                    if number and not re.match("^\d+$", number):
                        logger.error("Expected clusters value \"{}\" is not numerical.".format(number))

        # Required SDRF fields
        required_sdrf_names = get_controlled_vocabulary("required_singlecell_sdrf_fields", "atlas")
        for field in required_sdrf_names:
            if not (field.lower() in self.sdrf_values or get_name(field) in self.header_dict):
                logger.error("Required SDRF field \"{}\" not found.".format(field))

        # Valid SDRF values
        library_construction_terms = get_controlled_vocabulary("singlecell_library_construction", "atlas")
        sc_protocol_values = {}
        for i, c in enumerate(self.sdrf_header):
            # Check for supported library construction terms
            if re.search(r"library construction", self.normalise_header(c), flags=re.IGNORECASE):
                sc_protocol_values = {row[i] for row in self.sdrf}
                print(sc_protocol_values)
                if len(sc_protocol_values) > 1:
                    logger.warn("Experiment contains more than 1 single cell library construction protocol.")
                for protocol in sc_protocol_values:
                    if protocol.lower() not in library_construction_terms.get("all", []):
                        logger.error("Library construction protocol is not supported for Expression Atlas."
                                     .format(protocol))
            # Not all rows should be "not OK"
            elif re.search(r"single cell (well)? quality", self.normalise_header(c), flags=re.IGNORECASE):
                well_quality_values = {row[i] for row in self.sdrf}
                print(well_quality_values)
                if len(well_quality_values) == 1 and "not OK" in well_quality_values:
                    logger.error("Single cell quality values are all \"not OK\".")
            # Technical replicate group values must only contain letters and numbers
            elif re.search(r"technical replicate group", self.normalise_header(c), flags=re.IGNORECASE):
                tech_rep_values = {row[i] for row in self.sdrf}
                wrong_values = {v for v in tech_rep_values if not re.match(r"^[A-Za-z0-9]*$", v)}
                if len(wrong_values) > 0:
                    logger.error("Technical replicate group values can only contain letters and numbers: "
                                 .format(", ".join([str(v) for v in wrong_values])))
                # Technical replicate group column should not be empty
                if len(tech_rep_values) == 1 and "" in tech_rep_values:
                    logger.error("Technical replicate group values are all empty.")

        # SDRF terms required for droplet experiments
        for protocol in sc_protocol_values:
            if protocol.lower() in library_construction_terms.get("droplet"):
                droplet_terms = get_controlled_vocabulary("required_droplet_sdrf_fields", "atlas")
                for dt in droplet_terms:
                    if dt not in self.sdrf_values:
                        logger.error("Required SDRF droplet field \"{}\" not found.".format(dt))
                break

    def check_all(self, logger=logging.getLogger()):
        """Trigger all applicable checks"""

        self.run_general_checks(logger)
        if self.submission_type == "singlecell":
            self.run_singlecell_checks(logger)

    @staticmethod
    def normalise_header(field_name):
        """Strip field names such as Comment and make everything lowercase without spaces"""
        return get_value(field_name).lower().strip()

    def map_samples_to_files(self):
        """Make a dictionary with the mapping of samples to their associated raw data files"""
        # Compare source name values against files
        sample2datafile = defaultdict(list)
        # Find source name and data file index
        sample_node = "sourcename"
        if self.submission_type == "microarray":
            data_node = "arraydatafile"
        else:
            data_node = "scanname"
        layout_comment = "library_layout"
        # The fist columns for each node (expecting only one raw data column)
        sample_index = next(iter(self.header_dict.get(sample_node, [])), None)
        data_index = next(iter(self.header_dict.get(data_node, [])), None)
        layout_index = next(iter([i for i, x in enumerate(self.sdrf_header)
                                  if re.search(layout_comment, x, flags=re.IGNORECASE)]), None)
        # Collect all data files for each sample
        if sample_index is not None and data_index is not None:
            for row in self.sdrf:
                if layout_index is not None:
                    name = str(row[sample_index])+"~"+str(row[layout_index])
                    sample2datafile[name].append(row[data_index])
                else:
                    sample2datafile[row[sample_index]].append(row[data_index])

        return sample2datafile
