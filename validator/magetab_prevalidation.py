"""This module contains functions that assert the assumptions about the MAGE-TAB files
that the converter module relies on.
SDRF prevalidation parses the SDRF header and checks that certain nodes are present exactly once

"""
import re

from utils.converter_utils import get_name, get_value, get_controlled_vocabulary


def idf_prevalidation(idf_dict, submission_type, logger, atlas=False):
    """Check that all IDF fields and comments are from the allowed list and can be parsed properly.

    In IDF the field name need to be spelled exactly like defined in the MAGE-TAB spec
    (no variablilty with spaces and capitalisation) allowed for AE loading"""

    # Check if all fields can be parsed
    idf_fields = get_controlled_vocabulary("idf", "magetab")
    # Remove spaces (we don't care about them in order to parse correctly)
    known_fields = [get_name(field) for field in idf_fields]
    # Add controlled comment fields (values in square brackets)
    known_fields.extend(get_controlled_vocabulary("idf_comment_terms"))
    for idf_key in idf_dict:
        if idf_key not in known_fields:
            logger.error("Cannot parse IDF field \"{}\".".format(idf_key))

    # Check fields that can only contain one value
    max1 = ("mage-tabversion", "investigationtitle", "dateofexperiment", "publicreleasedate", "experimentdescription")
    for field in max1:
        if idf_dict.get(field) and len([x for x in idf_dict[field] if x]) > 1:
            logger.error("IDF field \"{}\" contains more than one value. This is not allowed".format(field))

    # Single cell IDF checks
    if atlas and submission_type == "singlecell":
        required_comments = get_controlled_vocabulary("required_singlecell_idf_comments", "atlas")
        print(required_comments)
        for comment in required_comments:
            if comment not in idf_dict:
                logger.error("Comment \"{}\" not found in IDF. Required for Single Cell Atlas.".format(comment))


def sdrf_prevalidation(sdrf_list, header, header_dict, submission_type, logger, atlas=False):
    """Perform basic checks on the SDRF, making sure that all expected nodes and protocols are present,
    and that the basic assumptions about the relationships between samples and extracts are correct.
    """

    # Strip whitespace
    header_names = [get_name(h) for h in header]
    node_positions = []
    samples = []
    extracts = []
    le = []

    # Source Name
    if not present_exactly_once("sourcename", header_names):
        logger.error("Source Name node was not found or more than once.")
    else:
        _add_first_occurance("sourcename", header_dict, node_positions)
        # For later check of sample to extract relationship
        samples = [row[header_dict.get("sourcename")[0]] for row in sdrf_list]

    # Extract Name
    if not present_exactly_once("extractname", header_names):
        logger.error("Extract Name node was not found or more than once.")
    else:
        _add_first_occurance("extractname", header_dict, node_positions)
        extracts = [row[header_dict.get("extractname")[0]] for row in sdrf_list]

    # Labeled Extract Name
    if submission_type == "microarray":
        if not present_exactly_once("labeledextractname", header_names):
            logger.error("Labeled Extract node was not found or more than once.")
        else:
            _add_first_occurance("labeledextractname", header_dict, node_positions)
            le = [row[header_dict.get("labeledextractname")[0]] for row in sdrf_list]

    # Assay Name
    if not present_exactly_once("assayname", header_names):
        logger.error("Assay Name node was not found or more than once.")
    else:
        _add_first_occurance("assayname", header_dict, node_positions)

    # Array Data File
    if not (present_exactly_once("arraydatafile", header_names) or
            present_exactly_once("arraydatamatrixfile", header_names) or
            present_exactly_once("scanname", header_names)):
        logger.error("No raw data node was found or more than one.")

    # Check for duplicated characteristics terms
    characteristics = [get_value(x).lower() for x in header if re.match("characteristics", x, flags=re.IGNORECASE)]
    duplicated_categories = [c for c in characteristics if not present_exactly_once(c, characteristics)]
    if duplicated_categories:
        logger.error("The following characteristics categories are present more than once: {}".format(
            ", ".join(duplicated_categories)))

    # Protocols should be between all major nodes
    ref_positions = header_dict.get("protocolref")
    # Add first node for each data file node type to node positions, where we want a protocol
    _add_first_occurance("arraydatafile", header_dict, node_positions)
    _add_first_occurance("scanname", header_dict, node_positions)
    _add_first_occurance("arraydatamatrixfile", header_dict, node_positions)
    _add_first_occurance("derivedarraydatafile", header_dict, node_positions)
    _add_first_occurance("derivedarraydatamatrixfile", header_dict, node_positions)
    # Go through nodes up until the one before last
    for i, node_pos in enumerate(node_positions[:-1]):
        # Get all protocol positions that fulfill the criteria of being between the current and the next node
        next_node_pos = node_positions[i+1]
        prots_in_range = [prot for prot in ref_positions if node_pos < prot < next_node_pos]
        if not prots_in_range:
            logger.warn("There is no Protocol REF connecting \"{}\" and \"{}\".".format(
                header[node_pos], header[next_node_pos]))

    # There must not be more samples than extracts
    if samples and extracts:
        if len(samples) > len(extracts):
            logger.error("Found more samples ({}) than extracts ({}), please check the relationship.".format(
                len(samples), len(extracts)))

    # There must not be more labeled extracts than samples
    if submission_type == "microarray" and samples and le:
        if len(samples) > len(le):
            logger.error("Found more samples ({}) than labeled extracts ({}), please check the relationship.".format(
                len(samples), len(le)))


def cross_magetab_validation(sdrf_list, header, header_dict, idf_dict, submission_type, logger, atlas=False):

    # Factor values must match

    # REF values have to be defined in the IDF
    sdrf_ref_values = set()
    for header in header_dict:
        if re.search(r"ref$", header):
            for index in header_dict[header]:
                # Skipping indexes that are followed by "Term Source REF"
                if index+1 in header_dict.get("termsourceref", []):
                    continue
                ref_values = {row[index] for row in sdrf_list}
                for ref in ref_values:
                    sdrf_ref_values.add(ref)
                    if header == "protocolref" and ref not in idf_dict["protocolname"]:
                        logger.error("Protocol REF \"{}\" is not defined in the IDF.".format(ref))
                    elif header == "termsourceref" and ref not in idf_dict["termsourcename"]:
                        logger.error("Term Source REF \"{}\" is not defined in the IDF.".format(ref))

    # Warn about defined but unused REFs
    for ref_value in idf_dict.get("protocolname", []):
        if ref_value and ref_value not in sdrf_ref_values:
            logger.warn("Protocol \"{}\" is defined in the IDF but not used.".format(ref_value))


def present_exactly_once(item, term_list):
    return term_list.count(item) == 1


def _add_first_occurance(item, lookupdict, target):
    if item in lookupdict:
        if isinstance(lookupdict[item], list):
            target.append(lookupdict[item][0])
