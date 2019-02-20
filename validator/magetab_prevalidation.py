"""This module contains functions that assert the assumptions about the MAGE-TAB files
that the converter module relies on.
SDRF prevalidation parses the SDRF header and checks that certain nodes are present exactly once

"""
import re

from converter.parsing import get_name, get_value


def sdrf_prevalidation(sdrf_list, header, header_dict, is_microarray, logger):

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
    if is_microarray:
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
        logger.error("The following characteristics categories is present more than once: {}".format(
            ", ".join(duplicated_categories)))

    # Protocols must be between all major nodes
    ref_positions = header_dict.get("protocolref")
    # Add first node for each data file node type to node positions, where we want a protocol
    _add_first_occurance("arraydatafile", header_dict, node_positions)
    _add_first_occurance("scanname", header_dict, node_positions)
    _add_first_occurance("arraydatamatrixfile", header_dict, node_positions)
    _add_first_occurance("derivedarraydatafile", header_dict, node_positions)
    _add_first_occurance("derivedarraydatamatrixfile", header_dict, node_positions)
    print(node_positions)
    # Go through nodes up until the one before last
    for i, node_pos in enumerate(node_positions[:-1]):
        # Get all protocol positions that fulfill the criteria of being between the current and the next node
        next_node_pos = node_positions[i+1]
        prots_in_range = [prot for prot in ref_positions if node_pos < prot < next_node_pos]
        print(prots_in_range)
        if not prots_in_range:
            logger.error("There is no Protocol REF connecting \"{}\" and \"{}\".".format(
                header[node_pos], header[next_node_pos]))

    # There must not be more samples than extracts
    if samples and extracts:
        if len(samples) > len(extracts):
            logger.error("Found more samples ({}) than extracts ({}), please check the relationship.".format(
                len(samples), len(extracts)))

    # There must not be more labeled extracts than samples
    if is_microarray and samples and le:
        if len(samples) > len(le):
            logger.error("Found more samples ({}) than labeled extracts ({}), please check the relationship.".format(
                len(samples), len(le)))


def present_exactly_once(item, term_list):
    return term_list.count(item) == 1


def _add_first_occurance(item, lookupdict, target):
    if item in lookupdict:
        if isinstance(lookupdict[item], list):
            target.append(lookupdict[item][0])
