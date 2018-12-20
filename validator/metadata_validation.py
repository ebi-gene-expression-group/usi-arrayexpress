
from converter import datamodel
from utils.converter_utils import ontology_lookup

# Input is a Submission object generated from MAGE-TAB or JSON


def run_protocol_checks(sub: datamodel.Submission, logger):
    """Run checks on protocol objects and return list of error codes."""

    protocols = sub.protocol

    codes = []
    names = set()
    p_types = set()
    allowed_types = ontology_lookup("protocol_types")
    mandatory = [label for label, attrib in allowed_types.items()
                 if attrib["exp_type"] == "all" and
                 (attrib["mandatory"] == "ma" or attrib["mandatory"] == "seq")]
    exclusive = [label for label, attrib in allowed_types.items()
                 if attrib["exp_type"] == "all" and
                 attrib["mandatory"] == "one of"]
    found_exclusive = False

    if not protocols:
        logger.error("Experiment has no protocols. At least one expected.")
        codes.append("PROT-E01")
        return codes
    for p in protocols:
        if p.alias:
            # Protocol names should be unique.
            if p.alias in names:
                logger.error("Protocol name \"{}\" is not unique.".format(p.alias))
                codes.append("PROT-E04")
            names.add(p.alias)
        # Protocol must have a name
        else:
            logger.error("Protocol found with no name. Not checking it further.")
            codes.append("PROT-E02")
            continue
        if p.description:
            # Protocol description should be longer than 50 characters
            if len(p.description) < 50:
                logger.warning("Protocol \"{}\" is shorter than 50 characters.".format(p.alias))
                codes.append("PROT-W01")
        # Protocol must have description
        else:
            logger.error("Protocol \"{}\" has no description.".format(p.alias))
            codes.append("PROT-E03")
        if p.protocol_type:
            # Protocol type must be from controlled vocabulary (EFO)
            p_types.add(p.protocol_type.value)
            if p.protocol_type.value not in allowed_types:
                logger.error("Protocol \"{}\" has a type that is not from controlled vocabulary/EFO: "
                             "\"{}\"".format(p.alias, p.protocol_type))
                codes.append("PROT-E05")
            if p.protocol_type.value in exclusive:
                found_exclusive = True
        else:
            # Protocol must have a protocol type
            logger.warn("Protocol \"{}\" has no protocol type.".format(p.alias))
            codes.append("PROT-W02")

    # Mandatory protocol types (for all experiment types) must be present
    for p_type in mandatory:
        if p_type not in p_types:
            logger.error("A {} must be included.".format(p_type))
            codes.append("PROT-E06")

    # Every experiment must have at least one growth/treatment/sample collection protocol
    if not found_exclusive:
        logger.error("A growth, treatment or sample collection protocol must be included.")
        codes.append("PROT-E07")

    return codes













