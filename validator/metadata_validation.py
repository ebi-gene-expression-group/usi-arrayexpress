
from converter import datamodel
from utils.converter_utils import ontology_term
from utils.common_utils import get_term_descendants

# Input is a Submission object generated from MAGE-TAB or JSON


def run_protocol_checks(sub: datamodel.Submission, logger):
    """Run checks on protocol objects and return list of error codes."""

    protocols = sub.protocol

    codes = []
    names = set()
    p_types = set()
    allowed_types = ontology_term("protocol_types")
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
                             "\"{}\"".format(p.alias, p.protocol_type.value))
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


def run_sample_checks(sub: datamodel.Submission, logger):
    """Run checks on sample objects and factor values and return list of error codes."""

    samples = sub.sample
    factors = [f.value for f in sub.study.experimental_factor]
    units = set()
    characteristics = []
    codes = []

    if not samples:
        logger.error("Experiment has no samples. At least one expected.")
        codes.append("SAMP-E01")
        return codes
    for s in samples:
        # Sample must have a name
        if not s.alias:
            logger.error("Sample found with no name. Not checking it further.")
            codes.append("SAMP-E02")
            continue
        # Sample must have organism/taxon annotation
        if not s.taxon:
            logger.error("Sample \"{}\" has no organism specified.".format(s.alias))
            codes.append("SAMP-E03")
        # Collecting units and categories
        for a, a_attr in s.attributes.items():
            if a_attr.unit and a_attr.unit.value not in units:
                units.add(a_attr.unit.value)
            if a not in characteristics:
                characteristics.append(a)

    # Check units
    unit_term = ontology_term("unit")
    allowed_units = get_term_descendants(unit_term["ontology"], unit_term["uri"], logger)
    for unit_label in units:
        if unit_label not in allowed_units:
            logger.error("Unit \"{}\" is not from approved list (EFO term).".format(unit_label))
            codes.append("SAMP-E04")

    # Check that factors defined in study are found in sample attributes
    undefined_factors = [f for f in factors if f not in characteristics]
    if len(undefined_factors) > 0:
        logger.error("The following factors are declared but not annotated: {}".format(", ".join(undefined_factors)))
        codes.append("SAMP-E05")

    return codes










