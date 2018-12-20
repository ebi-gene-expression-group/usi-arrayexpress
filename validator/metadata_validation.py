
from converter import datamodel


# Input is a Submission object generated from MAGE-TAB or JSON


def run_protocol_checks(sub: datamodel.Submission, logger):
    """Run checks on protocol objects and return list of error codes."""

    protocols = sub.protocol

    codes = []
    names = set()

    if not protocols:
        logger.error("Experiment has no protocols. At least one expected.")
        codes.append("PROT-E01")
        return codes
    for p in protocols:
        if p.alias:
            # Protocol names should be unique.
            if p.alias in names:
                logger.error("Protocol name {} is not unique.".format(p.alias))
                codes.append("PROT-E04")
            names.add(p.alias)
        # Protocol must have a name
        else:
            logger.error("Protocol found with no name. Not checking it further.")
            codes.append("PROT-E02")
            continue
        if p.description:
            # Protocol description should be longer than 50 characters.
            if len(p.description) < 50:
                logger.warning("Protocol {} is shorter than 50 characters.".format(p.alias))
                codes.append("PROT-W01")
        # Protocol must have description
        else:
            logger.error("Protocol {} has no description.".format(p.alias))
            codes.append("PROT-E03")

        if p.protocol_type:
            # Protocol type must be from controlled vocabulary (EFO)
            pass
        else:
            logger.warn("Protocol {} has no protocol type.".format(p.alias))
            codes.append("PROT-W02")

    return codes










