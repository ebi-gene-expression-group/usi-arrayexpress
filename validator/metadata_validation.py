import re

from converter import datamodel
from utils import converter_utils
from utils.converter_utils import ontology_term, get_controlled_vocabulary, is_accession
from utils.common_utils import get_term_descendants, get_ena_library_terms_via_usi


REGEX_DATE_FORMAT = re.compile("([12]\d{3}-(0[1-9]|1[0-2])-(0[1-9]|[12]\d|3[01]))")
REGEX_DOI_FORMAT = re.compile("^10\.\d{4,9}\/\S+$")


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
            codes.append("PROT-E07")

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
    organisms = set()
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
        else:
            organisms.add(s.taxon)
        # Collecting units and categories
        for a, a_attr in s.attributes.items():
            if a_attr.unit and a_attr.unit.value:
                if a_attr.unit.value not in units:
                    units.add(a_attr.unit.value)
            if a not in characteristics:
                characteristics.append(a)

    # Check organism name is in taxonomy
    for o in organisms:
        taxon_id = converter_utils.get_taxon(o)
        logger.debug("Found taxon ID: ", taxon_id)
        if not isinstance(taxon_id, int):
            logger.error("Organism \"{}\" was not found in NCBI taxonomy.".format(o))
            codes.append("SAMP-E08")

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

    # Check that factor values vary
    unique_factor_values = {}
    for f in factors:
        if f not in undefined_factors:
            # Get all values for a given factor
            factor_values = [s.attributes[f].value.rstrip() for s in samples]
            # Filter duplicated values to get the number of unique entries
            unique_factor_values[f.lower()] = converter_utils.remove_duplicates(factor_values)
    non_factors = [f_name for f_name, values in unique_factor_values.items() if len(values) < 2]
    good_factors = [f_name for f_name, values in unique_factor_values.items() if len(values) > 1]
    # Go through factors and check for special cases that are exempt from the rule
    for f in non_factors:
        # Special case dose
        if f == "dose":
            # Error if both dose + compound/irradiade do not vary
            if "compound" in non_factors or "irradiate" in non_factors:
                logger.error("For factor values including dose, at least one must vary.")
                codes.append("SAMP-E07")
            # Allow dose to not vary if compound/irradiate do
            elif "compound" in good_factors or "irradiate" in good_factors:
                continue
        # For compound/irradiate, we already check this above. Need to skip here to suppress second error message
        elif f == "compound" and "dose" in non_factors:
            continue
        elif f == "irradiate" and "dose" in non_factors:
            continue
        # Special case immunoprecipitate
        elif f == "immunoprecipitate":
            logger.info("Found factor \"immunoprecipitate\". This doesn't need to vary.")
        else:
            logger.error("Factor value \"{}\" does not vary.".format(f))
            codes.append("SAMP-E06")

    return codes


def run_study_checks(sub: datamodel.Submission, logger):
    """Run checks on study object and return list of error codes."""

    study = sub.study
    codes = []

    # Title
    if not study.title:
        logger.error("Study does not have a title.")
        codes.append("STUD-E01")
    elif len(study.title) > 255:
        logger.warn("Study title may be too long. Max 255 characters are allowed.")
        codes.append("STUD-W01")

    # Description
    if not study.description:
        logger.error("Study does not have a description. Experiment description must be specified.")
        codes.append("STUD-E02")

    # Experiment types
    if not study.experiment_type:
        logger.error("Study does not have any experiment types.")
        codes.append("STUD-E03")
    else:
        allowed_types = ontology_term("experiment_type")
        all_types = allowed_types["microarray"] + allowed_types["sequencing"] + allowed_types["singlecell"]
        for et in study.experiment_type:
            if et not in all_types:
                logger.error("Experiment type \"{}\" is not an allowed experiment type.".format(et))
                codes.append("STUD-E04")

    # Protocol refs
    if len(study.protocolrefs) < 1:
        logger.error("At least one protocol must be used in an experiment")
        codes.append("STUD-E07")

    # Experimental factors
    if not study.experimental_factor:
        logger.error("Study does not have any experimental variables. At least one must be included.")
        codes.append("STUD-E08")

    # Experimental design
    if study.experimental_design:
        design_term = ontology_term("study_design")
        allowed_designs = get_term_descendants(design_term["ontology"], design_term["uri"], logger)
        for dt in study.experimental_design:
            if dt.value not in allowed_designs:
                logger.error("Experimental design \"{}\" is not an allowed term.".format(dt.value))
                codes.append("STUD-E09")

    # Date format
    if study.date_of_experiment:
        if not REGEX_DATE_FORMAT.match(study.date_of_experiment):
            logger.error("Date of experiment must be in YYYY-MM-DD format.")
            codes.append("STUD-E05")

    # Accession format of related experiments
    rel_exp_label = "RelatedExperiment"
    if rel_exp_label in study.comments:
        for acc in study.comments[rel_exp_label]:
            if not is_accession(acc):
                logger.error("Related experiment \"{}\" does not match allowed accession pattern.".format(acc))
                codes.append("STUD-E06")

    return codes


def run_project_checks(sub: datamodel.Submission, logger):
    """Run checks on project object and return list of error codes."""

    project = sub.project
    codes = []
    found_submitter = False
    found_submitter_details = False

    # Contacts
    if not project.contacts:
        logger.error("No contacts found. At least one contact must be included.")
        codes.append("PROJ-E01")
    else:
        # Roles
        role_term = ontology_term("role")
        allowed_roles = get_term_descendants(role_term["ontology"], role_term["uri"], logger)
        for i, c in enumerate(project.contacts):
            if c.roles:
                for r in c.roles:
                    role_value = r.lower().rstrip()
                    if role_value not in allowed_roles:
                        logger.warning("Contact role \"{}\" is not an allowed term.".format(role_value))
                        codes.append("PROJ-E05")
                    elif role_value == "submitter":
                        found_submitter = True
                        if c.email and c.affiliation:
                            found_submitter_details = True
            if not c.lastName:
                logger.error("A contact must have last name specified: {}.".format(c))
                codes.append("PROJ-E02")
        # At least one contact must have role "submitter"
        if not found_submitter:
            logger.error("At least one contact must have role \"submitter\".")
            codes.append("PROJ-E03")
        # At least one submitter contact needs email and affiliation
        if not found_submitter_details:
            logger.error("At least one contact with role \"submitter\" must have email and affiliation specified.")
            codes.append("PROJ-E04")

    # Format of PubMed ID and DOI
    if project.publications:
        for pub in project.publications:
            if pub.pubmedId:
                try:
                    int(pub.pubmedId)
                except ValueError:
                    logger.error("PubMed ID must be numerical. Got \"{}\".".format(pub.pubmedId))
                    codes.append("PROJ-E06")
            if pub.doi:
                if not REGEX_DOI_FORMAT.match(pub.doi.rstrip()):
                    logger.error("Publication DOI \"{}\" does not match expected pattern.".format(pub.doi))
                    codes.append("PROJ-E07")

    # Release date
    if project.releaseDate:
        if not REGEX_DATE_FORMAT.match(project.releaseDate):
            logger.error("Release date \"{}\" is not in YYYY-MM-DD format.".format(project.releaseDate))
            codes.append("PROJ-E09")
    else:
        logger.error("No release date found. Project must have release date specified.")
        codes.append("PROJ-E08")

    return codes


def run_assay_checks(sub: datamodel.Submission, logger):

    assays = sub.assay
    exptype = sub.info["submission_type"]
    codes = []

    # Get controlled terms from ENA's assay schema
    if exptype == "sequencing":
        library_terms = get_ena_library_terms_via_usi(logger)

    if not assays:
        logger.error("Experiment has no assays. At least one expected.")
        codes.append("ASSA-E01")
        return codes

    for a in assays:
        additional_attributes = a.get_attributes()
        # Assay must have name
        if not a.alias:
            logger.error("Assay \"{}\" does not have a name specified. Not checking it.".format(a))
            codes.append("ASSA-E02")
            continue
        # Technology type
        if not a.technology_type:
            logger.error("Assay \"{}\" does not have technology type specified.".format(a.alias))
            codes.append("ASSA-E03")
        elif exptype == "microarray" and a.technology_type != "array assay":
            logger.error("Technology Type must be 'array assay' in a microarray submission. "
                         "Found \"{}\".".format(a.technology_type))
            codes.append("ASSA-E04")
        elif exptype == "sequencing" and a.technology_type != "sequencing assay":
            logger.error("Technology Type must be 'sequencing assay' in a sequencing submission. "
                         "Found \"{}\".".format(a.technology_type))
            codes.append("ASSA-E05")

        # Microarray checks for label and array design
        if exptype == "microarray":
            # Label
            if not a.label:
                logger.error("Microarray assay \"{}\" does not have label attribute specified.".format(a.alias))
                codes.append("ASSA-E06")
            # Array design
            if not a.array_design:
                logger.error("Microarray assay \"{}\" does not have array design specified.".format(a.alias))
                codes.append("ASSA-E07")
            elif not is_accession(a.array_design, "ARRAYEXPRESS"):
                logger.error("Array design \"{}\" is not a valid ArrayExpress accession number.".format(a.array_design))
                codes.append("ASSA-E08")

        # Sequencing checks
        elif exptype == "sequencing":
            # Absense of MA fields
            if "label" in additional_attributes:
                logger.error("Found sequencing assay \"{}\" with 'label' attribute.".format(a.alias))
                codes.append("ASSA-E09")
            if "array_design" in additional_attributes:
                logger.error("Found sequencing assay \"{}\" with 'array design' attribute.".format(a.alias))
                codes.append("ASSA-E10")
            # ENA library terms must match against controlled vocabulary
            for term, cv in library_terms.items():
                # Assuming here that the names of the fields are exactly the same as in the assay attributes
                value = getattr(a, term)
                if value not in cv:
                    logger.error("Value \"{}\" for {} does not match against ENA's controlled vocabulary".format(value, term))
                    codes.append("ASSA-E11")

    return codes

