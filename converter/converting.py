"""Module to convert experiment metadata in MAGE-TAB (IDF/SDRF) format to USI submittable JSON format and
helper functions for converting from common data model to USI JSON objects"""

import codecs
import os
import re

from collections import OrderedDict, defaultdict

from converter.datamodel import Attribute, Project, Study, Protocol, Sample, MicroarrayAssay, SeqAssay, DataFile, \
    AssayData, Submission
from converter.parsing import parse_idf, parse_sdrf
from utils.common_utils import create_logger
from utils.converter_utils import is_accession, get_efo_url, strip_extension, write_json_file


def generate_usi_project_object(project):
    project_object = OrderedDict()

    project_object["alias"] = project.alias
    project_object["title"] = project.title
    project_object["description"] = project.description
    project_object["contacts"] = [attrib2dict(contact) for contact in project.contacts]
    project_object["publications"] = [attrib2dict(pub) for pub in project.publications]
    project_object["releaseDate"] = project.releaseDate

    return project_object


def generate_usi_study_object(study, sub_info):
    study_object = OrderedDict()
    study_object["alias"] = study.alias
    study_object["title"] = study.title
    study_object["description"] = study.description

    study_object["studyType"] = "FunctionalGenomics"

    study_object["projectRef"] = generate_usi_ref_object(study.projectref, sub_info)
    study_object["protocolRefs"] = study.protocolrefs
    study_object["protocolRefs"] = [generate_usi_ref_object(ref, sub_info) for ref in study.protocolrefs]

    study_attributes = defaultdict(list)
    # The factors and design types need "value" entry
    for factor in study.experimental_factor:
        attrib_obj = generate_usi_attribute_entry(factor)
        study_attributes["experimental_factor"].extend(attrib_obj)

    for design in study.experimental_design:
        study_attributes["experimental_design"].extend(generate_usi_attribute_entry(design))

    for et in study.experiment_type:
        study_attributes["experiment_type"].extend(generate_usi_attribute_entry(et))

    # Optional attributes
    if study.date_of_experiment:
        study_attributes["date_of_experiment"] = generate_usi_attribute_entry(study.date_of_experiment)

    study_object["attributes"] = study_attributes

    return study_object


def generate_usi_protocol_object(protocol):

    protocol_object = OrderedDict()

    protocol_object["alias"] = protocol.alias
    protocol_object["description"] = protocol.description

    attribute_names = protocol.get_ae_attributes()
    protocol_object['attributes'] = OrderedDict()
    for a in attribute_names:
        protocol_object["attributes"][a] = generate_usi_attribute_entry(getattr(protocol, a))

    return protocol_object


def generate_usi_sample_object(sample):

    sample_object = OrderedDict()
    sample_object['alias'] = sample.alias
    sample_object['taxon'] = sample.taxon
    sample_object['taxonId'] = sample.taxonId

    if sample.description:
        sample_object['description'] = sample

    attributes = {}
    for a_name, a_attrib in sample.attributes.items():
        attributes[a_name] = generate_usi_attribute_entry(a_attrib)
    sample_object['attributes'] = attributes

    return sample_object


def generate_usi_assay_object(assay, study_info):

    assay_object = OrderedDict()

    if assay.accession:
        assay_object['alias'] = assay.accession
        assay_object['title'] = assay.alias
    else:
        assay_object['alias'] = assay.alias

    attributes = assay.get_attributes()
    assay_object['attributes'] = OrderedDict()
    for category in attributes:
        assay_object['attributes'][category] = generate_usi_attribute_entry(getattr(assay, category))

    assay_object['studyRef'] = generate_usi_ref_object(study_info.get('alias'), study_info)

    sample_refs = assay.sampleref
    # Can be one or more samples linked to one assay, but schema is expecting an array
    if not isinstance(sample_refs, list):
        sample_refs = [sample_refs]
    assay_object['sampleUses'] = [generate_usi_ref_object(ref, study_info) for ref in sample_refs]

    protocol_refs = assay.protocolrefs
    if not isinstance(protocol_refs, list):
        protocol_refs = [protocol_refs]
    assay_object["protcolUses"] = [generate_usi_ref_object(ref, study_info) for ref in protocol_refs]

    return assay_object


def generate_usi_data_object(assay_data, sub_info):
    """

    :param assay_data:
    :type:
    :param sub_info:
    :return:
    """
    ad_object = OrderedDict()

    ad_object["alias"] = assay_data.alias
    ad_object["data_type"] = assay_data.data_type

    ad_object["files"] = []
    files = assay_data.files
    for fo in files:
        file_object = attrib2dict(fo)
        ad_object["files"].append(file_object)

    ad_object["AssayRefs"] = [generate_usi_ref_object(x, sub_info) for x in assay_data.assayrefs]

    return ad_object


def generate_usi_attribute_entry(attribute_info):
    """
    Expected input is a dictionary key and value pair. The value can be a dictionary itself.
    In that case, we'll test if we find units and ontology terms and add those accordingly.
    Schema definition:
     "attributes": {
            "description": "Attributes for describing a submittable.",
            "type": "object",
            "properties": {},
            "patternProperties": {
                "^.*$": {
                    "type": "array",
                    "minItems": 1,
                    "items": {
                        "type": "object",  (this line not in USI template spec!)
                        "properties": {
                            "value": { "type": "string", "minLength": 1 },
                            "units": { "type": "string", "minLength": 1 },
                            "terms": {
                                "type": "array",
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "url": {"type": "string", "format": "uri" }
                                    },
                                    "required": ["url"]
                                }
                            }
                        },
                        "required": [ "value" ]
    """

    # This the attribute entry is always a list even though we only expect one value per attribute
    attribute_object = list()
    if isinstance(attribute_info, dict):
        # Initialise the dictionary (with the value lookup) as the first item in the list
        attribute_object.append(OrderedDict([("value", attribute_info.get("value"))]))
        if attribute_info.get("unit"):
            unit_info = attribute_info.get("unit")
            attribute_object[0]["units"] = unit_info.get("value")
            if unit_info.get("term_accession"):
                # USI model does support >1 term URIs for a value but we can't distinguish between
                # the value term and the unit term. For now will only include an ontology URI for
                # a value term (see below).
                pass
        if attribute_info.get("term_accession"):
            term_acc = attribute_info.get("term_accession")
            attribute_object[0]["terms"] = [{"url": get_efo_url(term_acc)}]
    elif isinstance(attribute_info, Attribute):
        attribute_object.append(OrderedDict([("value", attribute_info.value)]))
        if attribute_info.unit:
            unit_info = attribute_info.unit
            attribute_object[0]["units"] = unit_info.value
        if attribute_info.term_accession:
            attribute_object[0]["terms"] = [{"url": get_efo_url(attribute_info.term_accession)}]
    else:
        attribute_object.append({"value": attribute_info})

    return attribute_object


def attrib2dict(ob):
    """Get all attributes of an object (ob) and turn them into a dictionary.
    the attributes will become the keys of the dict and the object values the dict values."""

    attrib_dict = OrderedDict()
    for a in ob.__dict__:
        if getattr(ob, a):
            attrib_dict[a] = getattr(ob, a)
    return attrib_dict


def generate_usi_ref_object(alias, sub_info, accession=None):
    """
    Schema definition:
    "submittableRef": {
        "type": "object",
        "properties": {
        "alias": { "type": "string", "minLength": 1 },
        "accession": { "type": "string", "minLength": 1 },
        "team": { "type": "string", "minLength": 1 }
        },
        "anyOf": [
            { "required": [ "alias", "team" ] },
            { "required": [ "accession" ] }
    """

    ref_object = dict()

    if accession and is_accession(accession):
        ref_object['accession'] = accession
    elif is_accession(alias):
        ref_object['accession'] = alias
    else:
        ref_object['alias'] = alias
        ref_object['team'] = sub_info.get('team')

    return ref_object


SDRF_FILE_NAME_REGEX = r"^\s*SDRF\s*File"
DATA_DIRECTORY = "unpacked"


def mtab2usi_conversion(idf_file_path):
    """
    Run data transformation from a set of IDF/SDRF files to a set of USI JSON files

    This function ties together the two parts of
    - reading in MAGE-TAB data into common data model via the function "data_objects_from_magetab"
    - writing JSON files from data in common data model via the function "datamodel2json_conversion"

    :param idf_file_path: string, file path to IDF file
    :return: None
    """
    process_name = "mtab2usi_conversion"

    current_dir = os.path.dirname(idf_file_path)
    idf_file_name = os.path.basename(idf_file_path)
    sdrf_file_path = None

    # Create logger
    logger = create_logger(current_dir, process_name, idf_file_name)

    # Figure out the name and location of sdrf files
    with codecs.open(idf_file_path, 'rU', encoding='utf-8') as f:
        # U flag makes it portable across in unix and windows (\n and \r\n are treated the same)
        for line in f:
            if re.search(SDRF_FILE_NAME_REGEX, line):
                sdrf_file_name = line.split("\t")[1].strip()
                if os.path.exists(current_dir + DATA_DIRECTORY):
                    sdrf_file_path = os.path.join(current_dir, DATA_DIRECTORY, sdrf_file_name)
                else:
                    sdrf_file_path = os.path.join(current_dir, sdrf_file_name)

    logger.debug("Found SDRF file: {}".format(sdrf_file_path))

    sub = data_objects_from_magetab(idf_file_path, sdrf_file_path)

    datamodel2json_conversion(sub, current_dir, logger)


def data_objects_from_magetab(idf_file_path, sdrf_file_path):
    """
    Parse IDF/SDRF files and transform metadata to common datamodel

    :param idf_file_path: string, path to IDF file
    :param sdrf_file_path: string, path to SDRF file
    :return: Submission class object
    """

    study_info, protocols = parse_idf(idf_file_path)
    samples, extracts, le, assays, raw_data, processed_data, is_microarray = parse_sdrf(sdrf_file_path)

    # For MAGE-TAB files we don't have USI submission info might need to store these somewhere once we get this
    idf_file_name = os.path.basename(idf_file_path)
    sub_info = {"alias": re.sub("\.idf\.txt$", "", idf_file_name),
                "accession": study_info.get("accession"),
                "team": "my-super-test-team",
                "metadata": idf_file_path}

    # Project
    project_object = Project.from_magetab(study_info)

    # Study
    study_object = Study.from_magetab(study_info)

    # Protocols
    protocol_objects = []
    for p in protocols:
        protocol = Protocol.from_magetab(p)
        protocol_objects.append(protocol)

    # Samples
    sample_objects = []
    for sample in samples.values():
        new_sample = Sample.from_magetab(sample)
        sample_objects.append(new_sample)

    # Assays
    assay_objects = []

    if is_microarray:
        linked_extracts = []
        for le_name, le_attributes in le.items():
            for extract_name, extract_attributes in extracts.items():
                if extract_name in le_attributes["extract_ref"]:
                    # Assuming there is only one extract per le
                    linked_extracts = extract_attributes
                    break

            # Get all assays referencing this extract
            linked_assays = []
            for assay_name, assay_attributes in assays.items():
                if le_name in assay_attributes["extract_ref"]:
                    linked_assays.append(assay_attributes)

            new_assay = MicroarrayAssay.from_magetab(le_attributes, linked_extracts, linked_assays)
            assay_objects.append(new_assay)
    # Sequencing assays
    else:
        for extract_name, extract_attributes in extracts.items():
            # Get all assays referencing this extract
            linked_assays = []
            for assay_name, assay_attributes in assays.items():
                if extract_name in assay_attributes["extract_ref"]:
                    linked_assays.append(assay_attributes)

            new_assay = SeqAssay.from_magetab(extract_attributes, linked_assays, protocols)
            assay_objects.append(new_assay)

    # Assay data
    ad_objects = []

    file_groups = OrderedDict()
    for f_name, f_attrib in raw_data.items():
        if len(f_attrib.get("assay_ref")) > 1:
            # For matrix files, one object per file
            file_groups[strip_extension(f_name)] = [f_attrib]

        elif len(f_attrib.get("assay_ref")) == 1:
            a_ref = f_attrib.get("assay_ref")[0]
            # Get the other files with the same assay ref
            file_groups[a_ref] = [f_attrib for f_attrib in raw_data.values() if a_ref in f_attrib.get("assay_ref")]

    for name, group in file_groups.items():
        file_objects = [DataFile.from_magetab(f_attrib) for f_attrib in group]
        assay_data = AssayData.from_magetab(name, file_objects, group)
        ad_objects.append(assay_data)

    # Analysis (processed data)
    print(processed_data)
    analysis_objects = []



    sub = Submission(sub_info,
                     project_object,
                     study_object,
                     protocol_objects,
                     sample_objects,
                     assay_objects,
                     ad_objects,
                     analysis_objects)

    return sub


def datamodel2json_conversion(submission, working_dir, logger):
    """
    Take metadata in common datamodel and write JSON files

    :param submission: object, Submission class object that holds metadata of the whole experiment
    :param working_dir: string, directory to write files to
    :param logger: object, log handler
    :return: None
    """

    # Dict to store USI objects to write to file:
    json_objects = {
        "project": generate_usi_project_object(submission.project),
        "study": generate_usi_study_object(submission.study, submission.info),
        "protocol": [generate_usi_protocol_object(p) for p in submission.protocol],
        "sample": [generate_usi_sample_object(s) for s in submission.sample],
        "assay": [generate_usi_assay_object(a, submission.info) for a in submission.assay],
        "assay_data": [generate_usi_data_object(ad, submission.info) for ad in submission.assay_data]
    }

    # Write individual JSON files
    for submittable_type, objects in json_objects.items():
        logger.info("Writing JSON file for {}.".format(submittable_type))
        write_json_file(working_dir, objects, submittable_type, submission.info)