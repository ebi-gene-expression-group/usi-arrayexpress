"""Helper functions for converting from MAGE-TAB to USI JSON and back"""

import codecs
import json
from collections import OrderedDict, defaultdict
from utils.converter_utils import is_accession
from converter.datamodel import Attribute


def generate_usi_project_object(project):
    project_object = OrderedDict()

    project_object["alias"] = project.alias
    project_object["title"] = project.title
    project_object["description"] = project.description
    project_object["contacts"] = project.contacts
    project_object["publications"] = project.publications
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
        file_object = OrderedDict()
        for a in fo.__dict__:
            if getattr(fo, a):
                file_object[a] = getattr(fo, a)
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
            # TODO: Function that looks up term accessions and returns EFO/OLS URI
            # Using term accession for now
            attribute_object[0]["terms"] = [{"url": attribute_info.get("term_accession")}]
    elif isinstance(attribute_info, Attribute):
        attribute_object.append(OrderedDict([("value", attribute_info.value)]))
        if attribute_info.unit:
            unit_info = attribute_info.unit
            attribute_object[0]["units"] = unit_info.value
        if attribute_info.term_accession:
            attribute_object[0]["terms"] = [{"url": attribute_info.term_accession}]
    else:
        attribute_object.append({"value": attribute_info})

    return attribute_object


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


def usi_object_file_name(object_type, study_info):

    if study_info.get('accession'):
        return "{}_{}.json".format(study_info.get('accession'), object_type)
    elif study_info.get('alias'):
        return "{}_{}.json".format(study_info.get('alias'), object_type)
    else:
        print('ERROR: No study name found in study_info.')


def write_json_file(json_object, object_type, sub_info):
    json_file_name = usi_object_file_name(object_type, sub_info)
    with codecs.open(json_file_name, 'w', encoding='utf-8') as jf:
        json.dump(json_object, jf)



