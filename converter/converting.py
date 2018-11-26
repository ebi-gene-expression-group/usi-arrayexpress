"""Helper functions for converting from MAGE-TAB to USI JSON and back"""

import codecs
import json
from collections import OrderedDict, defaultdict
from utils.converter_utils import is_accession

"""
# Cannot import "common" from fgsubs; it has too many dependencies, so copying this here
def get_url(url):
    data = None
    try:
        r = urllib3.urlopen(url)
        data = json.load(r)
    except urllib2.HTTPError as e:
        print('HTTPError when retrieving url: "' + url + '" : ' + str(e.code))
    except urllib2.URLError as e:
        print('URLError when retrieving url: "' + url + '" : ' + str(e.reason))
    except httplib.HTTPException as e:
        print('HTTPException when retrieving url: "' + url + '"')
    else:
        return data


def get_taxon(organism):
    taxon_to_id = {}
    if organism != '':
        # TODO: Fix urllib.urlencode instead of organism.replace(' ','%20')
        organism_data = get_url(
            'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/suggest-for-search/' + organism.replace(
                ' ', '%20') + '?limit=200')
        if not organism_data:
            if re.search(r" and | \+ ", organism):
                # It looks as if we have more than one organism mixed in one sample - in the case assign the 'mixed
                # sample' taxon_id (c.f. https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1427524) - as
                # per instructions on https://www.ebi.ac.uk/seqdb/confluence/display/GXA/Curation+Look-up
                return 1427524
            else:
                print("Failed to retrieve organism data from ENA taxonomy service for: " + organism)
        else:
            for organism_entry in organism_data:
                if organism_entry['scientificName'].lower() == organism.lower():
                    taxon_id = int(organism_entry['taxId'])
                    return taxon_id
"""


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
    for factor_dict in study.experimental_factor:
        factor_dict["value"] = factor_dict["experimental_factor"]
        attrib_obj = generate_usi_attribute_entry(factor_dict)
        study_attributes["experimental_factor"].extend(attrib_obj)

    for design_dict in study.experimental_design:
        design_dict["value"] = design_dict["experimental_design"]
        study_attributes["experimental_design"].extend(generate_usi_attribute_entry(design_dict))

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
            if unit_info.get("term accession"):
                # USI model does support >1 term URIs for a value but we can't distinguish between
                # the value term and the unit term. For now will only include an ontology URI for
                # a value term (see below).
                pass
        if attribute_info.get("term accession"):
            # TODO: Function that looks up term accessions and returns EFO/OLS URI
            # Using term accession for now
            attribute_object[0]["terms"] = [{"url": attribute_info.get("term accession")}]
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



