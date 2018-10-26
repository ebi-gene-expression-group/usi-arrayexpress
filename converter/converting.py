"""Helper functions for converting from MAGE-TAB to USI JSON and back"""

import re
import json

from collections import OrderedDict
import json

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

def generate_usi_assay_object(assay, study_info):

    assay_object = OrderedDict()

    assay_object['alias'] = assay.alias


    assay_object['attributes'] = []

    assay_object['studyRef'] = generate_usi_ref_object(study_info.get('alias'), study_info)

    sample_refs = assay.samplerefs
    # Can be one or more samples linked to one assay, but schema is expecting an array
    if not isinstance(sample_refs, list):
        sample_refs = [sample_refs]
    assay_object['sampleUses'] = [generate_usi_ref_object(ref, study_info) for ref in sample_refs]


    return assay_object


def usi_object_file_name(object_type, study_info):

    if study_info.get('accession'):
        return "{}_{}.json".format(study_info('accession'), object_type)
    elif study_info.get('alias'):
        return "{}_{}.json".format(study_info.get('alias'), object_type)
    else:
        print('ERROR: No study name found in study_info.')


def generate_usi_attribute_object(attribute_info):
    """ "attributes": {
            "description": "Attributes for describing a submittable.",
            "type": "object",
            "properties": {},
            "patternProperties": {
                "^.*$": {
                    "type": "array",
                    "minItems": 1,
                    "items": {
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

    pass


def generate_usi_ref_object(ref, study_info):
    """ Schema definitions

    "assayRefs": {
            "description": "Reference(s) to assay(s).",
            "type": "array",
            "items": { "$ref": "#/definitions/submittableRef" }
        }

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
        ]
    }
    "required": [ "alias" ]"""

    ref_object = dict()
    ref_object['alias'] = str(ref)
    ref_object['team'] = study_info.get('team')

    if re.match('^[A-Z]-[A-Z]{4}-[0-9]+', ref):
        ref_object['accession'] = str(ref)

    return ref_object


