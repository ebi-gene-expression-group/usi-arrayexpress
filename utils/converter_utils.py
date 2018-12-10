import json
import pkg_resources
import re
import os
import codecs

from utils.eutils import esearch

USI_JSON_DIRECTORY = "usijson"


def read_json_file(filename):
    try:
        with open(filename) as fh:
            data = json.load(fh, encoding="utf-8")
            return data
    except IOError as err:
        print("Cannot import file: %s" % err)
    except ValueError as j_err:
        print("Cannot read JSON file: %s" % j_err)
        raise


def usi_object_file_name(object_type, study_info):

    if study_info.get('accession'):
        return "{}_{}.json".format(study_info.get('accession'), object_type)
    elif study_info.get('alias'):
        return "{}_{}.json".format(study_info.get('alias'), object_type)
    else:
        print('ERROR: No study name found in study_info.')


def write_json_file(wd, json_object, object_type, sub_info):

    json_file_name = usi_object_file_name(object_type, sub_info)
    json_file_path = os.path.join(wd, USI_JSON_DIRECTORY, json_file_name)
    os.makedirs(os.path.dirname(json_file_path), exist_ok=True)
    with codecs.open(json_file_path, 'w', encoding='utf-8') as jf:
        json.dump(json_object, jf)


def ontology_lookup(category):
    """Read the json with expected EFO terms and return the dict for the given category."""
    resource_package = __name__
    resource_path = "ontology_terms.json"
    all_terms = json.loads(pkg_resources.resource_string(resource_package, resource_path))

    return all_terms[category]


def get_controlled_vocabulary(category):
    """Read the json with controlled vocab and return the dict for the given category."""
    resource_package = __name__
    resource_path = "term_translations.json"
    all_terms = json.loads(pkg_resources.resource_string(resource_package, resource_path))

    return all_terms[category]


def remove_duplicates(ref_list):
    """Return a copy of a list with all duplicated values removed."""
    return list(set(ref_list))


def is_accession(accession, archive=None):
    """Return True if the input is a valid accession format from specified EBI archives.
    With the optional parameter the test can be performed against a specific archive only."""

    regex_lookup = {
        "ARRAYEXPRESS": "^[A-Z]-[A-Z]{4}-[0-9]+",
        "BIOSAMPLES": "^SAM[END][AG]?[0-9]+",
        "ENA": "^ER[RXSP][0-9]+$",
        "BIOSTUDIES": "^S-[A-Z]+[0-9]+$"}

    regex_ebi_accession = "|".join(regex_lookup.values())

    if archive:
        try:
            regex = regex_lookup.get(archive)
            return re.match(regex, accession)
        except KeyError:
            print("Not a valid accession type: {}".format(archive))
    else:
        return re.match(regex_ebi_accession, accession)


# To store organisms that we have already looked-up in the taxonomy (this is slow...)
organism_lookup = {}


def get_taxon(organism):
    """Return the NCBI taxonomy ID for a given species name."""
    if organism and organism not in organism_lookup:

        db = 'taxonomy'
        a = esearch(db=db, term=organism)
        try:
            taxon_id = int(a['esearchresult']['idlist'][0])
            organism_lookup[organism] = taxon_id
            return taxon_id
        except IndexError:
            if re.search(r" and | \+ ", organism):
                # It looks as if we have more than one organism mixed in one sample - in the case assign the 'mixed
                # sample' taxon_id (c.f. https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1427524) - as
                # per instructions on https://www.ebi.ac.uk/seqdb/confluence/display/GXA/Curation+Look-up
                return 1427524
            else:
                print("Failed to retrieve organism data from ENA taxonomy service for: " + organism)
        #except KeyError:
        #    time.sleep(10)
        #    return esearch(db, organism)
    else:
        return organism_lookup.get(organism)


def strip_extension(filename):
    """Take a filename string as input and strip off the file extension and any patterns to be ignored"""
    extensions = ('\.fastq\.gz$', '\.fq\.gz$', '\.txt\.gz$', '\.fastq\.bz2$', '\.[a-zA-Z0-9]+$', )
    ignore = ('_001$', )
    filebase = None

    for ext in extensions:
        if re.search(ext, filename):
            filebase = re.sub(ext, '', filename)
            break
    for ip in ignore:
        if re.search(ip, filebase):
            filebase = re.sub(ip, '', filebase)
            return filebase
    if filebase:
        return filebase
    else:
        return filename


