import codecs
import json
import os
import re
import pkg_resources

from collections import OrderedDict

from utils.eutils import esearch
from converter.parsing import read_sdrf_file, get_value

USI_JSON_DIRECTORY = "usijson"
SDRF_FILE_NAME_REGEX = r"^\s*SDRF\s*File"
DATA_DIRECTORY = "unpacked"


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


def ontology_term(category):
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
        "ARRAYEXPRESS": "^\s*[A-Z]-[A-Z]{4}-[0-9]+\s*$",
        "BIOSAMPLES": "^\s*SAM[END][AG]?[0-9]+\s*$",
        "ENA": "^\s*ER[RXSP][0-9]+\s*$",
        "BIOSTUDIES": "^\s*S-[A-Z]+[0-9]+\s*$",
        "PRIDE": "^\s*PXD[0-9]+\s*$"}

    regex_ebi_accession = "|".join(regex_lookup.values())

    if archive:
        try:
            regex = regex_lookup.get(archive)
            return re.match(regex, accession)
        except KeyError:
            print("Not a valid archive type: {}".format(archive))
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


def get_efo_url(term_accession):
    """Return URI for a given EFO ontology term

    To simplify this, it seems there is an easy pattern to the URL structure:
    ebi/efo for EFO terms, purl.obolibrary.org/obo/ for imported terms
    """

    url_friendly_term = re.sub(':', '_', term_accession)
    if term_accession.startswith("EFO"):
        return "http://www.ebi.ac.uk/efo/" + url_friendly_term
    else:
        return "http://purl.obolibrary.org/obo/" + url_friendly_term


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


def attrib2dict(ob):
    """Get all attributes of an object (ob) that have a value and turn them into a dictionary.
    The attributes will become the keys of the dict and the object values the dict values."""

    attrib_dict = OrderedDict([(a, getattr(ob, a)) for a in vars(ob) if getattr(ob, a)])

    return attrib_dict


def get_sdrf_path(idf_file_path, logger):
    """Read IDF and get the SDRF file name, look for the SDRF in the data directory (i.e. "unpacked")
    or in the same directory as the IDF.

    :param idf_file_path: full or relative path to IDF file
    :param logger: log handler
    :return: path to SDRF file
    """

    current_dir = os.path.dirname(idf_file_path)
    sdrf_file_path = None
    # Figure out the name and location of sdrf files
    with codecs.open(idf_file_path, 'rU', encoding='utf-8') as f:
        # U flag makes it portable across in unix and windows (\n and \r\n are treated the same)
        for line in f:
            if re.search(SDRF_FILE_NAME_REGEX, line):
                sdrf_file_name = line.split('\t')[1].strip()
                if os.path.exists(current_dir + DATA_DIRECTORY):
                    sdrf_file_path = os.path.join(current_dir, DATA_DIRECTORY, sdrf_file_name)
                else:
                    sdrf_file_path = os.path.join(current_dir, sdrf_file_name)
    logger.debug("Generated SDRF file path: {}".format(sdrf_file_path))
    if not os.path.exists(sdrf_file_path):
        logger.error("SDRF file {} does not exist".format(sdrf_file_path))

    return sdrf_file_path


def guess_submission_type_from_sdrf(sdrf_path):
    """ Guess the basic experiment type (microarray or sequencing) from SDRF header"""

    sdrf_data, header, header_dict = read_sdrf_file(sdrf_path)
    if 'arraydesignref' in header_dict or 'labeledextractname' in header_dict:
        return "microarray"
    elif "comment" in header_dict:
        for comment_index in header_dict.get("comment"):
            if get_value(header[comment_index]) == "library construction" \
                    or get_value(header[comment_index]) == "single cell isolation":
                return "singlecell"
    elif "technologytype" in header_dict:
        index = header_dict.get("technologytype")
        if sdrf_data[0][index] == "array assay":
            return "microarray"
        elif sdrf_data[0][index] == "sequencing assay":
            return "sequencing"


def guess_submission_type_from_idf(idf_path):
    pass