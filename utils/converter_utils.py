import codecs
import json
import logging
import os
import pkg_resources
import re

from collections import OrderedDict, defaultdict

from utils.eutils import esearch


USI_JSON_DIRECTORY = "usijson"
SDRF_FILE_NAME_REGEX = r"^\s*SDRF\s*File"
DEFAULT_DATA_DIRECTORY = "unpacked"


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
    return get_controlled_vocabulary(category, "ontology")


def get_controlled_vocabulary(category, resource="translations"):
    """Read the json with controlled vocab and return the dict for the given category.
    The resource parameter specifies which file to read."""
    resource_package = __name__
    if resource == "ontology":
        resource_path = "ontology_terms.json"
    elif resource == "magetab":
        resource_path = "magetab_fields.json"
    else:
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


def get_taxon(organism, logger=logging.getLogger()):
    """Return the NCBI taxonomy ID for a given species name."""

    if organism and organism not in organism_lookup:
        # If we have more than one organism mixed in one sample - in the case assign the 'mixed
        # sample' taxon_id (c.f. https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1427524)
        if re.search(r" and | \+ ", organism):
            return 1427524
        logger.info("Looking up species in NCBI taxonomy. Please wait...")
        db = 'taxonomy'
        a = esearch(db=db, term=organism)
        try:
            taxon_id = int(a['esearchresult']['idlist'][0])
            organism_lookup[organism] = taxon_id
            return taxon_id
        except Exception as e:
            logger.error("Failed to retrieve organism data from ENA taxonomy service for {} due to {}".format(organism, str(e)))
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


def get_sdrf_path(idf_file_path, logger, data_dir):
    """Read IDF and get the SDRF file name, look for the SDRF in the data directory (i.e. "unpacked")
    or in the same directory as the IDF.

    :param idf_file_path: full or relative path to IDF file
    :param logger: log handler
    :param data_dir: path to folder with SDRF
    :return: path to SDRF file
    """

    current_dir = os.path.dirname(idf_file_path)
    sdrf_file_path = ""
    if not data_dir:
        data_dir = DEFAULT_DATA_DIRECTORY
    # Figure out the name and location of sdrf files
    with codecs.open(idf_file_path, 'rU', encoding='utf-8') as f:
        # U flag makes it portable across in unix and windows (\n and \r\n are treated the same)
        for line in f:
            if re.search(SDRF_FILE_NAME_REGEX, line):
                sdrf_file_name = line.split('\t')[1].strip()
                data_path = os.path.join(current_dir, data_dir)
                if os.path.exists(data_path):
                    sdrf_file_path = os.path.join(data_path, sdrf_file_name)
                else:
                    sdrf_file_path = os.path.join(current_dir, sdrf_file_name)
    logger.debug("Generated SDRF file path: {}".format(sdrf_file_path))
    if not os.path.exists(sdrf_file_path):
        logger.error("SDRF file {} does not exist".format(sdrf_file_path))

    return sdrf_file_path


def guess_submission_type_from_sdrf(sdrf_data, header, header_dict):
    """ Guess the basic experiment type (microarray or sequencing) from SDRF"""

    if 'arraydesignref' in header_dict or 'labeledextractname' in header_dict:
        return "microarray"
    elif "comment" in header_dict:
        for comment_index in header_dict.get("comment"):
            if get_value(header[comment_index]) == "library construction" \
                    or get_value(header[comment_index]) == "single cell isolation":
                return "singlecell"
    if "technologytype" in header_dict:
        index = header_dict.get("technologytype")
        if len(index) > 0:
            index = index[0]
            if sdrf_data[0][index] == "array assay":
                return "microarray"
            elif sdrf_data[0][index] == "sequencing assay":
                return "sequencing"


def guess_submission_type_from_idf(idf_dict):
    """Based on the experiment type, we can try to infer the basic experiment type
    This returns the type of the first experiment type found. We cannot handle mixed type experiments.
    """
    print(idf_dict)
    if "AEExperimentType" in idf_dict:
        all_types = get_controlled_vocabulary("experiment_type", "ontology")
        print(all_types)
        for exptype in idf_dict["AEExperimentType"]:
            if exptype in all_types["sequencing"]:
                return "sequencing"
            elif exptype in all_types["microarray"]:
                return "microarray"
            elif exptype in all_types["singlecell"]:
                return "singlecell"


def guess_submission_type(idf_file, sdrf_file, logger):
    """Read IDF/SDRF to get submission type"""

    idf_dict = read_idf_file(idf_file)
    sdrf_data, header, header_dict = read_sdrf_file(sdrf_file)
    submission_type = guess_submission_type_from_sdrf(sdrf_data, header, header_dict)
    if not submission_type:
        submission_type = guess_submission_type_from_idf(idf_dict)
    logger.info("Found experiment type: {}".format(submission_type))
    return submission_type, idf_dict,


def get_name(header_string):
    """Return the first part of an SDRF header in lower case and without spaces."""
    no_spaces = header_string.replace(' ', '')
    field_name = no_spaces.split('[')[0]
    return field_name.lower()


def get_value(header_string):
    """Return the value within square brackets of an SDRF header."""
    field_value = header_string.split('[')[-1]
    return field_value.strip(']')


def read_sdrf_file(sdrf_file):
    """
    Read SDRF file and return the table content as nested list,
    the header row as list, and a dictionary of the fields and their indexes
    :param sdrf_file: string, path to SDRF file
    """

    with codecs.open(sdrf_file, encoding='utf-8') as sf:
        header = sf.readline().rstrip().split('\t')
        sdrf_raw = sf.readlines()

    sdrf_list = [x.rstrip('\n').split('\t') for x in sdrf_raw]
    header_dict = defaultdict(list)
    for i, field in enumerate(header):
        short_name = get_name(field)
        header_dict[short_name].append(i)

    return sdrf_list, header, header_dict


def read_idf_file(idf_file):
    """This function reads in an IDF file and determines whether it is a normal or a merged file.
    It then returns the data as a dictionary with the field names as keys and values as list."""
    idf_dict = OrderedDict()
    with codecs.open(idf_file, encoding='utf-8') as fi:
        idf_raw = fi.readlines()

        for row in idf_raw:
            idf_row = row.rstrip('\n').split('\t')
            # Skip empty lines or lines full of just empty spaces/tabs
            if len(''.join(idf_row).strip()) == 0:
                continue
            if re.search(r"^\[SDRF\]", idf_row[0]):
                # idf_file is a combined idf/sdrf file - stop when you get to the beginning of sdrf section
                break
            # In the case of a merged idf/sdrf file - skip the beginning of the idf section)
            if not re.search(r"^\[IDF\]", idf_row[0]):
                # Store label in separate variable
                row_label = idf_row.pop(0)
                # For comments get the value inside square brackets
                if re.search('comment', row_label, flags=re.IGNORECASE):
                    row_label = get_value(row_label)
                # For normal Row labels remove whitespaces
                else:
                    row_label = get_name(row_label)
                # Skip rows that have a label but no values
                if not len(''.join(idf_row).strip()) == 0:
                    # Some comments are duplicated
                    if row_label in idf_dict:
                        # in that case add to existing entry
                        idf_dict[row_label].extend(idf_row)
                    else:
                        # Store values in idf_dict
                        idf_dict[row_label] = idf_row
    return idf_dict