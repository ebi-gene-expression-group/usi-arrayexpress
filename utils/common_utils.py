import logging
import os.path
import urllib
import requests
import json
import sys

from datetime import datetime

from utils.converter_utils import get_term_from_url, get_ontology_from_term


def create_logger(working_dir, process_name, object_name, log_level=20, logger_name=""):

    # Need to give getLogger file_path (name) argument - so that it creates a unique logger for a given file_path
    log_file_name = "{}_{}_{}.log".format(process_name, object_name, datetime.now().strftime('%Y-%m-%d'))
    log_file_path = os.path.join(working_dir, log_file_name)
    logger = logging.getLogger(logger_name)

    # Logger level: 10=DEBUG, 20=INFO, 30=WARNING, 40=ERROR, 50=CRITICAL)
    # If not specified, default is INFO
    logger.setLevel(log_level)

    # This handler is for writing to the log file
    hdlr = logging.FileHandler(log_file_path)
    formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s: %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)

    # Create second handler for simultaneous console logging (with shorter prompt)
    stream_hdlr = logging.StreamHandler()
    formatter = logging.Formatter('%(name)s %(levelname)-8s %(message)s')
    stream_hdlr.setFormatter(formatter)
    logger.addHandler(stream_hdlr)

    return logger


def query_ols(api_url, param, logger):
    """Basic function to query OLS API"""

    base_url = "https://www.ebi.ac.uk/ols/api/"
    url = base_url + api_url
    data = download_json(logger, url, param)
    return data


def url_encode_for_ols(term_url):
    """For OLS API query special characters in term URLs need to be double-encoded,
    e.g. http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0000001
    """
    return urllib.parse.quote(urllib.parse.quote(term_url, safe=""))


def get_ontology_from_term_url(term_url):
    """Return the ontology for a given ontology term URL

    First try finding the term in EFO, if it is not there,
    use the first bit of the term accession

    :param term_url: URL for an ontology term in OLS
    :return: the name of the ontology that the term comes from
    """

    # Try first to look up term in EFO
    url_encoded = url_encode_for_ols(term_url)
    api_url = "ontologies/efo/terms/{}".format(url_encoded)
    data = query_ols(api_url, {}, logging.getLogger())
    if data:
        try:
            return data.get("ontology_prefix")
        except KeyError:
            logging.error("Failed to receive valid response from OLS searching for {}.".format(term_url))
    # If we haven't found anything in EFO, generate the prefix of the term accession
    else:
        return get_ontology_from_term(get_term_from_url(term_url))


def get_term_descendants(ontology, term_url, logger):
    """
    Use OLS API to retrieve all child terms (descendants) of a given term URL

    :param ontology: name of ontology in OLS
    :param term_url: the url of the term
    :param logger: for logging of errors

    :return: set of child terms labels
    """

    logger.propagate = True
    efo_children = set()

    url_encoded = url_encode_for_ols(term_url)
    param = {'size': 200}
    api_url = "ontologies/{}/terms/{}/descendants".format(ontology, url_encoded)

    data = query_ols(api_url, param, logger)

    if data:
        try:
            for d in data["_embedded"]["terms"]:
                efo_children.add(d["label"])
            return efo_children
        except KeyError:
            logger.error("Failed to receive valid response from {}.".format(api_url))


def get_term_parent(ontology, term):
    """Return the label of the parent term of a given ontology (EFO) term."""

    term_url = ols_lookup(ontology, term)
    url_encoded = url_encode_for_ols(term_url)
    print(term_url)
    api_url = "ontologies/{}/terms/{}/parents".format(ontology, url_encoded)
    data = query_ols(api_url, {}, logging.getLogger())
    if data:
        try:
            for d in data["_embedded"]["terms"]:
                # Return the first hit
                return d["label"]
        except KeyError:
            logger.error("Failed to receive valid response from {}.".format(api_url))


def ols_lookup(ontology, term):
    """Get ontology term URL for a given term label."""

    # Just in case we already have a url
    term_encoded = url_encode_for_ols(term)
    api_url = "search?q={{{}}}&ontology={}".format(term_encoded, ontology)
    print(api_url)
    data = query_ols(api_url, {}, logging.getLogger())
    if data:
        try:
            for d in data["response"]["docs"]:
                # Return the first hit
                return d["iri"]
        except KeyError:
            logging.error("Failed to receive valid response from OLS searching for {}.".format(term))


def get_ena_library_terms_via_usi(logger):
    """Read ENA's controlled vocabulary using USI's API and
    return dictionary of the field names that have a set of allowed values (enum)."""

    url = "https://submission-dev.ebi.ac.uk/api/dataTypes/sequencingExperiments"
    data = download_json(logger, url)
    if data:
        return {field: description["items"]["properties"]["value"]["enum"]
                for field, description in data["validationSchema"]["properties"]["attributes"]["properties"].items()
                if description["items"]["properties"]["value"].get("enum")}


def get_ena_instrument_terms_via_usi(logger):
    """Read ENA's controlled vocabulary using USI's API and
    return list of the instrument models that are allowed (enum)."""

    url = "https://submission-dev.ebi.ac.uk/api/dataTypes/sequencingExperiments"
    data = download_json(logger, url)
    instruments = []
    if data:
        for x in data["validationSchema"]["properties"]["attributes"]["oneOf"]:
            instruments.extend(x["properties"]["instrument_model"]["items"]["properties"]["value"]["enum"])
        return instruments


def download_json(logger, url, parameters=None):
    """Basic function to retrieve URL and return JSON object."""

    logger.debug("Calling: " + url)
    r = requests.get(url, params=parameters)
    if r.status_code != 200:
        logger.error("Failed to receive response from {}. Got error: {}.".format(url, r.status_code))
    else:
        return json.loads(r.text)


def file_exists(input_file):
    if not os.path.exists(input_file):
        print("Invalid input. File does not exist: {}".format(input_file))
        sys.exit()


if __name__ == '__main__':
    # For testing ontology term retrieval
    wd = os.path.dirname(os.path.realpath(__file__))
    logger = create_logger(wd, "testing" "common_utils", 10)
    uri = "http://purl.obolibrary.org/obo/UO_0000000"
    ontology = "EFO"
    terms = get_term_descendants(ontology, uri, logger)



