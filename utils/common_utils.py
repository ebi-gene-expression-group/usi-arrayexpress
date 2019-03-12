import logging
import os.path
import urllib
import requests
import json
import sys

from datetime import datetime


def create_logger(working_dir, process_name, object_name, log_level=20, logger_name=__name__):

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

    base_url = "https://www.ebi.ac.uk/ols/api/ontologies"

    url_encoded = urllib.parse.quote(urllib.parse.quote(term_url, safe=""))
    param = {'size': 200}
    api_url = "{}/{}/terms/{}/descendants".format(base_url, ontology, url_encoded)

    data = download_json(logger, api_url, param)

    if data:
        try:
            for d in data["_embedded"]["terms"]:
                efo_children.add(d["label"])
            return efo_children
        except KeyError:
            logger.error("Failed to receive valid response from {}.".format(api_url))


def get_ena_library_terms_via_usi(logger):
    """Read ENA's controlled vocabulary using USI's API and
    return dictionary of the field names with the allowed values."""

    url = "https://submission-dev.ebi.ac.uk/api/dataTypes/sequencingExperiments"
    data = download_json(logger, url)
    if data:
        return {field: description["items"]["properties"]["value"].get("enum", [])
                for field, description in data["validationSchema"]["properties"]["attributes"]["properties"].items()}


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



