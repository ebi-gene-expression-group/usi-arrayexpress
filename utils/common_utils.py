import logging
import os.path
import urllib
import requests
import json
import sys

from datetime import datetime


def create_logger(working_dir, process_name, object_name, log_level=20):

    # Need to give getLogger file_path (name) argument - so that it creates a unique logger for a given file_path
    log_file_name = "{}_{}_{}.log".format(process_name, object_name, datetime.now().strftime('%Y-%m-%d'))
    log_file_path = os.path.join(working_dir, log_file_name)
    logger = logging.getLogger(log_file_path)

    # Logger level: 10=DEBUG, 20=INFO, 30=WARNING, 40=ERROR, 50=CRITICAL)
    # If not specified, default is INFO
    logger.setLevel(log_level)

    hdlr = logging.FileHandler(log_file_path)
    formatter = logging.Formatter('%(asctime)s %(levelname)s : %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)

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

    logger.debug("Calling: " + api_url)
    r = requests.get(api_url, params=param)
    data = json.loads(r.text)

    if r.status_code != 200:
        logger.error("Failed to receive response from {}. Got error: {}.".format(api_url, r.status_code))
    else:
        try:
            for d in data["_embedded"]["terms"]:
                efo_children.add(d["label"])
            return efo_children
        except KeyError:
            logger.error("Failed to receive valid response from {}.".format(api_url))


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



