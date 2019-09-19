import logging
import os.path
import sys
import urllib
from datetime import datetime


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


def file_exists(input_file, logger=None):
    if not os.path.exists(input_file):
        log_message = "Invalid input. File does not exist: {}".format(input_file)
        if logger:
            logger.error(log_message)
        else:
            print(log_message)
        sys.exit()
    return True


def dir_exists(dir_name):
    if not os.path.isdir(dir_name):
        print("Invalid input. Directory does not exist: {}".format(dir_name))
        sys.exit()
    return True


def is_valid_url(url, logger=None):
    """Check if a given URL exists without downloading the page/file

    For HTTP and HTTPS URLs, urllib.requests returns a http.client.HTTPResponse object,
    for FTP URLs it returns a urllib.response.addinfourl object
    """
    try:
        r = urllib.request.urlopen(url)
        logger.debug("Checking {}... Done.".format(url))
        if r:
            return True
    except urllib.error.URLError:
        return False