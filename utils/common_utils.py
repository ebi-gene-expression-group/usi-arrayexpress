import logging
from datetime import datetime
import os.path


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

