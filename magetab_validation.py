"""
This script takes an IDF file as input and runs validation of the metadata in the common datamodel.
"""

import argparse
import os
import logging

from utils.common_utils import create_logger
from utils.converter_utils import get_sdrf_path
from converter.converting import data_objects_from_magetab

import validator.metadata_validation as mv


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="name of MAGE-TAB IDF file")

    args = parser.parse_args()

    return args.idf


def main():
    process_name = "magetab_validation"

    idf_file = parse_args()

    # Create logger
    current_dir, idf_file_name = os.path.split(idf_file)
    logger = create_logger(current_dir, process_name, idf_file_name)
    # Add handler to direct output to screen in addition to log file
    stream_hdlr = logging.StreamHandler()
    logger.addHandler(stream_hdlr)

    # Get path to SDRF file
    sdrf_file_path = get_sdrf_path(idf_file, logger)

    # Read in MAGE-TAB and convert to common data model
    sub = data_objects_from_magetab(idf_file, sdrf_file_path)

    # Validate metadata in common data model
    mv.run_project_checks(sub, logger)
    mv.run_study_checks(sub, logger)
    mv.run_protocol_checks(sub, logger)
    mv.run_sample_checks(sub, logger)


if __name__ == '__main__':
    main()
