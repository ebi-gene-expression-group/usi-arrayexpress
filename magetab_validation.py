#!/usr/bin/env python

"""
This script takes an IDF file as input and runs validation of the metadata in the common datamodel.
"""

import argparse
import os
import logging

from utils.common_utils import create_logger
from utils.converter_utils import get_sdrf_path
from converter.parsing import read_sdrf_file
from converter.converting import data_objects_from_magetab

import validator.magetab_prevalidation as pre
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

    # Perform prevalidation checks on MAGE-TAB format
    is_microarray = True
    sdrf_data, header, header_dict = read_sdrf_file(sdrf_file_path)
    pre.sdrf_prevalidation(sdrf_data, header, header_dict, is_microarray, logger)

    # Read in MAGE-TAB and convert to common data model
    sub = data_objects_from_magetab(idf_file, sdrf_file_path)

    # Collect error codes
    error_codes = []

    # Validate metadata in common data model
    error_codes.extend(mv.run_project_checks(sub, logger))
    error_codes.extend(mv.run_study_checks(sub, logger))
    error_codes.extend(mv.run_protocol_checks(sub, logger))
    error_codes.extend(mv.run_sample_checks(sub, logger))

    if error_codes:
        logger.info("Validation finished with the following error codes: \n{}".format("\n".join(error_codes)))
    else:
        logger.info("Validation was successful!")


if __name__ == '__main__':
    main()
