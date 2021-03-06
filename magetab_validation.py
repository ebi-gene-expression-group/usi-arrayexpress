#!/usr/bin/env python

"""
This script takes an IDF file as input and runs validation of the metadata in the common datamodel.
"""

import argparse
import os

from utils.common_utils import create_logger, file_exists
from utils.converter_utils import get_sdrf_path, guess_submission_type_from_sdrf, guess_submission_type_from_idf, \
    read_sdrf_file, read_idf_file
from converter.magetab2dm import data_objects_from_magetab

import validator.magetab_prevalidation as pre
import validator.metadata_validation as mv


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="Path to the MAGE-TAB IDF file")
    parser.add_argument('-d', "--data_dir", default="",
                        help="Path to the directory with SDRF and data files")
    parser.add_argument('-v', '--verbose', action='store_const', const=10, default=20,
                        help="Option to output detailed logging (debug level).")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-sc', '--singlecell', action='store_const', const="singlecell", dest='submission_type',
                       help="Force submission type to be 'singlecell'")
    group.add_argument('-seq', '--sequencing', action='store_const', const="sequencing", dest='submission_type',
                       help="Force submission type to be 'sequencing'")
    group.add_argument('-ma', '--microarray', action='store_const', const="microarray", dest='submission_type',
                       help="Force submission type to be 'microarray'")

    args = parser.parse_args()

    return args


def main():
    process_name = "magetab_validation"

    args = parse_args()
    idf_file, data_dir, logging_level = args.idf, args.data_dir, args.verbose
    submission_type = args.submission_type

    # Exit if IDF file doesn't exist
    file_exists(idf_file)

    # Create logger
    current_dir, idf_file_name = os.path.split(idf_file)
    logger = create_logger(current_dir, process_name, idf_file_name, logger_name="Validation", log_level=logging_level)

    # Get path to SDRF file
    sdrf_file_path = get_sdrf_path(idf_file, logger, data_dir)

    # Read IDF/SDRF and get submission type
    idf_dict = read_idf_file(idf_file)
    sdrf_data, header, header_dict = read_sdrf_file(sdrf_file_path)

    # Set submission type
    if not submission_type:
        submission_type = guess_submission_type_from_sdrf(sdrf_data, header, header_dict)
        if not submission_type:
            submission_type = guess_submission_type_from_idf(idf_dict)
        logger.info("Detected submission type: {}".format(submission_type))
    else:
        logger.info("Setting submission type to \"{}\"".format(submission_type))

    # Logger for prevalidation (create new to show different styling)
    mtab_logger = create_logger(current_dir, process_name, idf_file_name, logger_name="MAGE-TAB", log_level=logging_level)

    # Perform prevalidation checks on MAGE-TAB format
    pre.idf_prevalidation(idf_dict, mtab_logger)
    pre.sdrf_prevalidation(sdrf_data, header, header_dict, submission_type, mtab_logger)

    # Read in MAGE-TAB and convert to common data model
    sub = data_objects_from_magetab(idf_file, sdrf_file_path, submission_type)

    # Collect error codes
    error_codes = []

    # Logger for metadata validation output
    metadata_logger = create_logger(current_dir, process_name, idf_file_name, logger_name="Metadata", log_level=logging_level)

    # Validate metadata in common data model
    error_codes.extend(mv.run_project_checks(sub, metadata_logger))
    error_codes.extend(mv.run_study_checks(sub, metadata_logger))
    error_codes.extend(mv.run_protocol_checks(sub, metadata_logger))
    error_codes.extend(mv.run_sample_checks(sub, metadata_logger))
    error_codes.extend(mv.run_assay_checks(sub, metadata_logger))
    error_codes.extend(mv.run_file_checks(sub, metadata_logger))
    if submission_type == "singlecell":
        error_codes.extend(mv.run_singlecell_checks(sub, metadata_logger))

    if error_codes:
        logger.info("Validation finished with the following error codes: \n{}".format("\n".join(set(error_codes))))
    else:
        logger.info("Validation was successful!")


if __name__ == '__main__':
    main()
