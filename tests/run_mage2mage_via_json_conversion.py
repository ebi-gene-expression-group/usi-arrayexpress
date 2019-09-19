""" Script for testing the MAGE-TAB to JSON and JSON to MAGE-TAB converter.
The input is a set of MAGE-TAB files which is being read into the datamodel and written back out as JSON envelope.
Then the JSON envelope is converted back to MAGE-TAB.

Mandatory input is the IDF file path. The SDRF is expected to be located in the same folder.
Optional parameters are the output directory (-o) and verbose logging (-v).
"""

import argparse
import json
import os

import pkg_resources

import validator.metadata_validation as mv

from converter.dm2json import datamodel2json_conversion
from converter.json2dm import JSONConverter
from utils.common_utils import create_logger, file_exists, dir_exists
from utils.converter_utils import get_sdrf_path, guess_submission_type, read_json_file, usi_object_file_name, \
    dict_to_vertical_table
from converter.magetab2dm import data_objects_from_magetab
from converter.dm2magetab import generate_idf, generate_sdrf, write_sdrf_file

process_name = "mage2mage-test"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="Path to the MAGE-TAB IDF file")
    parser.add_argument('-o', "--outdir",
                        help="Path where to write the new MAGE-TAB and JSON files")
    parser.add_argument('-v', '--verbose', action='store_const', const=10, default=20,
                        help="Option to output detailed logging (debug level)")
    parser.add_argument('-l', '--logdir',
                        help="The path where to write the log. Default is same directory as IDF file.")

    return parser.parse_args()


def main():
    args = parse_args()

    # Check if input file exists
    idf_file = args.idf
    file_exists(idf_file)

    # Create logger
    current_dir, idf_file_name = os.path.split(idf_file)
    outdir = current_dir
    if args.outdir and dir_exists(args.outdir):
        outdir = args.outdir

    # Central logger for all input files
    logdir = current_dir
    if args.logdir and dir_exists(args.logdir):
        logdir = args.logdir

    logger = create_logger(logdir, process_name, "all_files",
                           log_level=args.verbose, logger_name="MAGE-to-MAGE")
    logger.info("New Experiment: {}".format(idf_file_name))

    # Get path to SDRF file
    sdrf_file_path = get_sdrf_path(idf_file, logger)

    # Get correct submission type
    sub_type, idf_dict = guess_submission_type(idf_file, sdrf_file_path, logger)

    # Read in MAGE-TAB and convert to common data model
    sub = data_objects_from_magetab(idf_file, sdrf_file_path, sub_type)

    # Run metadata validation in data model
    error_codes = []
    error_codes.extend(mv.run_project_checks(sub, logger))
    error_codes.extend(mv.run_study_checks(sub, logger))
    error_codes.extend(mv.run_protocol_checks(sub, logger))
    error_codes.extend(mv.run_sample_checks(sub, logger))
    error_codes.extend(mv.run_assay_checks(sub, logger))
    error_codes.extend(mv.run_file_checks(sub, logger))
    if sub.info.get("submission_type") == "singlecell":
        error_codes.extend(mv.run_singlecell_checks(sub, logger))

    # Write JSON envelope file
    datamodel2json_conversion(sub, outdir, logger, write_envelope=True)

    # Read JSON envelope file
    json_file = os.path.join(outdir, usi_object_file_name("envelope", sub.info))
    logger.debug("JSON file name: {}".format(json_file))
    json_data = read_json_file(json_file)

    # Convert from JSON to data model
    mapping_file = pkg_resources.resource_string("converter", "config/mapping_ae-usi_to_datamodel.json")
    mapping = json.loads(mapping_file)
    ae_converter = JSONConverter(mapping, import_key="ae")
    sub2 = ae_converter.convert_usi_sub(json_data, source_file_name=json_file)

    # Run metadata validation again
    error_codes = []
    error_codes.extend(mv.run_project_checks(sub2, logger))
    error_codes.extend(mv.run_study_checks(sub2, logger))
    error_codes.extend(mv.run_protocol_checks(sub2, logger))
    error_codes.extend(mv.run_sample_checks(sub2, logger))
    error_codes.extend(mv.run_assay_checks(sub2, logger))
    error_codes.extend(mv.run_file_checks(sub2, logger))
    if sub.info.get("submission_type") == "singlecell":
        error_codes.extend(mv.run_singlecell_checks(sub2, logger))

    # Write back to MAGE-TAB

    # New file paths
    new_idf_file = os.path.join(outdir, os.path.basename(idf_file) + "_new.txt")
    new_sdrf_file = os.path.join(outdir, os.path.basename(sdrf_file_path) + "_new.txt")

    # Generate IDF dictionary
    logger.debug("Generating IDF file")
    idf = generate_idf(sub2)
    # Write out a new IDF file
    dict_to_vertical_table(idf, new_idf_file, logger)

    # Generate SDRF: Output is a pandas dataframe
    logger.debug("Generating SDRF file")
    raw_out = generate_sdrf(sub2)

    # Rename the columns to the new header list, created by applying a function
    # to "de-uniquify" the header fields, and write new SDRF file
    write_sdrf_file(raw_out, new_sdrf_file, logger)


if __name__ == '__main__':
    main()
