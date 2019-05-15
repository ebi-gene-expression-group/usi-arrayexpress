""" Script for testing the MAGE-TAB writer, where the input is a set of MAGE-TAB files
which is being read into the datamodel and written back out as MAGE-TAB.

Mandatory input is the IDF file path. The SDRF is expected to be located in the same folder.
Optional parameters are the output directory (-o) and verbose logging (-v).
"""

import argparse
import os

import utils.converter_utils
from utils.common_utils import create_logger, file_exists
from utils.converter_utils import get_sdrf_path, guess_submission_type
from converter.converting import data_objects_from_magetab
from converter.datamodel2magetab import generate_idf, generate_sdrf, write_sdrf_file

process_name = "magetab_writer"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="Path to the MAGE-TAB IDF file")
    parser.add_argument('-o', "--outdir",
                        help="Path where to write the new MAGE-TAB files")
    parser.add_argument('-v', '--verbose', action='store_const', const=10, default=20,
                        help="Option to output detailed logging (debug level)")

    return parser.parse_args()


def main():
    args = parse_args()

    # Check if input file exists
    idf_file = args.idf
    file_exists(idf_file)

    # Create logger
    current_dir, idf_file_name = os.path.split(idf_file)
    logger = create_logger(current_dir, process_name, idf_file_name,
                           log_level=args.verbose, logger_name="MAGE-to-MAGE")

    # Get path to SDRF file
    sdrf_file_path = get_sdrf_path(idf_file, logger)

    # Get correct submission type
    sub_type, idf_dict = guess_submission_type(idf_file, sdrf_file_path, logger)

    # Read in MAGE-TAB and convert to common data model
    sub = data_objects_from_magetab(idf_file, sdrf_file_path, sub_type)

    # New file paths
    if args.outdir:
        new_idf_file = os.path.join(args.outdir, os.path.basename(idf_file) + "_new.txt")
        new_sdrf_file = os.path.join(args.outdir, os.path.basename(sdrf_file_path) + "_new.txt")
    else:
        new_idf_file = idf_file + "_new.txt"
        new_sdrf_file = sdrf_file_path + "_new.txt"

    # Generate IDF dictionary
    idf = generate_idf(sub)
    # Write out a new IDF file
    utils.converter_utils.dict_to_vertcial_table(idf, new_idf_file, logger)

    # Generate SDRF: Output is a pandas dataframe
    raw_out = generate_sdrf(sub)

    # Rename the columns to the new header list, created by applying a function
    # to "de-uniquify" the header fields, and write new SDRF file
    write_sdrf_file(raw_out, new_sdrf_file, logger)


if __name__ == '__main__':
    main()
