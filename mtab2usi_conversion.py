#!/usr/bin/env python

"""
This script takes an IDF file as input and runs metadata conversion from MAGE-TAB to USI-JSON format.
"""

import argparse
from os.path import isdir, split

from utils.common_utils import create_logger, file_exists
from utils.converter_utils import get_sdrf_path, guess_submission_type
from converter.dm2json import datamodel2json_conversion
from converter.magetab2dm import data_objects_from_magetab


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="name of MAGE-TAB IDF file")
    parser.add_argument('-o', '--outdir',
                        help="Path where to write the JSON file(s)")
    parser.add_argument('-e', '--envelope', action='store_true',
                        help="Option to output only one submission envelope type JSON file")

    args = parser.parse_args()

    return args


def main():
    process_name = "mtab2usi_conversion"

    args = parse_args()
    idf_file = args.idf
    file_exists(idf_file)

    # Output directory
    current_dir, idf_file_name = split(idf_file)
    outdir = current_dir
    if args.outdir and isdir(args.outdir):
        outdir = args.outdir

    # Create logger
    logger = create_logger(current_dir, process_name, idf_file_name)

    # Get path to SDRF file
    sdrf_file_path = get_sdrf_path(idf_file, logger)

    # Get submission type
    submission_type, idf_data = guess_submission_type(idf_file, sdrf_file_path, logger)

    # Read in MAGE-TAB and convert to common data model
    sub = data_objects_from_magetab(idf_file, sdrf_file_path, submission_type)

    # Dump data in common data model as USI-JSON files
    datamodel2json_conversion(sub, outdir, logger, envelope=args.envelope)


if __name__ == '__main__':
    main()
