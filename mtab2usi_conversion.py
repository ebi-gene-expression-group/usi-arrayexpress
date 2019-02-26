#!/usr/bin/env python

"""
This script takes an IDF file as input and runs metadata conversion from MAGE-TAB to USI-JSON format.
"""

import argparse
import os

from utils.common_utils import create_logger, file_exists
from utils.converter_utils import get_sdrf_path
from converter.converting import data_objects_from_magetab, datamodel2json_conversion


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="name of MAGE-TAB IDF file")

    args = parser.parse_args()

    return args.idf


def main():
    process_name = "mtab2usi_conversion"

    idf_file = parse_args()
    file_exists(idf_file)

    # Create logger
    current_dir, idf_file_name = os.path.split(idf_file)
    logger = create_logger(current_dir, process_name, idf_file_name)

    # Get path to SDRF file
    sdrf_file_path = get_sdrf_path(idf_file, logger)

    # Read in MAGE-TAB and convert to common data model

    sub = data_objects_from_magetab(idf_file, sdrf_file_path)

    # Dump data in common data model as USI-JSON files
    datamodel2json_conversion(sub, current_dir, logger)


if __name__ == '__main__':
    main()
