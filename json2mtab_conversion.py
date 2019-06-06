#!/usr/bin/env python

"""
This script takes a JSON submission envelope file (with all submittables)
and runs metadata conversionfrom USI-JSON to MAGE-TAB format.
"""

import argparse
import os

from converter import json2dm
from validator.json_schema_validation import validate_submission_json
from utils.common_utils import file_exists, create_logger
from utils.converter_utils import read_json_file


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json',
                        help="Path to the USI-JSON file")
    parser.add_argument('-s', '--schema',
                        help="Path to the JSON schema file")
    parser.add_argument('-v', '--verbose', action='store_const', const=10, default=20,
                        help="Option to output detailed logging (debug level).")
    args = parser.parse_args()

    return args


def main():
    process_name = "json2mtab_conversion"

    args = parse_args()
    json_file = args.json

    # Exit if IDF file doesn't exist
    file_exists(json_file)

    #logger = create_logger(os.path.dirname(json_file), process_name, json_file)

    validate_submission_json(json_file)

    #sub = json2dm.data_objects_from_json(json_data, json_file)
    #print(sub.sample)


if __name__ == '__main__':
    main()
