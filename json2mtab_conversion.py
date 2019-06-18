#!/usr/bin/env python

"""
This script takes a JSON submission envelope file (with all submittables)
and runs metadata conversion from USI-JSON to MAGE-TAB format.
"""

import argparse
from os import path

from converter import json2dm
from validator.json_schema_validation import validate_submission_json
from utils.common_utils import file_exists, create_logger
from utils.converter_utils import read_json_file


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json',
                        help="Path to the USI-JSON file")
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

    # Create logger for JSON errors
    json_logger = create_logger(path.dirname(json_file), process_name, path.basename(json_file),
                                logger_name="JSON")
    # Validate the submission JSON against the full ArrayExpress submission schema
    validate_submission_json(json_file, logger=json_logger)

    json_data = read_json_file(json_file)
    #sub = json2dm.data_objects_from_json(json_data, json_file)
    #print(sub.sample)

    wd = path.dirname(path.realpath(__file__))
    mapping_file = path.join(wd, "converter", "config", "mapping_ae-usi_to_datamodel.json")
    ae_converter = json2dm.JSONConverter(mapping_file, import_key="ae")
    sub = ae_converter.convert(json_data)


if __name__ == '__main__':
    main()
