#!/usr/bin/env python

"""
This script takes a JSON file and runs validation against a JSON schema. If no schema file is given,
 it checks against the ArrayExpress submission schema.
"""

import argparse

from validator.json_schema_validation import validate_submission_json
from utils.common_utils import file_exists


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json',
                        help="Path to JSON file")
    parser.add_argument('-s', '--schema',
                        help="Path to the JSON schema file")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    json_file = args.json
    file_exists(json_file)

    try:
        if args.schema:
            schema_file = args.schema
            file_exists(schema_file)
            # Validate the JSON against the provided schema
            validate_submission_json(json_file, schema_file=schema_file)
        else:
            # Validate the JSON against the full ArrayExpress submission schema
            validate_submission_json(json_file)
    except Exception as e:
        print("ERROR: Cannot read or validate the JSON input\n{}".format(e))


if __name__ == '__main__':
    main()
