import argparse
import logging
from os import path

import validator.metadata_validation as mv
from converter import json2dm
from utils.common_utils import file_exists
from utils.converter_utils import read_json_file


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json',
                        help="Path to the USI-JSON file")

    parser.add_argument('-v', '--verbose', action='store_const', const=10, default=20,
                        help="Option to output detailed logging (debug level).")
    args = parser.parse_args()

    return args


args = parse_args()
    json_file = args.json

    json_logger = logging.getLogger("AE-metadata-check")


    # Exit if IDF file doesn't exist
    file_exists(json_file)

    # Convert data to data model
    json_data = read_json_file(json_file)
    wd = path.dirname(path.realpath(__file__))
    mapping_file = path.join(wd, "converter", "config", "mapping_ae-usi_to_datamodel.json")
    mapping = read_json_file(mapping_file)
    ae_converter = json2dm.JSONConverter(mapping, import_key="ae")
    sub = ae_converter.convert_usi_sub(json_data, source_file_name=json_file)

    # Validate metadata in common data model
    error_codes = []
    error_codes.extend(mv.run_project_checks(sub, json_logger))
    error_codes.extend(mv.run_study_checks(sub, json_logger))
    error_codes.extend(mv.run_protocol_checks(sub, json_logger))
    error_codes.extend(mv.run_sample_checks(sub, json_logger))
    error_codes.extend(mv.run_assay_checks(sub, json_logger))
    error_codes.extend(mv.run_file_checks(sub, json_logger))
    if sub.info.get("submission_type") == "singlecell":
        error_codes.extend(mv.run_singlecell_checks(sub, json_logger))
