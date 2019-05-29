"""Module to run JSON schema validation of a USI submission metadata file against the ArrayExpress submission schema"""

import os
import jsonschema

import json_schemas
from utils.converter_utils import read_json_file


def load_submission_schema():
    """Load the ArrayExpress submission schema as JSON"""

    # Construct path to schema and load as json
    schema_path = os.path.dirname(json_schemas.__file__)
    schema_file = os.path.join(schema_path, "arrayexpress_submission.json")
    schema = read_json_file(schema_file)

    return schema_file, schema


def validate_submission_json(submission_json):

    schema_file, schema = load_submission_schema()

    # Create validator with resolver to help locate the referenced submittable schemas with absolute path
    resolver = jsonschema.RefResolver("file://" + schema_file, schema)
    validator = jsonschema.Draft4Validator(schema, resolver=resolver)

    # Check validity of schema first
    # We shouldn't need this but for now, running this test before calling the checks on the object
    validator.check_schema(schema)

    # Validate the submission JSON against the submission schema
    errors = validator.iter_errors(submission_json)
    tree = jsonschema.ErrorTree(errors)
    print("Total number of errors: " + str(tree.total_errors))
    print(tree["assayData"][0]["files"][0]["checksum"].errors)
    print(tree["protocols"][0].errors)
    recurse_through_errors(tree)


def recurse_through_errors(error_tree):
    """The error tree arranges the error messages in the same structure as the validated instance.
    This allows to get the traceback of which exact instance is failing validation."""

    if error_tree.errors:
        # Found a non-empty error dict, this is the lowest level
        print(error_tree.errors)
    else:
        # Go a level deeper and look for errors
        for branch in error_tree:
            print(branch)
            recurse_through_errors(error_tree[branch])


