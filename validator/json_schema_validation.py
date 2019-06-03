"""Module to run JSON schema validation of a USI submission metadata file against the ArrayExpress submission schema"""

import os
import jsonschema

import json_schemas
from utils.converter_utils import read_json_file


def load_submission_schema():
    """Find and load the ArrayExpress submission schema as JSON"""

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

    # Print out all error messages with where and why details
    for e in errors:
        print(format_json_error_message(e))


def recurse_through_errors(error_tree, traceback):
    """The error tree arranges the error messages in the same structure as the validated instance.
    This allows to get the traceback of which exact instance is failing validation."""

    if error_tree.errors:
        # Found a non-empty error dict, this is the lowest level
        print(error_tree.errors)
    else:
        # Go a level deeper and look for errors
        for branch in error_tree:
            traceback.append(branch)
            print(branch)
            print(traceback)
            recurse_through_errors(error_tree[branch], traceback)


def format_traceback_path(path_element_list):
    """Transform a list of JSON instance path elements and return a single string for the error message.
    The traceback can contain number elements which need to be converted to string."""
    path_as_string = [str(p) for p in path_element_list]
    return " > ".join(path_as_string)


def format_json_error_message(error_object):
    """Input is a jsonschema.ValidationError object. Output is a printable string for the error log.

    Using the following attributes from ValidationError:
    absolute_path: deque of all the parents of the instance with the error
    message: the problem
    validator: the type of JSON validation keyword that caused the error (e.g. oneOf, additionalProperties, required)
    context: gives the list of sub-errors in case of oneOf/anyOf
    """
    message = error_object.message
    if error_object.validator == "oneOf" or error_object.validator == "anyOf":
        # Making the error message shorter to not repeat the complete instance
        # And adding all the errors in the sub-schemas that made the oneOf/any check fail
        message = "Instance is not valid under any of the given schemas:" + "\n\t" + \
                  "\n\t".join([format_json_error_message(sub_error) for sub_error in error_object.context])

    error_string = "Error in {}: ({} error) {}".format(
        format_traceback_path(error_object.absolute_path), error_object.validator,  message)

    return error_string

