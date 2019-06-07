"""This tests the loading of schema files and working validation of a JSON file against a schema"""

import os
from unittest import TestCase

from validator.json_schema_validation import load_arrayexpress_submission_schema, validate_submission_json


class TestLoadingSchema(TestCase):

    def test_load_submission_schema(self):
        # Call schema loading function that reads the schema files in a fixed location
        schema_file, schema = load_arrayexpress_submission_schema()
        # Check that we have non null output
        self.assertIsNotNone(schema_file)
        self.assertIsNotNone(schema)

    def test_reading_schemas(self):
        # Call the schema loading function
        schema_file, schema = load_arrayexpress_submission_schema()
        # Check the schema content can be parsed
        self.assertIn("samples", schema["properties"])


class TestValidateJson(TestCase):

    def setUp(self):
        # This test uses a very simple test JSON and test schema
        wd = os.path.dirname(os.path.realpath(__file__))
        self.schema_file = os.path.join(wd, 'test_data', 'simple_schema.json')
        self.test_data_file = os.path.join(wd, 'test_data', 'simple_data.json')

    def test_validate_submission_json(self):
        # This function currently has no output, just testing it runs without errors
        validate_submission_json(self.test_data_file, self.schema_file)
