"""This tests the loading of schema files and working validation of a submission against a schema"""

from unittest import TestCase

from validator.json_schema_validation import load_submission_schema


class TestLoadingSchema(TestCase):

    def test_load_submission_schema(self):
        # Call schema loading function that reads the schema files in a fixed location
        schema_file, schema = load_submission_schema()
        # Check that we have non null output
        self.assertIsNotNone(schema_file)
        self.assertIsNotNone(schema)

    def test_reading_schemas(self):
        # Call the schema loading function
        schema_file, schema = load_submission_schema()
        # Check the schema content can be parsed
        self.assertIn("samples", schema["properties"])


