"""This test fetches all JSON files from the usijson directory and validates them
against the schema with matching name in the json_schemas folder."""

import unittest
import os
import jsonschema
import json
import re

from converter.magetab2dm import mtab2usi_conversion


class TestJsonSchemasAndObjects(unittest.TestCase):

    def setUp(self):
        test_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(test_dir)
        self.schema_dir = os.path.join(base_dir, "json_schemas")
        self.mage_dir = os.path.join(test_dir, "test_data")
        self.json_dir = os.path.join(test_dir, "test_data", "usijson")
        # Re-write JSON files
        magetab_files = [os.path.join(self.mage_dir, f) for f in os.listdir(self.mage_dir) if f.endswith(".idf.txt")]
        for idf_file in magetab_files:
            mtab2usi_conversion(idf_file)

    def test_validate_json_schemas(self):

        schema_files = [os.path.join(self.schema_dir, f) for f in os.listdir(self.schema_dir) if f.endswith(".json")]

        for sf in schema_files:
            print("Testing schema: {}".format(sf))
            with open(sf) as js_file:
                schema = json.load(js_file)
                validator = jsonschema.Draft4Validator(schema)
                validator.check_schema(schema)

    def test_validate_json_object_against_schema(self):

        test_objects = [os.path.join(self.json_dir, f) for f in os.listdir(self.json_dir) if f.endswith(".json")]

        types = ("project", "study", "protocol", "sample", "assay_data", "analysis")
        for test_file_name in test_objects:
            with open(test_file_name) as tf:
                test_object = json.load(tf)
                for schema_name in types:
                    # Find the fitting schema
                    if re.search(schema_name, test_file_name):
                        print("Testing {} against {} schema".format(test_file_name, schema_name))
                        with open(os.path.join(self.schema_dir, "arrayexpress_{}_schema.json".format(schema_name))) as so:
                            schema_object = json.load(so)
                            # The file contains a list and only the individual elements are described by the schema
                            if isinstance(test_object, list):
                                for entry in test_object:
                                    jsonschema.validate(entry, schema_object)
                            else:
                                jsonschema.validate(test_object, schema_object)


if __name__ == '__main__':
    unittest.main()
