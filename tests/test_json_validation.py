import os
from unittest import TestCase

from validator.json_schema_validation import validate_submission_json


class TestValidateJson(TestCase):

    def setUp(self):
        wd = os.path.dirname(os.path.realpath(__file__))
        self.schema_file = os.path.join(wd, 'test_data', 'simple_test_schema.json')
        self.test_data_file = os.path.join(wd, 'test_data', 'simple_test_json.json')

    def test_validate_submission_json(self):
        validate_submission_json(self.test_data_file, self.schema_file)
