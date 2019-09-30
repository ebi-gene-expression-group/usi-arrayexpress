
import re
from collections import OrderedDict

from datamodel.submission import Submission
from datamodel.sample import Sample
from datamodel.protocol import Protocol
from datamodel.study import Study
from datamodel.project import Project
from datamodel.data import AssayData, Analysis
from datamodel.assay import SeqAssay, SingleCellAssay, MicroarrayAssay
from datamodel.components import Attribute, Unit
from utils.converter_utils import guess_submission_type_from_study, get_term_from_url
from utils.common_utils import get_ontology_from_term_url, get_term_parent


class JSONConverter:

    def __init__(self, mapping, import_key="ae"):
        """
        The conversion relies on the instructions in the mapping file. This is a JSON schema representation
        of the datamodel that describes the import for each attribute.
        The following keys are expected:
        "import: the name of the import route, e.g. "ae" for a USI ArrayExpress submission
        Each import element should have:
            "path": describes the way through the JSON to find the element to convert, as a list of element names
            "method": contain the names of the function that should be applied to a given object The "path"

        :param mapping: JSON schema object describing the import strategy for each element in the datamodel
        :param import_key: the key to select the import strategy from the mapping file
        """
        self.mapping = mapping
        self.import_key = import_key
        self.unit_types = {}  # Will be used to store already discovered unit types

    def convert_submission(self, envelope_json, submission_type=None, source_file_name=None):
        """
        Converter that takes a JSON as input and converts it to a Submission class object
        based on the specifications in the mapping file
        :param envelope_json: Input JSON with all submittable objects
        :param submission_type: microarray, sequencing or singlecell
        :param source_file_name: name of the original metadata file
        :return: Submission object
        """
        # We only take the first project in the list
        project_json = next(iter(envelope_json.get("projects", [])), {})
        project = Project(**self.convert_submittable(project_json, "project"))

        # We only take the first study in the list
        study_json = next(iter(envelope_json.get("studies", [])), {})
        study = Study(**self.convert_submittable(study_json, "study"))

        protocol_json = envelope_json.get("protocols", [])
        protocols = [Protocol(**self.convert_submittable(p, "protocol")) for p in protocol_json]

        samples_json = envelope_json.get("samples", [])
        samples = [Sample(**self.convert_submittable(s, "sample")) for s in samples_json]

        # To pick the right assay sub-type we need to know the submission type
        if not submission_type:
            submission_type = study.submission_type
            # Try to guess the experiment type from the experiment type
            if not submission_type:
                submission_type = guess_submission_type_from_study(study)
                if not submission_type:
                    raise Exception("Cannot identify submission type.")

        assays = []
        if submission_type == "microarray":
            assays = [MicroarrayAssay(**self.convert_submittable(a, "microarray_assay"))
                      for a in envelope_json.get("assays", [])]
        elif submission_type == "sequencing":
            assays = [SeqAssay(**self.convert_submittable(a, "sequencing_assay"))
                      for a in envelope_json.get("assays", [])]
        elif submission_type == "singlecell":
            assays = [SingleCellAssay(**self.convert_submittable(a, "singlecell_assay"))
                      for a in envelope_json.get("assays", [])]

        assay_data = [AssayData(**self.convert_submittable(ad, "assay_data"))
                      for ad in envelope_json.get("assayData", [])]

        analysis = [Analysis(**self.convert_submittable(a, "analysis"))
                    for a in envelope_json.get("analyses", [])]

        sub_info = {
            "team": envelope_json.get("submission", {}).get("team", {}).get("name"),
            "alias": envelope_json.get("submission", {}).get("id", {}),
            "submission_type": submission_type,
            "metadata": source_file_name
        }

        submission = Submission(sub_info, project, study, protocols, samples, assays, assay_data, analysis)
        return submission

    def convert_submittable(self, submittable_object, submittable_name):
        """
        Use the mapping instructions and convert the different attributes accordingly.
        :param submittable_object: part of a JSON submission envelope containing the details for a single "submittable"
        :param submittable_name: name of the submittable as referenced in the mapping file
        :return: dictionary of the datamodel class attributes with the values to be inserted
        """
        submittable_attributes = OrderedDict()

        # Go through attributes in the config
        for attribute, attribute_info in self.mapping.get(submittable_name, {}).items():
            # Get information how to convert
            convert_function = None
            mapping_info = attribute_info.get("import", {}).get(self.import_key, {})
            path = mapping_info.get("path", [])
            method = mapping_info.get("method", "")
            translation = mapping_info.get("translation", {})
            if method:
                convert_function = getattr(self, method)
            if path and convert_function:
                target_object = self.interpret_path(path, submittable_object)
                if isinstance(target_object, list):
                    if attribute_info.get("type") in ["string", "attribute_object"]:
                        # Take the first entry
                        submittable_attributes[attribute] = next(iter([convert_function(o, translation=translation)
                                                                       for o in target_object]))
                    else:
                        submittable_attributes[attribute] = [convert_function(o, translation=translation)
                                                             for o in target_object]
                elif target_object:
                    submittable_attributes[attribute] = convert_function(target_object, translation=translation)

        return submittable_attributes

    def interpret_path(self, path, json_object):
        """
        Return the object in the JSON reached by traversing through the elements in the path.
        :param path: list of element names
        :param json_object: JSON instance
        :return: JSON object
        """
        if len(path) == 1:
            return json_object.get(path[0])
        else:
            return self.interpret_path(path[1:], json_object.get(path[0], {}))

    # The following are converting functions that do the conversion from JSON sub-elements
    # For each data object the function in the mapping file (under import > method) is called
    
    def import_publication(self, element, translation={}):
        return self.convert_submittable(element, "publication")

    def import_contact(self, element, translation={}):
        return self.convert_submittable(element, "contact")

    def import_file(self, element, translation={}):
        return self.convert_submittable(element, "data_file")

    def import_lib_attribs(self, element, translation={}):
        return self.convert_submittable(element, "lib_attribs")

    def generate_attribute_from_json(self, element, translation={}):
        """
        Convert a USI-JSON formatted attribute to a datamodel.Attribute object
        :param element: the USI-JSON attribute object
        :param translation: (optional) a dictionary with translations for controlled terms
        :return: datamodel.Attribute object
        """
        term_accession = None
        term_source = None
        unit = None

        # Try to translate the value or if not in the look-up return the value string directly
        value = translation.get(element.get("value"), element.get("value"))

        terms = element.get("terms", [])
        if terms:
            # We only take the first term
            term = next(iter(terms))
            term_accession = get_term_from_url(term.get("url"))
            term_source = get_ontology_from_term_url(term.get("url"))
        unit_value = element.get("units")
        if unit_value:
            # Generate unit type as this is not a field in USI's unit model, also remove "derived" from the label
            unit_type = self.get_unit_type(unit_value)
            # USI does not support ontology annotations for unit terms, initialising with None
            unit = Unit(value=unit_value,
                        unit_type=unit_type,
                        term_accession=None,
                        term_source=None)

        return Attribute(value=value,
                         unit=unit,
                         term_accession=term_accession,
                         term_source=term_source)

    def get_unit_type(self, unit_value):
        """Look up a unit term in EFO and return the parent term as 'unit type'.
        Units which have already been looked up previously are stored in the converter's unit_types dictionary,
        to reduce number of requests being made."""
        if unit_value not in self.unit_types:
            unit_type = re.sub("\\s*derived\\s*", "", get_term_parent("efo", unit_value), 1)
            self.unit_types[unit_value] = unit_type
            return unit_type
        else:
            return self.unit_types[unit_value]

    @staticmethod
    def get_reference_value_from_json(element, translation={}):
        """References can contain an accession or a combination of alias and team.
        We prefer the accession if there is one, otherwise use the alias (which should be
        unique within the submission)."""
        if element and element.get("accession"):
            return element.get("accession")
        elif element and element.get("alias"):
            return element.get("alias")

    def get_reference_usage(self, element, translation={}):
        reference = element.get("protocolRef", element.get("sampleRef", {}))
        return self.get_reference_value_from_json(reference)

    @staticmethod
    def import_string(element, translation={}):
        """Return the string value with optional translation according to controlled vocabulary."""
        return translation.get(element, str(element))

    @staticmethod
    def get_string_from_attribute(element, translation={}):
        """Return the value of an attribute as string."""
        return translation.get(str(element.get("value")), str(element.get("value")))

    def generate_sample_attribute_dict(self, sample_attributes, translation={}):

        return OrderedDict([(category, self.generate_attribute_from_json(attribute[0]))
                            for category, attribute in sample_attributes.items()])

    def generate_attribute_from_string(self, element, translation={}):
        return Attribute(value=self.import_string(element, translation))

    def get_list_from_string(self, element, translation={}):
        return [self.import_string(element, translation)]

    @staticmethod
    def get_first_object_from_list(input_json,  translation={}):
        """Return the first object from a list or an empty dictionary if the list is empty."""
        return next(iter(input_json), [])

    @staticmethod
    def import_as_is(input_json, translation={}):
        """Return the object as is."""
        return input_json



