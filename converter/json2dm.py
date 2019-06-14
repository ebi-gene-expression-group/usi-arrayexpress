
import re
import os
from converter import datamodel
from utils.converter_utils import guess_submission_type_from_study, get_term_from_url, read_json_file
from utils.common_utils import get_ontology_from_term_url


def data_objects_from_json(json_data, json_file_path):
    """
    Read in a USI JSON submission envelope, convert metadata to the common data model and return a submission object

    :param json_data: USI submission envelope with all submittables for an AE experiment submission
    :return: Submission class object
    """

    # We expect one study, which should be guaranteed by the validation, anyway safely getting first element here
    study = datamodel.Study.from_json(next(iter(json_data.get("studies", []))))

    # Now we can fill in the submission info after we get the study object
    submission_info = {
        "team": json_data.get("submission", {}).get("team", {}).get("name"),
        "alias": re.sub("\.json", "", os.path.basename(json_file_path)),
        "metadata": json_file_path,
        "submission_type": guess_submission_type_from_study(study)
    }

    # Samples
    samples = [datamodel.Sample.from_json(s) for s in json_data.get("samples")]

    project = None

    protocols = []
    assays = []
    assay_data = []
    analysis = []

    print(submission_info)

    sub = datamodel.Submission(submission_info, project, study, protocols, samples, assays, assay_data, analysis)

    return sub


class JSONConverter:

    def __init__(self, mapping_file, import_key="ae"):
        self.mapping = read_json_file(mapping_file)
        self.import_key = import_key

        print(mapping_file)

    def convert(self, envelope_json):
        sub_info = {}

        # We only take the first project in the list
        project_json = next(iter(envelope_json.get("projects", [])))
        project = datamodel.Project.from_dict(self.convert_datamodel_object(project_json, "project"))
        print(project)

        # We only take the first study in the list
        study_json = next(iter(envelope_json.get("studies", [])))
        study = datamodel.Study.from_dict(self.convert_datamodel_object(study_json, "study"))

        protocol_json = envelope_json.get("protocols")
        protocols = [datamodel.Protocol.from_dict(self.convert_datamodel_object(p, "protocol")) for p in protocol_json]

        samples = []
        assays = []
        assay_data = []
        analysis = []

        print(project)
        print(study)
        print(protocols)

        submission = datamodel.Submission(sub_info, project, study, protocols, samples, assays, assay_data, analysis)
        return submission

    def convert_datamodel_object(self, submittable_object, submittable_name):
        """
        Use the mapping instructions and convert the different attributes accordingly.
        :param submittable_object: part of a JSON submission envelope containing the details for a single "submittable"
        :param submittable_name: name of the submittable as referenced in the mapping file
        :return: dictionary of the datamodel class attributes with the values to be inserted
        """
        submittable_attributes = {}
        print(submittable_object)
        # Go through attributes in the config
        for attribute, attribute_info in self.mapping.get(submittable_name, {}).items():
            print(attribute)
            # Get information how to convert
            convert_function = None
            mapping_info = attribute_info.get("import", {}).get(self.import_key, {})
            print(mapping_info)
            path = mapping_info.get("path", [])
            method = mapping_info.get("method", "")
            translation = mapping_info.get("translation", {})
            if method:
                convert_function = getattr(self, method)
            if path and convert_function:
                target_object = self.interpret_path(path, submittable_object)
                print(target_object)
                if isinstance(target_object, list):
                    if attribute_info.get("type") == "string":
                        # Take the first entry
                        submittable_attributes[attribute] = next(iter([convert_function(o, translation=translation)
                                                                       for o in target_object]))
                    else:
                        submittable_attributes[attribute] = [convert_function(o, translation=translation)
                                                             for o in target_object]
                elif target_object:
                    submittable_attributes[attribute] = convert_function(target_object, translation=translation)
        print(submittable_attributes)
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
            next_level = json_object.get(path.pop(0), {})
            return self.interpret_path(path, next_level)

    def import_publication(self, element, translation={}):
        return self.convert_datamodel_object(element, "publication")

    def import_contacts(self, element, translation={}):
        return self.convert_datamodel_object(element, "contact")

    @staticmethod
    def generate_attribute_from_json(element, translation={}):

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
            # USI does not support ontology annotations for unit terms
            unit = datamodel.Unit(unit_value, None, None, None)

        attribute = datamodel.Attribute(value,
                                        unit,
                                        term_accession,
                                        term_source)

        return attribute

    @staticmethod
    def get_reference_value_from_json(element, translation={}):
        """References can contain an accession or a combination of alias and team.
        We prefer the accession if there is one, otherwise use the alias (which should be
        unique within the submission)."""
        if element and element.get("accession"):
            return element.get("accession")
        elif element and element.get("alias"):
            return element.get("alias")
        else:
            return ""

    @staticmethod
    def import_string(element, translation={}):
        if translation:
            return translation.get(element, str(element))
        return str(element)



