
import re
import os
from converter import datamodel
from utils.converter_utils import guess_submission_type_from_study, get_term_from_url
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


def generate_attribute_from_json(element):

    term = next(iter(element.get("terms", [])))
    term_accession = get_term_from_url(term.get("url"))
    term_source = get_ontology_from_term_url(term.get("url"))

    unit = datamodel.Unit(element.get("units"), None, None, None)

    attribute = datamodel.Attribute(element.get("value"),
                                    unit,
                                    term_accession,
                                    term_source)

    return attribute
