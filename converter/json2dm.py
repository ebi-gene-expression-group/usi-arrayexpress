
import re
import os
from converter import datamodel
from utils.converter_utils import guess_submission_type_from_json


def data_objects_from_json(json_data, json_file_path):
    """
    Read in a USI JSON submission envelope, convert metadata to the common data model and return a submission object

    :param json_data: USI submission envelope with all submittables for an AE experiment submission
    :return: Submission class object
    """

    samples = [datamodel.Sample.from_json(s) for s in json_data.get("samples")]

    project = None
    study = None
    protocols = []
    assays = []
    assay_data = []
    analysis = []

    submission_info = {
        "team": json_data.get("submission", {}).get("team", {}).get("name"),
        "alias": re.sub("\.json", "", os.path.basename(json_file_path)),
        "metadata": json_file_path,
        "submission_type": guess_submission_type_from_json(study)
    }
    print(submission_info)

    sub = datamodel.Submission(submission_info, project, study, protocols, samples, assays, assay_data, analysis)

    return sub
