"""Module to convert experiment metadata in MAGE-TAB (IDF/SDRF) format to USI submittable JSON format"""


from converter.parsing import parse_idf, parse_sdrf
from converter.datamodel import *
from converter.converting import *

from utils.common_utils import create_logger
from utils.converter_utils import write_json_file, strip_extension


SDRF_FILE_NAME_REGEX = r"^\s*SDRF\s*File"
DATA_DIRECTORY = "unpacked"


def mtab2usi_conversion(idf_file_path):
    process_name = "mtab2usi_conversion"

    current_dir = os.path.dirname(idf_file_path)
    idf_file_name = os.path.basename(idf_file_path)
    sdrf_file_path = None

    # Create logger
    logger = create_logger(current_dir, process_name, idf_file_name)

    # Figure out the name and location of sdrf files
    with codecs.open(idf_file_path, 'rU', encoding='utf-8') as f:
        # U flag makes it portable across in unix and windows (\n and \r\n are treated the same)
        for line in f:
            if re.search(SDRF_FILE_NAME_REGEX, line):
                sdrf_file_name = line.split("\t")[1].strip()
                if os.path.exists(current_dir + DATA_DIRECTORY):
                    sdrf_file_path = os.path.join(current_dir, DATA_DIRECTORY, sdrf_file_name)
                else:
                    sdrf_file_path = os.path.join(current_dir, sdrf_file_name)

    logger.debug("Found SDRF file: {}".format(sdrf_file_path))

    sub = data_objects_from_magetab(idf_file_path, sdrf_file_path)

    datamodel2json_conversion(sub, current_dir, logger)


def data_objects_from_magetab(idf_file_path, sdrf_file_path):
    """
    Parse IDF/SDRF files and transform metadata to common datamodel

    :param idf_file_path: string, path to IDF file
    :param sdrf_file_path: string, path to SDRF file
    :return: Submission class object
    """

    study_info, protocols = parse_idf(idf_file_path)
    samples, extracts, le, assays, raw_data, processed_data, is_microarray = parse_sdrf(sdrf_file_path)

    # For MAGE-TAB files we don't have USI submission info might need to store these somewhere once we get this
    idf_file_name = os.path.basename(idf_file_path)
    sub_info = {"alias": re.sub("\.idf\.txt$", "", idf_file_name),
                "accession": study_info.get("accession"),
                "team": "my-super-test-team",
                "metadata": idf_file_path}

    # Project
    project_object = Project.from_magetab(study_info)

    # Study
    study_object = Study.from_magetab(study_info)

    # Protocols
    protocol_objects = []
    for p in protocols:
        protocol = Protocol.from_magetab(p)
        protocol_objects.append(protocol)

    # Samples
    sample_objects = []
    for sample in samples.values():
        new_sample = Sample.from_magetab(sample)
        sample_objects.append(new_sample)

    # Assays
    assay_objects = []

    if is_microarray:
        linked_extracts = []
        for le_name, le_attributes in le.items():
            for extract_name, extract_attributes in extracts.items():
                if extract_name in le_attributes["extract_ref"]:
                    # Assuming there is only one extract per le
                    linked_extracts = extract_attributes
                    break

            # Get all assays referencing this extract
            linked_assays = []
            for assay_name, assay_attributes in assays.items():

                if le_name in assay_attributes["extract_ref"]:
                    linked_assays.append(assay_attributes)

            new_assay = MicroarrayAssay.from_magetab(le_attributes, linked_extracts, linked_assays)
            assay_objects.append(new_assay)

    else:
        for extract_name, extract_attributes in extracts.items():

            # Get all assays referencing this extract
            linked_assays = []
            for assay_name, assay_attributes in assays.items():
                if extract_name in assay_attributes["extract_ref"]:
                    linked_assays.append(assay_attributes)

            print(extract_name, len(linked_assays))

            new_assay = SeqAssay.from_magetab(extract_attributes, linked_assays, protocols)
            assay_objects.append(new_assay)

    # Assay data
    ad_objects = []

    file_groups = OrderedDict()
    for f_name, f_attrib in raw_data.items():
        if len(f_attrib.get("assay_ref")) > 1:
            # For matrix files, one object per file
            file_groups[strip_extension(f_name)] = [f_attrib]

        elif len(f_attrib.get("assay_ref")) == 1:
            a_ref = f_attrib.get("assay_ref")[0]
            # Get the other files with the same assay ref
            file_groups[a_ref] = [f_attrib for f_attrib in raw_data.values() if a_ref in f_attrib.get("assay_ref")]

    for name, group in file_groups.items():
        file_objects = [DataFile.from_magetab(f_attrib) for f_attrib in group]
        assay_data = AssayData.from_magetab(name, file_objects, group)
        ad_objects.append(assay_data)

    sub = Submission(sub_info,
                     project_object,
                     study_object,
                     protocol_objects,
                     sample_objects,
                     assay_objects,
                     ad_objects,
                     [])

    return sub


def datamodel2json_conversion(submission, working_dir, logger):
    """
    Take metadata in common datamodel and write JSON files

    :param submission: object, Submission class object that holds metadata of the whole experiment
    :param working_dir: string, directory to write files to
    :param logger: object, log handler
    :return: None
    """
    
    # Dict to store USI objects to write to file:
    json_objects = {
        "project": generate_usi_project_object(submission.project),
        "study": generate_usi_study_object(submission.study, submission.info),
        "protocol": [generate_usi_protocol_object(p) for p in submission.protocol],
        "sample": [generate_usi_sample_object(s) for s in submission.sample],
        "assay": [generate_usi_assay_object(a, submission.info) for a in submission.assay],
        "assay_data": [generate_usi_data_object(ad, submission.info) for ad in submission.assay_data]
    }

    # Write individual JSON files
    for submittable_type, objects in json_objects.items():
        logger.info("Writing JSON file for {}.".format(submittable_type))
        write_json_file(working_dir, objects, submittable_type, submission.info)

