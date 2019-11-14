"""A script to define the USI data types insert them into the database.

Environment variables are needed to specify the following parameters:
DBCON: MongoDB connection string
SUBS_DATABASE: The database name to write to
DATA_TYPE_COLLECTION: The name of the collection for storing the datatypes
"""


import json
import re
from os import environ
from pymongo import MongoClient
from pprint import pprint


def load_schema_for_mongo(schema_file):
    """The version of MongoDB that is used doesn't support special characters like $ and . at the start
     of a keyword. These need to be replaced before inserting the JSON schema into the database."""
    try:
        with open(schema_file) as sf:
            json_string = json.dumps(json.load(sf))
            mongo_string = re.sub(r"\"\$", "\"#dollar#", re.sub(r"\"\.", "\"#dot#", json_string))
            return json.loads(mongo_string)
    except IOError:
        print("Schema file {} not found or not readable.".format(schema_file))


new_data_types = {
    "functionalGenomicsSamples": {"displayNameSingular": "Functional Genomics Sample",
                                  "displayNamePlural": "Functional Genomics Samples",
                                  "description": "",
                                  "validationSchema": load_schema_for_mongo("arrayexpress_sample_schema.json"),
                                  "refRequirements": [{"refClassName": "",
                                                       "dataTypeIdForReferencedDocument": "",
                                                       "additionalRequiredValidationAuthors": []}],
                                  "requiredValidationAuthors": ["core",
                                                                "JsonSchema",
                                                                "BioSamples"],
                                  "optionalValidationAuthors": ["Ena",
                                                                "Taxonomy"],
                                  "submittableClassName": "uk.ac.ebi.subs.repository.model.Sample",
                                  "archive": "BioSamples"},

    "functionalGenomicsStudies": {"displayNameSingular": "Functional Genomics Study",
                                  "displayNamePlural": "Functional Genomics Studies",
                                  "description": "",
                                  "validationSchema": load_schema_for_mongo("arrayexpress_study_schema.json"),
                                  "refRequirements": [{"refClassName": "uk.ac.ebi.subs.data.component.ProjectRef",
                                                       "dataTypeIdForReferencedDocument": "projects",
                                                       "additionalRequiredValidationAuthors": []},
                                                      {"refClassName": "uk.ac.ebi.subs.data.component.ProtocolRefs",
                                                       "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                                                       "additionalRequiredValidationAuthors": []}
                                                      ],
                                  "requiredValidationAuthors": ["core",
                                                                "JsonSchema"],
                                  "optionalValidationAuthors": [],
                                  "submittableClassName": "uk.ac.ebi.subs.repository.model.Study",
                                  "archive": "ArrayExpress"},

    "functionalGenomicsProtocols": {"displayNameSingular": "Functional Genomics Protocol",
                                    "displayNamePlural": "Functional Genomics Protocols",
                                    "description": "Protocols describing the steps of a functional genomics experiment",
                                    "validationSchema": load_schema_for_mongo("arrayexpress_protocol_schema.json"),
                                    "refRequirements": [],
                                    "requiredValidationAuthors": ["core",
                                                                  "JsonSchema"],
                                    "optionalValidationAuthors": [],
                                    "submittableClassName": "uk.ac.ebi.subs.repository.model.Protocol",
                                    "archive": "ArrayExpress"},

    "functionalGenomicsMicroarrayAssays": {"displayNameSingular": "Microarray Assay",
                                           "displayNamePlural": "Microarray Assays",
                                           "description": "",
                                           "validationSchema": load_schema_for_mongo("arrayexpress_microarray_assay_schema.json"),
                                           "refRequirements": [
                                               {"refClassName": "uk.ac.ebi.subs.data.component.StudyRef",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsStudies",
                                                "additionalRequiredValidation": []},
                                               {"refClassName": "uk.ac.ebi.subs.data.component.ProtocolUse",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                                                "additionalRequiredValidation": []},
                                               {"refClassName": "uk.ac.ebi.subs.data.component.SampleUse",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsSamples",
                                                "additionalRequiredValidation": []}],
                                           "requiredValidationAuthors": ["core",
                                                                         "JsonSchema"],
                                           "optionalValidationAuthors": [],
                                           "submittableClassName": "uk.ac.ebi.subs.repository.model.Assay",
                                           "archive": "ArrayExpress"},

    "functionalGenomicsSequencingAssays": {"displayNameSingular": "Functional Genomics Sequencing Assay",
                                           "displayNamePlural": "Functional Genomics Sequencing Assays",
                                           "description": "",
                                           "validationSchema": load_schema_for_mongo("arrayexpress_sequencing_assay_schema.json"),
                                           "refRequirements": [
                                               {"refClassName": "uk.ac.ebi.subs.data.component.StudyRef",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsStudies",
                                                "additionalRequiredValidation": []},
                                               {"refClassName": "uk.ac.ebi.subs.data.component.ProtocolUse",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                                                "additionalRequiredValidation": []},
                                               {"refClassName": "uk.ac.ebi.subs.data.component.SampleUse",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsSamples",
                                                "additionalRequiredValidation": []}],
                                           "requiredValidationAuthors": ["core",
                                                                         "JsonSchema"],
                                           "optionalValidationAuthors": [],
                                           "submittableClassName": "uk.ac.ebi.subs.repository.model.Assay",
                                           "archive": "ArrayExpress"},

    "functionalGenomicsSingleCellAssays": {"displayNameSingular": "Single Cell Sequencing Assay",
                                           "displayNamePlural": "Single Cell Sequencing Assays",
                                           "description": "",
                                           "validationSchema": load_schema_for_mongo("arrayexpress_singlecell_assay_schema.json"),
                                           "refRequirements": [
                                               {"refClassName": "uk.ac.ebi.subs.data.component.StudyRef",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsStudies",
                                                "additionalRequiredValidation": []},
                                               {"refClassName": "uk.ac.ebi.subs.data.component.ProtocolUse",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                                                "additionalRequiredValidation": []},
                                               {"refClassName": "uk.ac.ebi.subs.data.component.SampleUse",
                                                "dataTypeIdForReferencedDocument": "functionalGenomicsSamples",
                                                "additionalRequiredValidation": []}],
                                           "requiredValidationAuthors": ["core",
                                                                         "JsonSchema"],
                                           "optionalValidationAuthors": [],
                                           "submittableClassName": "uk.ac.ebi.subs.repository.model.Assay",
                                           "archive": "ArrayExpress"},

    "functionalGenomicsAssayData": {"displayNameSingular": "Single Cell Sequencing Assay",
                                    "displayNamePlural": "Single Cell Sequencing Assays",
                                    "description": "",
                                    "validationSchema": load_schema_for_mongo("arrayexpress_assay_data_schema.json"),
                                    "refRequirements": [
                                       {"refClassName": "uk.ac.ebi.subs.data.component.AssayRef",
                                        "dataTypeIdForReferencedDocument": "functionalGenomicsStudies",
                                        "additionalRequiredValidation": []}
                                       ],
                                    "requiredValidationAuthors": ["core",
                                                                  "JsonSchema",
                                                                  "FileReference"],
                                    "optionalValidationAuthors": [],
                                    "submittableClassName": "uk.ac.ebi.subs.repository.model.AssayData",
                                    "archive": "ArrayExpress"},
    "functionalGenomicsAnalysis": {},
    "functionalGenomicsFiles": {}
    }


# Get MongoDB connection string from environment variable $DBCON
db_url = environ.get("DBCON")


# connect to MongoDB, change the << MONGODB URL >> to reflect your own connection string
client = MongoClient(db_url)
db = getattr(client, environ.get("SUBS_DATABASE"))


files_entry = getattr(db, environ.get("DATA_TYPE_COLLECTION")).find_one({"_id": "files"})

pprint(files_entry)




if __name__ == "__main__":
    #print(load_schema_for_mongo("arrayexpress_sample_schema.json"))
    #print(re.sub(r"\$", "#dollar#", "{$schema: test string}"))
    pass
