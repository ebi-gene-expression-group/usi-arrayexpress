"""A script to define the USI data types and insert them into the database.

Environment variables are needed to specify the following parameters:
DBCON: MongoDB connection string
SUBS_DATABASE: The database name to write to
DATA_TYPE_COLLECTION: The name of the collection for storing the data types
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

    # FG SAMPLE data type on hold, using normal plus checklists
    # "functionalGenomicsSamples": {
    #     "displayNameSingular": "Functional Genomics Sample",
    #     "displayNamePlural": "Functional Genomics Samples",
    #     "description": "",
    #     "validationSchema": load_schema_for_mongo("arrayexpress_sample_schema.json"),
    #     "refRequirements": [{"refClassName": "",
    #                        "dataTypeIdForReferencedDocument": "",
    #                        "additionalRequiredValidationAuthors": []}],
    #     "requiredValidationAuthors": ["core",
    #                                 "JsonSchema",
    #                                 "BioSamples"],
    #     "optionalValidationAuthors": ["Ena",
    #                                 "Taxonomy"],
    #     "submittableClassName": "uk.ac.ebi.subs.repository.model.Sample",
    #     "archive": "BioSamples"},

    "functionalGenomicsProtocols": {
        "displayNameSingular": "Functional Genomics Protocol",
        "displayNamePlural": "Functional Genomics Protocols",
        "description": "Protocols describing the steps of a functional genomics experiment",
        "validationSchema": load_schema_for_mongo("arrayexpress_protocol_schema.json"),
        "refRequirements": [],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Protocol",
        "archive": "ArrayExpress"
    },

    "functionalGenomicsProjects": {
        "displayNameSingular": "Functional Genomics Project",
        "displayNamePlural": "Functional Genomics Projects",
        "description": "Scientific projects including functional genomics studies",
        "validationSchema": load_schema_for_mongo("arrayexpress_project_schema.json"),
        "refRequirements": [],
        "requiredValidationAuthors": [
            "BioStudies",
            "JsonSchema"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Project",
        "archive": "BioStudies"
    },

    "functionalGenomicsMicroarrayStudies": {
        "displayNameSingular": "Functional Genomics Microarray Study",
        "displayNamePlural": "Functional Genomics Microarray Studies",
        "description": "",
        "validationSchema": load_schema_for_mongo("arrayexpress_study_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.ProjectRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsProjects",
                "additionalRequiredValidationAuthors": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.ProtocolRefs",
                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                "additionalRequiredValidationAuthors": []
            }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Study",
        "archive": "ArrayExpress"
    },

    "functionalGenomicsSequencingStudies": {
        "displayNameSingular": "Functional Genomics Sequencing Study",
        "displayNamePlural": "Functional Genomics Sequencing Studies",
        "description": "",
        "validationSchema": load_schema_for_mongo("arrayexpress_study_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.ProjectRef",
                "dataTypeIdForReferencedDocument": "projects",
                "additionalRequiredValidationAuthors": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.ProtocolRefs",
                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                "additionalRequiredValidationAuthors": []
            }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "Ena"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Study",
        "archive": [
            "ArrayExpress",
            "Ena"
        ]
    },

    "functionalGenomicsMicroarrayAssays": {
        "displayNameSingular": "Microarray Assay",
        "displayNamePlural": "Microarray Assays",
        "description": "",
        "validationSchema": load_schema_for_mongo("arrayexpress_microarray_assay_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.StudyRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsMicroarrayStudies",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.ProtocolUse",
                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.SampleUse",
                "dataTypeIdForReferencedDocument": "samples",
                "additionalRequiredValidation": []
            }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Assay",
        "archive": "ArrayExpress"
    },

    "functionalGenomicsSequencingAssays": {
        "displayNameSingular": "Functional Genomics Sequencing Assay",
        "displayNamePlural": "Functional Genomics Sequencing Assays",
        "description": "",
        "validationSchema": load_schema_for_mongo("arrayexpress_sequencing_assay_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.StudyRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsSequencingStudies",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.ProtocolUse",
                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.SampleUse",
                "dataTypeIdForReferencedDocument": "samples",
                "additionalRequiredValidation": []
            }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "Ena"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Assay",
        "archive": [
            "ArrayExpress",
            "Ena"
        ]
    },

    "functionalGenomicsSingleCellAssays": {
        "displayNameSingular": "Single Cell Sequencing Assay",
        "displayNamePlural": "Single Cell Sequencing Assays",
        "description": "",
        "validationSchema": load_schema_for_mongo("arrayexpress_singlecell_assay_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.StudyRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsSequencingStudies",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.ProtocolUse",
                "dataTypeIdForReferencedDocument": "functionalGenomicsProtocols",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.SampleUse",
                "dataTypeIdForReferencedDocument": "samples",
                "additionalRequiredValidation": []
            }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "Ena"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Assay",
        "archive": [
            "ArrayExpress",
            "Ena"
        ]
    },

    "functionalGenomicsMicroarrayAssayData": {
        "displayNameSingular": "Microarray Assay Data",
        "displayNamePlural": "Microarray Assay Data",
        "description": "",
        "validationSchema": load_schema_for_mongo("arrayexpress_assay_data_schema.json"),
        "refRequirements": [
           {
               "refClassName": "uk.ac.ebi.subs.data.component.AssayDataRef",
               "dataTypeIdForReferencedDocument": "functionalGenomicsMicroarrayAssayData",
               "additionalRequiredValidation": []}
           ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "Ena",
            "FileReference"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.AssayData",
        "archive": "ArrayExpress"
    },

    "functionalGenomicsSequencingAssayData": {
        "displayNameSingular": "Functional Genomics Sequencing Assay Data",
        "displayNamePlural": "Functional Genomics Sequencing Assay Data",
        "description": "Data files from one sequencing run or one lane of a functional genomics sequencing experiment",
        "validationSchema": load_schema_for_mongo("arrayexpress_assay_data_schema.json"),
        "refRequirements": [
           {
               "refClassName": "uk.ac.ebi.subs.data.component.AssayDataRef",
               "dataTypeIdForReferencedDocument": "functionalGenomicsSequencingAssayData",
               "additionalRequiredValidation": []}
           ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "Ena"
            "FileReference"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.AssayData",
        "archive": [
            "ArrayExpress",
            "Ena"
        ]
    },

    "functionalGenomicsSingleCellAssayData": {
        "displayNameSingular": "Functional Genomics Single Cell Assay Data",
        "displayNamePlural": "Functional Genomics Single Cell Assay Data",
        "description": "Data files from one sequencing run or one lane of a single cell experiment",
        "validationSchema": load_schema_for_mongo("arrayexpress_assay_data_schema.json"),
        "refRequirements": [
           {
               "refClassName": "uk.ac.ebi.subs.data.component.AssayDataRef",
               "dataTypeIdForReferencedDocument": "functionalGenomicsSingleCellAssayData",
               "additionalRequiredValidation": []}
           ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "Ena"
            "FileReference"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.AssayData",
        "archive": [
            "ArrayExpress",
            "Ena"
        ]
    },
    "functionalGenomicsMicroarrayAnalysis": {
        "displayNameSingular": "Functional Genomics Microarray Analysis",
        "displayNamePlural": "Functional Genomics Microarray Analysis",
        "description": "Processed data files from a functional genomics microarray experiment",
        "validationSchema": load_schema_for_mongo("arrayexpress_analysis_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.AssayDataRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsMicroarrayAssayData",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.AssayRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsMicroarrayAssays",
                "additionalRequiredValidation": []
            }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "FileReference"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Analysis",
        "archive": "ArrayExpress"
    },

    "functionalGenomicsSequencingAnalysis": {
        "displayNameSingular": "Functional Genomics Sequencing Analysis",
        "displayNamePlural": "Functional Genomics Sequencing Analysis",
        "description": "Processed data files from a functional genomics sequencing experiment",
        "validationSchema": load_schema_for_mongo("arrayexpress_analysis_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.AssayDataRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsSequencingAssayData",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.AssayRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsSequencingAssays",
                "additionalRequiredValidation": []
            }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "FileReference"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Analysis",
        "archive": "ArrayExpress"
    },

    "functionalGenomicsSingleCellAnalysis": {
        "displayNameSingular": "Functional Genomics Single Cell Analysis",
        "displayNamePlural": "Functional Genomics Single Cell Analysis",
        "description": "Processed data files from a single cell experiment",
        "validationSchema": load_schema_for_mongo("arrayexpress_analysis_schema.json"),
        "refRequirements": [
            {
                "refClassName": "uk.ac.ebi.subs.data.component.AssayDataRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsSingleCellAssayData",
                "additionalRequiredValidation": []
            },
            {
                "refClassName": "uk.ac.ebi.subs.data.component.AssayRef",
                "dataTypeIdForReferencedDocument": "functionalGenomicsSingleCellAssays",
                "additionalRequiredValidation": []
           }
        ],
        "requiredValidationAuthors": [
            "core",
            "JsonSchema",
            "FileReference"
        ],
        "optionalValidationAuthors": [],
        "submittableClassName": "uk.ac.ebi.subs.repository.model.Analysis",
        "archive": "ArrayExpress"
        },
    }


# Get MongoDB connection string from environment variable $DBCON
db_url = environ.get("DBCON")

# Connecting to MongoDB
client = MongoClient(db_url)
db = getattr(client, environ.get("SUBS_DATABASE"))

# For testing the connection
files_entry = getattr(db, environ.get("DATA_TYPE_COLLECTION")).find_one({"_id": "files"})
pprint(files_entry)


if __name__ == "__main__":
    #print(load_schema_for_mongo("arrayexpress_sample_schema.json"))
    #print(re.sub(r"\$", "#dollar#", "{$schema: test string}"))
    pass
