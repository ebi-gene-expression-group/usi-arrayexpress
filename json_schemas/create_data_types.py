
import json
import re
from os import environ
from pymongo import MongoClient
from pprint import pprint


def load_schema_for_mongo(schema_file):
    try:
        with open(schema_file) as sf:
            json_string = json.dumps(json.load(sf))
            print(json_string)
            new_string = re.sub(r"\$", "#dollar#", json_string)
            return json.loads(new_string)
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
    "functionalGenomicsProtocols": {},
    "functionalGenomicsMicroarrayAssays": {},
    "functionalGenomicsSequencingAssays": {},
    "functionalGenomicsSingleCellAssays": {},
    "functionalGenomicsAssayData": {},
    "functionalGenomicsAnalysis": {},
    "functionalGenomicsFiles": {}
    }


# Get MongoDB connection string from environment variable $DBCON
db_url = environ.get("DBCON")


# connect to MongoDB, change the << MONGODB URL >> to reflect your own connection string
client = MongoClient(db_url)
db = client.submissions_dev


files_entry = db.dataType.find({"_id": "files"})

pprint(files_entry)
# Issue the serverStatus command and print the results
#serverStatusResult = db.command("serverStatus")
#pprint(serverStatusResult)






if __name__ == "__main__":
    print(load_schema_for_mongo("arrayexpress_sample_schema.json"))
    print(re.sub(r"\$", "#dollar#", "{$schema: test string}"))
