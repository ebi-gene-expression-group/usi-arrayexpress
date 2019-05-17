# usi-arrayexpress

This repository contains modules for transformation of experimental metadata between [MAGE-TAB](fged.org/projects/mage-tab/) and JSON format used in the [Unified Submissions Interface (USI)](https://github.com/EMBL-EBI-SUBS). 


## Installation and requirements

Python 3.6 interpreter<br>
Package requirements:git 
* requests 2.20.1
* jsonschema 2

Add usi-arrayexpress directory to PYTHONPATH environment variable


## Converter

The converter functionality so far can read in MAGE-TAB files and generate a set of USI JSON files. This workflow can be tested with the mtab2usi_conversion.py script, e.g. `python mtab2usi_conversion tests/test_data/E-MTAB-4250.idf.txt`<br>
This will read in the IDF file and the SDRF file that is specified in the IDF. It transforms the metadata into a Python class data model that is roughly based on the [USI submissions data model](https://github.com/EMBL-EBI-SUBS/subs-data-model).<br>
 The output JSON files are created in a sub-folder, in the location of the IDF file. The JSON structure is based on the [USI JSON schemas](https://github.com/EMBL-EBI-SUBS/validation-schemas) modified to accommodate ArrayExpress specific metadata fields. 
 
 
 ## JSON schemas
 
 Prototypes for the JSON schemas describing the metadata required for ArrayExpress submissions. 
 
 
 ## Validator
 
 The validator submits the metadata that has been converted into the Python class data model to several checks that make sure all mandatory information has been included. 
 The types of checks are roughly based on the checks in the [MAGE-TAB checker](https://github.com/arrayexpress/magetabcheck) used by [Annotare](https://github.com/arrayexpress/annotare2), with a few additions so far. <br>
 A pre-validation module is used to check IDF and SDRF files before reading them into the datamodel. These checks include validation of MAGE-TAB fields.<br>
 The validation can be run on a given set of MAGE-TAB files using the magetab_validation.py script, e.g. `python magetab_validation.py tests/test_data/E-MTAB-4250.idf.txt`
 
