# usi-arrayexpress

This repository contains modules for transformation of experimental metadata between [MAGE-TAB](fged.org/projects/mage-tab/) and JSON format used in the [Unified Submissions Interface (USI)](https://github.com/EMBL-EBI-SUBS). 


## Installation and requirements

Python 3.6 interpreter<br>
Package requirements: requests 2.20.1 <br>
Add usi-arrayexpress directory to PYTHONPATH environment variable


## Converter

The converter functionality so far can read in MAGE-TAB files and generate a set of USI JSON files. This workflow can be tested with the mtab2usi_conversion.py script, e.g. `python mtab2usi_conversion tests/test_data/E-MTAB-4250.idf.txt`<br>
This will read in the IDF file and the SDRF file that is specified in the IDF. It transforms the metadata into a Python class data model that is roughly based on the [USI submissions data model](https://github.com/EMBL-EBI-SUBS/subs-data-model).<br>
 The output JSON files are created in a sub-folder, in the location of the IDF file. The JSON structure is based on the [USI JSON schemas](https://github.com/EMBL-EBI-SUBS/validation-schemas) modified to accommodate ArrayExpress specific metadata fields. 
 
 ## JSON schemas
 
 Prototypes for the JSON schemas describing the metadata required for ArrayExpress submissions. 