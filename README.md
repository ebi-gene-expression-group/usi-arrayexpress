# usi-arrayexpress

This repository contains modules for transformation of experimental metadata between [MAGE-TAB](http://fged.org/projects/mage-tab/) and JSON format used in the [Unified Submissions Interface (USI)](https://github.com/EMBL-EBI-SUBS). 


## Installation and requirements

Python 3.6 interpreter<br>
Package requirements:
* requests 2.20.1
* jsonschema 2.6.0
* pandas 0.24.2

Add usi-arrayexpress directory to PYTHONPATH environment variable


## Converter

The converter functionality can read in MAGE-TAB files and generate a set of USI JSON files. This workflow can be tested with the mtab2usi_conversion.py script, e.g. 
```
python mtab2usi_conversion.py tests/test_data/E-MTAB-4250.idf.txt
```
This will read in the IDF file and the SDRF file that is specified in the IDF. It transforms the metadata into a Python class data model that is roughly based on the [USI submissions data model](https://github.com/EMBL-EBI-SUBS/subs-data-model).<br>
 The output JSON files are created in a sub-folder, in the location of the IDF file. The JSON structure is based on the [USI JSON schemas](https://github.com/EMBL-EBI-SUBS/validation-schemas) modified to accommodate ArrayExpress specific metadata fields. 
 
 
 ### MAGE-TAB writer
 
 The datamodel2magetab converter module can take data stored in the common data model and write it out as MAGE-TAB files. Study, project and protocols metadata get combined in the IDF table, while sample, assay and file metadata are combined in the SDRF table. The test script `run_magetab_writer.py` takes MAGE-TAB files as input, converts the data to the data model and outputs new IDF/SDRF files. 
 
 
 ## JSON schemas
 
 Prototypes for the JSON schema describing the metadata required for ArrayExpress submissions.<br>
 JSON schema validation using [jsonschema package](https://github.com/Julian/jsonschema) can be run with `json_validation.py` script, e.g.
 ```
python json_validation.py tests/test_data/simple_data.json -s tests/test_data/simple_schema.json
 ```
 
 
 ## Validator
 
 The validator submits the metadata that has been converted into the Python class data model to several checks that make sure all mandatory information has been included. 
 The types of checks are roughly based on the checks in the [MAGE-TAB checker](https://github.com/arrayexpress/magetabcheck) used by [Annotare](https://github.com/arrayexpress/annotare2), with a few additions so far. <br>
 A pre-validation module is used to check IDF and SDRF files before reading them into the datamodel. These checks include validation of MAGE-TAB fields.<br>
 The validation can be run on a given set of MAGE-TAB files using the magetab_validation.py script, e.g. 
 ```
 python magetab_validation.py tests/test_data/E-MTAB-4250.idf.txt
 ```
  
  ### Single-cell Expression Atlas MAGE-TAB validator
  
  A separate MAGE-TAB pre-validation module is running checks that guarantee that the experiment can be processed for [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home). It reads metadata directly from the MAGE-TAB (not from the data model) since the specific format of the IDF/SDRF are important for correct analysis. 
  The checks can be invoked using the MAGE-TAB validation script
  ```
  python magetab_validation.py tests/test_data/E-MTAB-7704.idf.txt -a 
  ```
 
  