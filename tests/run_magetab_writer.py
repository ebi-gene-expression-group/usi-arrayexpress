# for testing
import utils.converter_utils
from utils.common_utils import create_logger, file_exists
from utils.converter_utils import get_sdrf_path, guess_submission_type
from converter.converting import data_objects_from_magetab
from converter import datamodel2magetab
from sys import argv
import os
import pandas


def column_name_to_magetab(header):
    """Transform unique column header back to MAGE-TAB style.
    Column headers are expected to be separated by "~~~". """
    header_list = header.split("~~~")
    return header_list[-1]


process_name = "magetab_writer"
idf_file = argv[1]
file_exists(idf_file)

# Create logger
current_dir, idf_file_name = os.path.split(idf_file)
logger = create_logger(current_dir, process_name, idf_file_name)

# Get path to SDRF file
sdrf_file_path = get_sdrf_path(idf_file, logger)
print(sdrf_file_path)

# Get correct submission type
sub_type, idf_dict = guess_submission_type(idf_file, sdrf_file_path, logger)

# Read in MAGE-TAB and convert to common data model
sub = data_objects_from_magetab(idf_file, sdrf_file_path, sub_type)

idf = datamodel2magetab.generate_idf(sub)

utils.converter_utils.dict_to_vertcial_table(idf, idf_file + "_new.txt")


# Output is a pandas dataframe
raw_out = datamodel2magetab.generate_sdrf(sub)


#utils.converter_utils.tuple_list_to_table(raw_out, sdrf_file_path + "_new.txt")

# Rename the columns to the new header list, created by applying a function to "de-uniquify" the header fields
datamodel2magetab.write_sdrf_file(raw_out, sdrf_file_path + "_new.txt", logger)



