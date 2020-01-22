
import argparse
import json
import logging
import re
import pkg_resources

from os import path

from converter import json2dm, dm2magetab
from utils.common_utils import file_exists
from utils.converter_utils import read_json_file, dict_to_vertical_table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json',
                        help="Path to the Annotare JSON file")
    parser.add_argument('-o', "--outdir",
                        help="Path where to write the MAGE-TAB files")
    return parser.parse_args()


def guess_submission_type_from_annotare_json(json_data):
    template_name = json_data.get("type")
    if re.search("SINGLE_CELL", template_name):
        return "singlecell"
    elif re.search("SEQUENCING", template_name):
        return "sequencing"
    elif re.search("ARRAY", template_name):
        return "microarray"


def main():
    args = parse_args()

    logger = logging.getLogger()

    # Check if input file exists
    json_file = args.json
    file_exists(json_file)

    json_data = read_json_file(json_file)

    # Determine experiment type
    sub_type = guess_submission_type_from_annotare_json(json_data)

    mapping = json.loads(pkg_resources.resource_string('datamodel',
                                                       path.join("config", "datamodel_mapping_config.json")))

    annotare_converter = json2dm.JSONConverter(mapping, import_key="annotare")

    sub = annotare_converter.convert_annotare_submission(json_data,
                                                         source_file_name=json_file,
                                                         submission_type=sub_type)

    # Generate IDF dictionary

    print(sub.project)
    idf = dm2magetab.generate_idf(sub)
    # Generate SDRF: Output is a pandas dataframe

    #sdrf = dm2magetab.generate_sdrf(sub)

    # New file paths
    prefix = "annotare_test"
    if args.outdir:
        new_idf_file = path.join(args.outdir, prefix + ".idf.txt")
        new_sdrf_file = path.join(args.outdir, prefix + ".sdrf.txt")
    else:
        new_idf_file = path.join(path.dirname(json_file), prefix + ".idf.txt")
        new_sdrf_file = path.join(path.dirname(json_file), prefix + ".sdrf.txt")

    # Write out a new IDF file
    print("Writing {}".format(new_idf_file))
    dict_to_vertical_table(idf, new_idf_file, logger)

    # Rename the columns to the new header list, created by applying a function
    # to "de-uniquify" the header fields, and write new SDRF file
   # dm2magetab.write_sdrf_file(sdrf, new_sdrf_file, logger)


if __name__ == '__main__':
    main()
