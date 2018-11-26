

from converter.datamodel import Assay, MicroarrayAssay
from utils.converter_utils import read_json_file


json_data = read_json_file("test_data/E-MTAB-4250.usi.json")

print(json_data)
assay_data = json_data.get("assays")

for x in assay_data:
    my_assay = Assay(x)
    print(my_assay.name)
    print(my_assay.technology_type)
    print(my_assay.get_technology_type())

