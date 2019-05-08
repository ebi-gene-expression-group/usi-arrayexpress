
from converter.converting import *
from converter.datamodel import SeqAssay, Protocol, MicroarrayAssay, DataFile, AssayData, Study, Project, Sample
from converter.parsing import parse_sdrf, parse_idf
from utils.converter_utils import strip_extension, write_json_file

wd = dir_path = os.path.dirname(os.path.realpath(__file__))
idf_path = os.path.join(wd, 'test_data', 'E-MTAB-6160.idf.txt')
sdrf_path = os.path.join(wd, 'test_data', 'E-MTAB-6160.sdrf.txt')
print(idf_path)
study_info2, protocols = parse_idf(idf_path)
samples, extracts, le, assays, raw_data, processed_data, is_microarray = parse_sdrf(sdrf_path)


# Test submission info dict
sub_info = {
    "team": "my-super-test-team",
    "alias": "submission1234",
    "accession": "E-MTAB-9606"
}


# Project
project = Project.from_magetab(study_info2)
print(project)
project_object = generate_usi_project_object(project)
write_json_file(project_object, "project", sub_info)

# Study
study = Study.from_magetab(study_info2)
print(study)
print(study.experimental_factor)
study_object = generate_usi_study_object(study, sub_info)
write_json_file(study_object, "study", sub_info)

# Protocols
protocol_objects = []
for p in protocols:
    protocol = Protocol.from_magetab(p)
    # This doesn't run, instead returns None!?
    protocol_object = generate_usi_protocol_object(protocol)
    protocol_objects.append(protocol_object)
write_json_file(protocol_objects, "protocol", sub_info)




print("Found Samples: {}".format(len(samples)))
print("Found Extracts: {}".format(len(extracts)))
print("Found Lableled Extracts: {}".format(len(le)))
print("Found Assays: {}".format(len(assays)))


# Samples
sample_objects = []
for sample in samples.values():
    new_sample = Sample.from_magetab(sample)
    print(new_sample.alias, new_sample.taxon, new_sample.taxonId)
    print(new_sample.attributes)
    sample_object = generate_usi_sample_object(new_sample)
    sample_objects.append(sample_object)
write_json_file(sample_objects, "sample", sub_info)


# Assays
assay_objects = []

if is_microarray:
    linked_extracts = []
    linked_assays = []
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
                print(assay_attributes)
                linked_assays.append(assay_attributes)

        print(le_name, linked_extracts.get("name"), len(linked_assays))

        new_assay = MicroarrayAssay.from_magetab(le_attributes, linked_extracts, linked_assays)

        print(new_assay.alias)
        print(new_assay.protocolrefs)
        print(new_assay.get_attributes())

        assay_object = generate_usi_assay_object(new_assay, sub_info)
        assay_objects.append(assay_object)
    print(new_assay)

    # Raw data
    for file_name, file_attributes in raw_data.items():
        new_file = DataFile.from_magetab(file_attributes)
        print(new_file)

else:
    for extract_name, extract_attributes in extracts.items():

        # Get all assays referencing this extract
        linked_assays = []
        for assay_name, assay_attributes in assays.items():
            if extract_name in assay_attributes["extract_ref"]:
                linked_assays.append(assay_attributes)

        print(extract_name, len(linked_assays))

        new_assay = SeqAssay.from_magetab(extract_attributes, linked_assays, protocols)

        print(new_assay.alias)
        print(new_assay.accession)
        print(new_assay.protocolrefs)
        print(new_assay.get_attributes())

        assay_object = generate_usi_assay_object(new_assay, sub_info)
        assay_objects.append(assay_object)


write_json_file(assay_objects, "assay", sub_info)


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
    ad_objects.append(generate_usi_data_object(assay_data, sub_info))

print(ad_objects)
write_json_file(ad_objects, "assay_data", sub_info)
