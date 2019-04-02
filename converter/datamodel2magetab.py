"""Module to convert metadata in the submission data model to MAGE-TAB files."""

from collections import OrderedDict
from utils import converter_utils
import codecs
import csv
import pandas as pd
import itertools



def generate_idf(sub):
    idf = OrderedDict([
        ("MAGE-TAB Version", "1.1"),
        ("Investigation Title", sub.study.title),
        ("Experiment Description", sub.study.description),
        ("Experimental Design", [d.value for d in sub.study.experimental_design]),
        ("Experimental Design Term Source REF", [d.term_source for d in sub.study.experimental_design]),
        ("Experimental Design Term Accession Number", [d.term_accession for d in sub.study.experimental_design]),
        ("Experimental Factor Name", [f.value for f in sub.study.experimental_factor]),
        ("Experimental Factor Type", [f.value for f in sub.study.experimental_factor]),
        ("Experimental Factor Term Source REF", [f.term_source for f in sub.study.experimental_factor]),
        ("Experimental Factor Term Accession Number", [f.term_accession for f in sub.study.experimental_factor]),
        ("Person Last Name", [p.lastName for p in sub.project.contacts]),
        ("Person First Name", [p.firstName for p in sub.project.contacts]),
        ("Person Mid Initials", [p.middleInitials for p in sub.project.contacts]),
        ("Person Email", [p.email for p in sub.project.contacts]),
        ("Person Phone", [p.phone for p in sub.project.contacts]),
        ("Person Fax", [p.fax for p in sub.project.contacts]),
        ("Person Address", [p.address for p in sub.project.contacts]),
        ("Person Affiliation", [p.affiliation for p in sub.project.contacts]),
        ("Person Roles", [';'.join(p.roles) for p in sub.project.contacts]),
        ("Date of Experiment", sub.study.date_of_experiment),
        ("Public Release Date", sub.project.releaseDate),
        ("PubMed ID", [p.pubmedId for p in sub.project.publications]),
        ("Publication DOI", [p.doi for p in sub.project.publications]),
        ("Publication Author List", [p.authors for p in sub.project.publications]),
        ("Publication Title", [p.articleTitle for p in sub.project.publications]),
        ("Publication Status", [p.status for p in sub.project.publications]),
        ("Publication Status Term Source REF", ["EFO" for p in sub.project.publications if p.status]),
        ("Publication Status Term Accession Number", []),
        ("Protocol Name", [p.alias for p in sub.protocol]),
        ("Protocol Type", [p.protocol_type.value for p in sub.protocol]),
        ("Protocol Term Source REF", [p.protocol_type.term_source for p in sub.protocol]),
        ("Protocol Term Accession Number", [p.protocol_type.term_accession for p in sub.protocol]),
        ("Protocol Description", [p.description for p in sub.protocol]),
        ("Protocol Hardware", [p.hardware for p in sub.protocol]),
        ("Protocol Software", [p.software for p in sub.protocol]),
        ("SDRF File", sub.info.get("accession") + ".sdrf.txt"),
        ("Term Source Name", ["EFO", "ArrayExpress"]),
        ("Term Source File", ["https://www.ebi.ac.uk/efo/", "https://www.ebi.ac.uk/arrayexpress/"]),
        ("Comment[AEExperimentType", [exptype for exptype in sub.study.experiment_type]),
        ("Comment[ArrayExpressAccession]", sub.info.get("accession"))
    ])

    if sub.info.get("submission_type") == "sequencing" or sub.info.submission_type == "singlecell":
        pass
        # TODO: Add comments for ENA accessions Comment[SecondaryAccession] and Comment[SequenceDataURI]

    return idf


def generate_sdrf(sub):
    """Transform metadata in data model to an SDRF table."""

    header = []
    rows = []
    all_sample_categories = []

    if sub.info.get("submission_type") == "microarray":

        for sample in sub.sample:

            sample_values = [
                ("Source Name", sample.alias),
                #("Characteristics", sample.attributes),
                ("Description", sample.description),
                ("Material Type", sample.material_type)
            ]
            # Expand sample attributes to characteristics columns (they can be different between different samples)
            for category, sample_attrib in sample.attributes.items():
                sample_values.extend(flatten_sample_attribute(category, sample_attrib, make_unique=True))

            # Get all assay objects that belong to this sample
            assays = [assay for assay in sub.assay if assay.sampleref == sample.alias]

            for assay in assays:

                assay_values = [
                    ("Extract Name", sample.alias),
                    ("labeling~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                      if sub.get_protocol(pref).protocol_type.value == "nucleic acid labeling protocol"]),
                    ("Labelled Extract Name", assay.alias),
                    ("hyb~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                      if sub.get_protocol(pref).protocol_type.value == "nucleic acid hybridization to array protocol"])
                ]

                # Get all assay data objects that belong to this assay
                data = [ad for ad in sub.assay_data if assay.alias in ad.assayrefs]

                for ad in data:
                    data_values = [
                        ("Assay Name", ad.alias),
                        ("Technology Type", assay.technology_type),
                        ("Array Design REF", assay.array_design),
                        ("design~Term Source REF", "Array Express"),
                        ("scanning~Protocol REF", [sub.get_protocol(pref).alias for pref in ad.protocolrefs
                                      if sub.get_protocol(pref).protocol_type.value == "array scanning and feature extraction protocol"]),
                    ]

                    # Get all data files
                    for f in ad.files:
                        if ad.data_type == "raw":
                            data_values.append(("Array Data File", f.name))
                        elif ad.data_type == "raw matrix":
                            data_values.append(("Array Data Matrix File", f.name))

                    # Add all lists together and transform to dictionary so that we can use pandas to write to table
                    rows.append((OrderedDict(sample_values), OrderedDict(assay_values), OrderedDict(data_values)))



    # TODO: Expand protocol refs as they can be variable length


    # Get all keys (categories) from the dictionaries/rows (in row-list position 1, tuple position 1)

    print(rows)
    sample_dicts = [row[0] for row in rows]
    data_frames = []
    for i in range(len(rows[0])):

        data_frames.append(pd.DataFrame.from_records([row[i] for row in rows ]))

    print(pd.concat(data_frames, axis=1))

    return pd.concat(data_frames, axis=1)

    #return rows


def flatten_unit(category, unit_object, make_unique=False, sep="~~~"):
    flat_list = []
    if unit_object.value and unit_object.unit_type:
        if make_unique:
            unit_header = category + sep + "Unit[{}]".format(unit_object.unit_type)
        else:
            unit_header = "Unit[{}]".format(unit_object.unit_type)
        flat_list.append((unit_header, unit_object.value))
        if unit_object.term_source:
            if make_unique:
                flat_list.append((category + "-unit" + sep + "Term Source REF", unit_object.term_source))
            else:
                flat_list.append(("Term Source REF", unit_object.term_source))
            if unit_object.term_accession:
                if make_unique:
                    flat_list.append((category + "-unit" + sep + "Term Accession Number", unit_object.term_accession))
                else:
                    flat_list.append(("Term Accession Number", unit_object.term_accession))
    return flat_list


def flatten_sample_attribute(category, attrib_object, make_unique=False, sep="~~~"):
    flat_list = []
    if attrib_object.value:
        header = "Characteristics[{}]".format(category)
        flat_list.append((header, attrib_object.value))
        if attrib_object.term_source:
            if make_unique:
                flat_list.append((category + sep + "Term Source REF", attrib_object.term_source))
            else:
                flat_list.append(("Term Source REF", attrib_object.term_source))
            if attrib_object.term_accession:
                if make_unique:
                    flat_list.append((category + sep + "Term Accession Number", attrib_object.term_accession))
                else:
                    flat_list.append(("Term Accession Number", attrib_object.term_accession))
        if attrib_object.unit:
            flat_list.extend(flatten_unit(category, attrib_object.unit, make_unique=make_unique))
    return flat_list


def sample_attribute_to_dict(category, attrib_object):
    attrib_dict = OrderedDict()
    if attrib_object.value:
        header = "Characteristics[{}]".format(category)
        attrib_dict[header] = attrib_object.value
        if attrib_object.term_source:
            attrib_dict["Term Source REF"] = attrib_object.term_source
            if attrib_object.term_accession:
                attrib_dict["Term Accession Number"] = attrib_object.term_accession
        if attrib_object.unit:
            attrib_dict["unit"] = unit_to_dict(category, attrib_object.unit)
    return attrib_dict


def unit_to_dict(category, unit_object):
    unit_dict = OrderedDict()
    if unit_object.value and unit_object.unit_type:
        unit_header = "Unit[{}]".format(unit_object.unit_type)
        unit_dict[unit_header] = unit_object.value
        if unit_object.term_source:
            unit_dict["Term Source REF"] = unit_object.term_source
            if unit_object.term_accession:
                unit_dict["Term Accession Number"] = unit_object.term_accession
    return unit_dict


def flatten_nested_zip(nested_zip, flat_list):
    """Take a list of nested tuples and return a single list with all non-empty values."""
    if isinstance(nested_zip, tuple):
        for item in nested_zip:
            if isinstance(item, tuple):
                return flatten_nested_zip(item, flat_list)
            # Removes None values
            elif item:
                flat_list.append(item)
    else:
        flat_list.append(nested_zip)
    return flat_list


def unique_keys_in_order(rows):
    """Take a list of ordered dictionaries and return a list of all keys without duplicates and in original order."""
    # Zip all keys into a nested tuple
    zipped_keys = []
    for d in rows:
        zipped_keys = list(itertools.zip_longest(list(d), zipped_keys))
    # Turn nested tuples into a flat list
    all_keys = []
    for x in zipped_keys:
        all_keys.extend(flatten_nested_zip(x, []))
    # Remove duplicates
    unique_keys = list(OrderedDict.fromkeys(all_keys))
    return unique_keys


