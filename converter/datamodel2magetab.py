"""Module to convert metadata in the submission data model to MAGE-TAB files."""

from collections import OrderedDict, defaultdict
from utils import converter_utils
import codecs
import csv
import pandas as pd
import itertools
import unittest


from utils.converter_utils import get_controlled_vocabulary


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

    rows = []

    submission_type = sub.info.get("submission_type")

    protocol_positions = get_protocol_positions(submission_type)

    if submission_type == "microarray":

        # For each node (sample, extract, assay etc.) start a list of tuples with category value pairs,
        # because each node block can have different attributes for each sample. Therefore all data points
        # are collected separately per node block and then merged at the end into one SDRF table.
        for sample in sub.sample:
            sample_values = [
                ("Source Name", sample.alias)
            ]
            # Expand sample attributes to characteristics columns (they can be different between different samples)
            for category, sample_attrib in sample.attributes.items():
                sample_values.extend(flatten_sample_attribute(category, sample_attrib, make_unique=True))

            if sample.description:
                sample_values.append(("Description", sample.description))
            if sample.material_type:
                sample_values.append(("Material Type", sample.material_type))

            # Get all assay objects that belong to this sample
            assays = [assay for assay in sub.assay if assay.sampleref == sample.alias]

            for assay in assays:

                source_to_extract_protocols = [
                    ("growth~~~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                               if sub.get_protocol(pref).protocol_type.value == "growth protocol"]),
                    ("treatment~~~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                                  if sub.get_protocol(pref).protocol_type.value == "treatment protocol"]),
                    ("dissection~~~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                                   if sub.get_protocol(pref).protocol_type.value == "dissection protocol"]),
                    ("collection~~~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                                   if sub.get_protocol(pref).protocol_type.value == "sample collection protocol"]),
                    ("conversion~~~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                                   if sub.get_protocol(pref).protocol_type.value == "conversion protocol"])
                ]

                assay_values = [
                    ("Extract Name", sample.alias),
                    ("labeling~~~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
                                      if sub.get_protocol(pref).protocol_type.value == "nucleic acid labeling protocol"]),
                    ("Labelled Extract Name", assay.alias),
                    ("hyb~~~Protocol REF", [sub.get_protocol(pref).alias for pref in assay.protocolrefs
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


    print(rows)

    data_frames = []
    for i in range(len(rows[0])):
        data_frames.append(pd.DataFrame.from_records([row[i] for row in rows]))

    raw_sdrf = pd.concat(data_frames, axis=1)

    return raw_sdrf

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
    print(attrib_object)
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


def get_protocol_positions(techtype):
    """Fetch all protocol types for microarray/sequencing studies and return a dictionary sorted by position in the SDRF.
    {1: [sample collection, growth, treatment], 2: [labeling], 3: [hybridization] ...}"""
    # Use the same protocols for singlecell as for sequencing
    if techtype == "singlecell":
        techtype = "sequencing"
    protocol_info = get_controlled_vocabulary("protocol_types", resource="ontology")
    protocol_positions = defaultdict(list)
    for p_type, p_info in protocol_info.items():
        if p_info["exp_type"] == "all" or p_info["exp_type"] == techtype:
            protocol_positions[p_info["position"]].append(p_type)
    return protocol_positions



