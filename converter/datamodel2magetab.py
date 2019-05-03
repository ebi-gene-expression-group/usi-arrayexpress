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

    submission_type = sub.info.get("submission_type")
    protocol_positions = get_protocol_positions(submission_type)
    rows = []

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
                sample_values.extend(flatten_sample_attribute(category, sample_attrib, "Characteristics"))

            if sample.description:
                sample_values.append(("Description", sample.description))
            if sample.material_type:
                sample_values.append(("Material Type", sample.material_type))

            # Get all assay objects that belong to this sample
            assays = [assay for assay in sub.assay if assay.sampleref == sample.alias]
            all_protocols = []

            for assay in assays:
                all_protocols.extend([sub.get_protocol(pref) for pref in assay.protocolrefs])

                extract_values = [("Extract Name", sample.alias)]

                le_values = [
                    ("Labelled Extract Name", assay.alias),
                    ("Label", assay.label)]

                # Get all assay data objects that belong to this assay
                data = [ad for ad in sub.assay_data if assay.alias in ad.assayrefs]

                for ad in data:
                    # For matrix files the Assay Name is inferred from the assay object (i.e. labeled extract name)
                    if ad.data_type == "raw matrix":
                        assay_name = assay.alias
                    else:
                        assay_name = ad.alias

                    assay_values = [
                        ("Assay Name", assay_name),
                        ("Technology Type", assay.technology_type),
                        ("Array Design REF", assay.array_design),
                        ("array-design~~~Term Source REF", "Array Express")]

                    all_protocols.extend([sub.get_protocol(pref) for pref in ad.protocolrefs])

                    # Get all data files
                    data_values = []
                    for f in ad.files:
                        if ad.data_type == "raw":
                            data_values.append(("Array Data File", f.name))
                        elif ad.data_type == "raw matrix":
                            data_values.append(("Array Data Matrix File", f.name))

                    # Factor values
                    factors = sub.study.experimental_factor

                    factor_values = []
                    # Look up factor in sample attributes and turn into ordered dict with unit/term columns
                    for f in factors:
                        factor_values.extend(flatten_sample_attribute(f.value, sample.attributes.get(f.value), "Factor Value"))

                    # Reformat protocol REFs for edges between nodes
                    protocol_refs = sort_protocol_refs_to_dict(protocol_positions, all_protocols)

                    # Add all lists together and transform to dictionary so that we can use pandas to write to table
                    rows.append([OrderedDict(sample_values), protocol_refs[1],
                                 OrderedDict(extract_values), protocol_refs[2],
                                 OrderedDict(le_values), protocol_refs[3],
                                 OrderedDict(assay_values), protocol_refs[5],
                                 OrderedDict(data_values),
                                 OrderedDict(factor_values)
                                 ])

    data_frames = []
    for i in range(len(rows[0])):
        data_frames.append(pd.DataFrame.from_records([row[i] for row in rows]))

    raw_sdrf = pd.concat(data_frames, axis=1)

    return raw_sdrf


def flatten_unit(category, unit_object, make_unique=True, sep="~~~"):
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


def flatten_sample_attribute(category, attrib_object, column_header, make_unique=True, sep="~~~"):
    flat_list = []
    if attrib_object.value:
        header = "{}[{}]".format(column_header, category)
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


def sort_protocol_refs_to_dict(protocol_positions, all_protocols, sep="~~~"):
    """
    Convert protocol references to a dictionary style to be transformed to columns
    :param protocol_positions: dictionary with the position as key and list of protocol types as value
    :param all_protocols: list of protocol objects
    :param sep: separator used to make dict keys (later column names unique)

    :return: nested dictionary
    {1: {11~~~Protocol REF: Protocol 1}, {12~~~Protocol REF: Protocol 2}
     2: {21~~~Protocol REF: Protocol 3},
     3: {31~~~Protocol REF: Protocol 4}}
    """
    protocol_dict = defaultdict(dict)

    for pos, p_types in protocol_positions.items():
        prefs_for_position = [p.alias for p in all_protocols if p.protocol_type.value in p_types]
        # Number of entries in the dict corresponds to the number of columns that will be created and
        # should be equal of the number of protocol refs for the same position
        column_number = 1
        for p in prefs_for_position:
            # Making the secondary key unique
            protocol_dict[pos][str(pos)+str(column_number)+sep+"Protocol REF"] = p
            column_number += 1

    return protocol_dict


