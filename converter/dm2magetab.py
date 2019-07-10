"""Module to convert metadata in the submission data model to MAGE-TAB files."""

import pandas as pd

from collections import OrderedDict, defaultdict

from utils.converter_utils import get_controlled_vocabulary


def generate_idf(sub):
    """Transform study/project/protocol metadata in data model to an IDF vertical table."""

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
        ("SDRF File", sub.info.get("accession", sub.info.get("alias")) + ".sdrf.txt"),
        ("Term Source Name", ["EFO", "ArrayExpress"]),
        ("Term Source File", ["https://www.ebi.ac.uk/efo/", "https://www.ebi.ac.uk/arrayexpress/"]),
        ("Comment[AEExperimentType]", [exptype for exptype in sub.study.experiment_type]),
        ("Comment[ArrayExpressAccession]", sub.info.get("accession"))
    ])

    if sub.info.get("submission_type") == "sequencing" or sub.info.get("submission_type") == "singlecell":
        pass
        # TODO: Add comments for ENA accessions Comment[SecondaryAccession] and Comment[SequenceDataURI]

    return idf


def generate_sdrf(sub):
    """Transform sample and file metadata in data model to an SDRF table."""

    submission_type = sub.info.get("submission_type")
    protocol_positions = get_protocol_positions(submission_type)
    rows = []

    # For each node (sample, extract, assay etc.) start a list of tuples with category value pairs,
    # because each node block can have different attributes for each sample. Therefore all data points
    # are collected separately per node block and then merged at the end into one SDRF table.
    for sample in sub.sample:
        row = []
        all_protocols = []

        sample_values = [("Source Name", sample.alias)]
        if sample.accession:
            sample_values.append(("Comment[BioSD_SAMPLE]", sample.accession))
        # Expand sample attributes to characteristics columns (they can be different between different samples)
        for category, sample_attrib in sample.attributes.items():
            sample_values.extend(flatten_sample_attribute(category, sample_attrib, "Characteristics"))
        if sample.description:
            sample_values.append(("Description", sample.description))
        if sample.material_type:
            sample_values.append(("Material Type", sample.material_type))

        print(sample_values)
        # Add source node
        row.extend([OrderedDict(sample_values)])
        print(row)

        # Get all assay objects that belong to this sample
        assays = [assay for assay in sub.assay if assay.sampleref == sample.alias]

        for assay in assays:
            # Reformat protocol REFs for edges between nodes
            all_protocols.extend([sub.get_protocol(pref) for pref in assay.protocolrefs])
            protocol_refs = sort_protocol_refs_to_dict(protocol_positions, all_protocols)
            print(protocol_refs)

            extract_values = [("Extract Name", sample.alias)]

            if submission_type == "microarray":
                le_values = [
                    ("Labelled Extract Name", assay.alias),
                    ("Label", assay.label)]
                # Add protocol refs, extract node and labeled extract node
                row.extend([protocol_refs[1], OrderedDict(extract_values),
                            protocol_refs[2], OrderedDict(le_values)])

            # submission type is sequencing or singlecell, get assay attributes and convert them to comments
            else:
                assay_attributes = assay.get_attributes()
                for aa in assay_attributes:
                    extract_values.append(("Comment[{}]".format(aa.upper()), getattr(assay, aa)))
                row.extend([protocol_refs[1], OrderedDict(extract_values)])

            # Get all assay data objects that belong to this assay
            data = [ad for ad in sub.assay_data if assay.alias in ad.assayrefs]
            for ad in data:
                # For matrix files the Assay Name is inferred from the assay object (i.e. labeled extract name)
                if ad.data_type == "raw matrix":
                    assay_name = assay.alias
                else:
                    assay_name = ad.alias

                all_protocols.extend([sub.get_protocol(pref) for pref in ad.protocolrefs])
                protocol_refs = sort_protocol_refs_to_dict(protocol_positions, all_protocols)
                assay_values = [("Assay Name", assay_name),
                                ("Technology Type", assay.technology_type)]

                if submission_type == "microarray":
                    assay_values.extend([("Array Design REF", assay.array_design),
                                         ("array-design~~~Term Source REF", "Array Express")])
                    # Add Assay node
                    row.extend([protocol_refs[3], OrderedDict(assay_values)])
                else:
                    row.extend([protocol_refs[4], OrderedDict(assay_values)])

                # Get all data files
                data_values = []
                for f in ad.files:
                    if ad.data_type == "raw":
                        data_values.append(("Array Data File", f.name))
                    elif ad.data_type == "raw matrix":
                        data_values.append(("Array Data Matrix File", f.name))

                    processed_data_values, factor_values = configure_analysis_and_factors(all_protocols, ad, assay,
                                                                                          sample, sub)
                    protocol_refs = sort_protocol_refs_to_dict(protocol_positions, all_protocols)
                    # Add data and factor nodes
                    rows.append(row +
                                [OrderedDict(data_values),
                                 protocol_refs[6],
                                 OrderedDict(processed_data_values),
                                 OrderedDict(factor_values)])

                    #append_sdrf_row(submission_type, rows, protocol_refs, sample_values, extract_values)

                    # Add all lists together to form the complete row of the SDRF. We do this
                    # at the level of raw data files, which means each raw data file gets a row.
                    # All lists are transformed to dictionaries so that we can use pandas to write a table.

                if not ad.files:
                    # Haven't found any raw data files, checking processed data and factors
                    processed_data_values, factor_values = configure_analysis_and_factors(all_protocols, assay, sample, sub)
                    protocol_refs = sort_protocol_refs_to_dict(protocol_positions, all_protocols)
                    row.extend([protocol_refs[6], OrderedDict(processed_data_values),
                                OrderedDict(factor_values)])
                    rows.append(row)

            if not data:
                # Haven't found any raw data, checking processed data and factors
                processed_data_values, factor_values = configure_analysis_and_factors(all_protocols, None, assay,
                                                                                      sample, sub)
                protocol_refs = sort_protocol_refs_to_dict(protocol_positions, all_protocols)
                row.extend([protocol_refs[6], OrderedDict(processed_data_values),
                            OrderedDict(factor_values)])
                rows.append(row)

        # Haven't found any assays, writing sample info only
        if not assays:
            processed_data_values, factor_values = configure_analysis_and_factors(all_protocols, None, None, sample, sub)
            protocol_refs = sort_protocol_refs_to_dict(protocol_positions, all_protocols)
            row.extend([protocol_refs[6], OrderedDict(processed_data_values),
                        OrderedDict(factor_values)])
            rows.append(row)

        # Append the complete row to the table

                    #if submission_type == "microarray":
                    #    rows.append([OrderedDict(sample_values), protocol_refs[1],
                    #                 OrderedDict(extract_values), protocol_refs[2],
                    #                 OrderedDict(le_values), protocol_refs[3],
                    #                 OrderedDict(assay_values), protocol_refs[5],
                    #                 OrderedDict(data_values), protocol_refs[6],
                    #                 OrderedDict(processed_data_values),
                    #                 OrderedDict(factor_values)])
                    #else:  # submission type is sequencing or singlecell
                    #    rows.append([OrderedDict(sample_values), protocol_refs[1],
                    #                 OrderedDict(extract_values), protocol_refs[4],
                    #                 OrderedDict(assay_values), protocol_refs[5],
                    #                 OrderedDict(data_values), protocol_refs[6],
                    #                 OrderedDict(processed_data_values),
                    #                 OrderedDict(factor_values)])
    print(rows)
    # This goes through the collection of ordered dictionaries and transforms them into pandas data frames,
    # while merging the nodes/attributes for different samples, e.g. all extract attributes from all samples together
    data_frames = []
    for i in range(len(rows[0])):
        data_frames.append(pd.DataFrame.from_records([row[i] for row in rows]))

    # Pandas concat merges the dictionaries for the different SDRF parts (nodes) together into one big table
    raw_sdrf = pd.concat(data_frames, axis=1)

    # Raw output still has "uniquified" column headers
    return raw_sdrf


def configure_analysis_and_factors(all_protocols, assay_data, assay, sample, sub):
    # Processed data
    processed_data = []
    if assay_data:
        processed_data = [px for px in sub.analysis if assay_data.alias in px.assaydatarefs]
    # Get all analysis objects that belong to this assay data object
    elif assay:
        processed_data = [px for px in sub.analysis if assay.alias in px.assayrefs]
    elif sample:
        processed_data = [px for px in sub.analysis if sample.alias in px.samplerefs]
    # Collect file names and turn into tuple list
    processed_data_values = []
    for px in processed_data:
        for f in px.files:
            if px.data_type == "processed":
                processed_data_values.append(("Derived Array Data File", f.name))
            elif px.data_type == "processed matrix":
                processed_data_values.append(("Derived Array Data Matrix File", f.name))
            if f.ftp_location:
                processed_data_values.append(("Comment[Derived ArrayExpress FTP file]", f.ftp_location))
        # Also add protocol references for how the processed data was generated from assay data
        all_protocols.extend([sub.get_protocol(pref) for pref in px.protocolrefs])
    # Factor values
    factors = sub.study.experimental_factor
    factor_values = []
    # Look up factor in sample attributes and turn into ordered dict with unit/term columns
    for f in factors:
        factor_values.extend(
            flatten_sample_attribute(f.value, sample.attributes.get(f.value), "Factor Value"))

    return processed_data_values, factor_values


def append_sdrf_row(submission_type, rows, protocol_refs, sample_values=[], extract_values=[], le_values=[],
                    assay_values=[], data_values=[], processed_data_values=[], factor_values=[]):

    if submission_type == "microarray":
        rows.append([OrderedDict(sample_values), protocol_refs[1],
                 OrderedDict(extract_values), protocol_refs[2],
                 OrderedDict(le_values), protocol_refs[3],
                 OrderedDict(assay_values), protocol_refs[5],
                 OrderedDict(data_values), protocol_refs[6],
                 OrderedDict(processed_data_values),
                 OrderedDict(factor_values)])
    else:  # submission type is sequencing or singlecell
        rows.append([OrderedDict(sample_values), protocol_refs[1],
                 OrderedDict(extract_values), protocol_refs[4],
                 OrderedDict(assay_values), protocol_refs[5],
                 OrderedDict(data_values), protocol_refs[6],
                 OrderedDict(processed_data_values),
                 OrderedDict(factor_values)])


def flatten_unit(category, unit_object, make_unique=True, sep="~~~"):
    """
    Transform a unit object with type and references to ontology terms into a list of key/value tuples
    with defined order (as in the SDRF).

    :param category: string, the type of attribute the unit belongs to, e.g. "age"
    :param unit_object: Unit object, the object holiding the unit's attributes
    :param make_unique: boolean, append name of category to make the key unique
    :param sep: string, separator between category name and SDRF header in unique header values
    :return: list of tuples
    """
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
    """
    Transform the attribute and unit objects with references to ontology terms
    into a list of key/value tuples with defined order (as in the SDRF).
    Works for Characteristics and Factor value types (this can be specified with "column header".

    :param category: string, the type of sample attribute
    :param attrib_object: Attribute object, the object holding the category's attributes
    :param column_header: string, Field name, either "Characteristics" or "Factor Value"
    :param make_unique: boolean, append name of category to make the key unique
    :param sep: string, separator between category name and SDRF header in unique header values
    :return: list of tuples
    """
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
    print(flat_list)
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

    :return: nested dictionary, with dict default value
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


def column_name_to_magetab(header, sep="~~~"):
    """Transform unique column header back to MAGE-TAB style.
    Column headers are expected to be separated by sep. """
    header_list = header.split(sep)
    return header_list[-1]


def write_sdrf_file(pandas_table, new_file_name, logger):
    """Write out SDRF tab-delimited text file from merged pandas table

    :param pandas_table: pandas data frame containing unique column headers for each SDRF column
    :param new_file_name: file path to write SDRF
    :param logger: log for errors
    :return: None

    """
    logger.debug("Renaming unique column headers back to MAGE-TAB format.")
    pandas_table.rename(column_name_to_magetab, axis="columns", inplace=True)

    logger.debug("Writing new SDRF {}".format(new_file_name))
    try:
        pandas_table.to_csv(new_file_name, sep='\t', encoding='utf-8', index=False)
    except Exception as e:
        logger.error("Failed to write SDRF: {}".format(str(e)))

