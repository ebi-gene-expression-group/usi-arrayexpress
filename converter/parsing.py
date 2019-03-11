from collections import defaultdict, OrderedDict

from utils.converter_utils import get_controlled_vocabulary, get_name, get_value, read_sdrf_file, read_idf_file


def get_protocol_refs(sdrf_row, header_dict, node_map, node2_index):
    """Collect the protocol refs between the given node and the previous one"""
    ref_list = []
    largest_index = 0
    if "protocolref" in header_dict:
        locations = header_dict["protocolref"]

        for nodes in node_map.values():
            node_end = nodes[0][1]
            if largest_index < node_end < node2_index:
                largest_index = nodes[0][1]

        ref_list = [sdrf_row[i] for i in locations if largest_index < i > node2_index]
    return ref_list


def get_comment_values(sdrf_row, header, node1_index, node2_index):
    """Collect all comment columns between the given nodes
    and store the values in the given row as a dictionary."""
    comments_dict = {}
    for i in range(node1_index, node2_index + 1):
        if get_name(header[i]) == 'comment':
            comments_dict[get_value(header[i])] = sdrf_row[i]
    return comments_dict


def get_node_positions(nodes, header):
    """Take an SDRF header row as list and return the start and end indexes of the nodes
    and the indexes of the preceeding protocol refs.

    Go through the header and find the start and end indexes for all nodes
    We can have more than one node of the same type, e.g. derived array data file.
    Duplicated nodes will create additional entries in the list value of the dict.

    The output is a dictionary of the node names as keys
    keys = node names
    values = a list of lists (one list for each node)
        each list has three values:
        index 0 = (type: int) header index of node name column, e.g. Source Name
        index 1 = (type: int) header index of the column before the next Protocol REF (or Factor Value at the end)
        index 3 = (type: list) header indexes of the Protocol REF column preceeding the node name column
    """

    node_breakpoints = defaultdict(list)
    protocols = []
    last_node = None
    node_open = False

    for i, h in enumerate(header):
        h = get_name(h)
        if h in nodes and not node_open:
            last_node = h
            node_breakpoints[last_node].append([i, None, protocols])
            node_open = True
            protocols = []
        elif h == "protocolref":
            protocols.append(i)
            if node_open:
                node_breakpoints[last_node][-1][1] = i - 1
                node_open = False
        # In case there is no protocol between the nodes
        elif h in nodes and node_open:
            node_breakpoints[last_node][-1][1] = i - 1
            last_node = h
            node_breakpoints[last_node].append([i, None, []])
        elif h == 'factorvalue' and node_open:
            node_breakpoints[last_node][-1][1] = i - 1
            node_open = False
    # To catch the last one in case there are no factor values
    if node_breakpoints[last_node][-1][1] is None:
        node_breakpoints[last_node][-1][1] = i

    return node_breakpoints


def parse_sdrf(sdrf_file):
    """
    Read SDRF data table and return dictionaries for the different nodes

    The function returns a dictionary for each node (samples, extracts, labeled extracts, assays,
    raw_data, processed_data). Each dict contains all unique entries and their attributes and
    links to other nodes.

    :param sdrf_file: string, path to SDRF file
    """
    sample_nodes = ("sourcename",)
    le_nodes = ("labeledextractname",)
    extract_nodes = ("extractname",)
    assay_nodes = ("assayname",)
    raw_data_nodes = ("arraydatafile", "arraydatamatrixfile", "scanname",)
    processed_data_nodes = ("derivedarraydatafile", "derivedarraydatamatrixfile",)

    nodes = sample_nodes + extract_nodes + le_nodes + assay_nodes + raw_data_nodes + processed_data_nodes

    # Read in the file and get data (separated from header), the header row, and a breakdown of header nodes/attributes
    sdrf_data, header, header_dict = read_sdrf_file(sdrf_file)

    # A map of the start, end and protocol refs of each node
    node_map = get_node_positions(nodes, header)

    samples = OrderedDict()
    extracts = OrderedDict()
    le = OrderedDict()
    assays = OrderedDict()
    raw_data = OrderedDict()
    processed_data = OrderedDict()

    # Recording the node names, to skip duplicate rows
    sample_names = []
    extract_names = []
    le_names = []
    assay_names = []

    for row in sdrf_data:

        sample_name = ""
        extract_name = ""
        le_name = ""
        assay_name = ""

        # Samples

        last_attribute = None
        last_unit = None
        last_termsource = None

        sample_attributes = {"name": "",
                             "characteristics": defaultdict(list),
                             "material_type": "",
                             "description": "",
                             "comments": {},
                             "factors": {}
                             }

        # Only one Source Name column is allowed
        if "sourcename" in node_map and len(node_map["sourcename"]) == 1:
            node_range = node_map["sourcename"][0]
            sample_name = row[node_range[0]]
            # Skipping the samples we have already seen
            if sample_name not in sample_names:
                sample_attributes["name"] = sample_name
                sample_names.append(sample_name)
                # Go through the header values in between the nodes
                for i in range(node_range[0] + 1, node_range[1] + 1):
                    # Get Characteristics
                    if "characteristics" in header_dict and i in header_dict["characteristics"]:
                        last_attribute = get_value(header[i])
                        sample_attributes["characteristics"][last_attribute] = {"value": row[i]}
                    # Add units
                    elif "unit" in header_dict and i in header_dict["unit"]:
                        if i - 1 in header_dict["characteristics"] and last_attribute:
                            last_unit = get_value(header[i])
                            sample_attributes["characteristics"][last_attribute]["unit"] = {"value": row[i]}
                            sample_attributes["characteristics"][last_attribute]["unit"]["unit_type"] = last_unit
                        else:
                            print("PARSER ERROR: [column {}] Unit found without Characteristics.".format(i + 1))
                    # Add Term Source REFs
                    elif "termsourceref" in header_dict and i in header_dict["termsourceref"]:
                        if "characteristics" in header_dict and i - 1 in header_dict["characteristics"] and last_attribute:
                            last_termsource = get_name(header[i])
                            sample_attributes["characteristics"][last_attribute]["term_source"] = row[i]
                        elif "unit" in header_dict and i - 1 in header_dict["unit"] and last_attribute and last_unit:
                            last_termsource = get_name(header[i])
                            sample_attributes["characteristics"][last_attribute]["unit"]["term_source"] = row[i]
                        else:
                            print(
                                "PARSER ERROR: [colum {}] \"Term source REF\" found without Characteristics or Unit".format(
                                    i + 1))
                    # Add Term Accession Number
                    elif "termaccessionnumber" in header_dict and i in header_dict["termaccessionnumber"]:
                        if "termsourceref" in header_dict and i - 1 in header_dict["termsourceref"] and last_termsource:
                            if "characteristics" in header_dict and i - 2 in header_dict["characteristics"] and last_attribute:
                                sample_attributes["characteristics"][last_attribute]["term_accession"] = row[i]
                            elif "unit" in header_dict and i - 2 in header_dict[
                                "unit"] and last_attribute and last_unit:
                                sample_attributes["characteristics"][last_attribute]["unit"]["term_accession"] = row[i]
                        else:
                            print(
                                "PARSER ERROR: [column {}] \"Term Accession Number\" found without Term Source REF".format(
                                    i + 1))
                    # Get Material Type
                    elif "materialtype" in header_dict and i in header_dict["materialtype"]:
                        sample_attributes["material_type"] = row[i]
                    # Get Description
                    elif "description" in header_dict and i in header_dict["description"]:
                        sample_attributes["description"] = row[i]
                    # Get Comments
                    sample_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1])

                samples[sample_name] = sample_attributes

        # Extract

        extract_attributes = {"name": "",
                              "comments": {},
                              "sample_ref": "",
                              "protocol_ref": []}

        if "extractname" in header_dict:
            # This creates separate extract dicts if there is more than one column. Should we allow that?
            node_range = node_map["extractname"][0]

            # Expecting the first value to the be node name
            extract_name = row[node_range[0]]

            if extract_name not in extract_names:
                extract_names.append(extract_name)
                extract_attributes["name"] = extract_name
                # Get comments
                extract_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1])
                # Keep reference to sample in that row
                extract_attributes["sample_ref"] = sample_name

                extract_attributes["protocol_ref"] = [row[i] for i in node_range[2]]
                extracts[extract_name] = extract_attributes

        # Labeled Extract

        le_attributes = {"name": "",
                         "label": "",
                         "comments": {},
                         "extract_ref": [],
                         "protocol_ref": []}

        if "labeledextractname" in header_dict:

            # Only one labeled extract (le) column is allowed
            node_range = node_map["labeledextractname"][0]
            le_name = row[node_range[0]]
            if le_name not in le_names:
                le_names.append(le_name)
                le_attributes["name"] = le_name
                # Get label (only one label column allowed)
                if "label" in header_dict:
                    le_attributes["label"] = row[header_dict["label"][0]]
                # Get comments
                le_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1])
                # Keep reference to sample in that row
                le_attributes["extract_ref"].append(extract_name)

                le_attributes["protocol_ref"] = [row[i] for i in node_range[2]]
                le[le_name] = le_attributes

        # Assay

        assay_attributes = {"name": "",
                            "technology_type": "",
                            "comments": {},
                            "extract_ref": [],
                            "protocol_ref": [],
                            "array_design": ""}

        if "assayname" in header_dict:
            node_range = node_map["assayname"][0]
            assay_name = row[node_range[0]]

            if assay_name not in assay_names:
                assay_attributes["name"] = assay_name
                assay_names.append(assay_name)

                for i in range(node_range[0] + 1, node_range[1] + 1):
                    if get_name(header[i]) == "technologytype":
                        assay_attributes["technology_type"] = row[i]
                    elif get_name(header[i]) == "arraydesignref":
                        assay_attributes["array_design"] = row[i]
                assay_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1])
                assay_attributes["protocol_ref"] = [row[i] for i in node_range[2]]
                assays[assay_name] = assay_attributes

                if le_name:
                    assay_attributes["extract_ref"].append(le_name)
                elif extract_name:
                    assay_attributes["extract_ref"].append(extract_name)
            # We allow more than 1 extract per assay for the two-colour case
            else:
                # Add the extract ref to the previous assay, using set to avoid duplicates
                if le_name:
                    assays[assay_name]["extract_ref"].append(le_name)
                elif extract_name:
                    assays[assay_name]["extract_ref"].append(extract_name)

        # Datafiles

        # Raw data
        parse_data_file_columns(raw_data_nodes, header_dict, header, node_map, row, raw_data,
                                sample_name, extract_name, le_name, assay_name)

        # Processed data

        parse_data_file_columns(processed_data_nodes, header_dict, header, node_map, row,
                                processed_data, sample_name, extract_name, le_name, assay_name)


        # Factor Values

        # For now we'll assume that each factor value corresponds to a biological replicate, i.e. sample,
        # and factor values get added to the sample attributes
        # This means this won't work for the edge cases where a technical parameter is used as factor,
        # e.g. comparison of the same sample sequenced with different sequencers

        factors = OrderedDict()

        for i in header_dict.get("factorvalue", []):
            # Get value and category
            factor_type = get_value(header[i])
            factors[factor_type] = {"value": row[i]}
            # Getting units (simplified)
            if i+1 in header_dict.get("unit", []):
                factors[factor_type]["unit"] = {"value": row[i+1],
                                                "unit_type": get_value(header[i + 1])}
        # Add to the sample attributes of the current sample
        # Note this will deliberately overwrite if there are different values within one sample (see comment above)
        samples[sample_name]["factors"] = factors

    return samples, extracts, le, assays, raw_data, processed_data


def parse_data_file_columns(data_nodes, header_dict, header, node_map, row, data_dict,
                            sample_name, extract_name, le_name, assay_name):
    """Parse data file columns and add file attributes to the respective dict.

    This function directly modifies the dictionary that is passed in.
    It works as part of the parse_sdrf function, using the same strategy
    to parse raw data and processed data nodes."""

    file_attributes = {
        "name": "",
        "data_type": "",
        "assay_ref": [],
        "le_ref": [],
        "extract_ref": [],
        "sample_ref": [],
        "protocol_ref": [],
        "comments": []
    }

    for rdn in data_nodes:
        if rdn in header_dict:
            # To make sure we can handle multiple columns of the same type
            for node_range in node_map[rdn]:
                node_index = node_range[0]
                file_name = row[node_index]
                # Initialise dict entry for each new file
                if file_name not in data_dict:
                    data_dict[file_name] = file_attributes
                    data_dict[file_name]["name"] = file_name
                    data_dict[file_name]["data_type"] = rdn
                    data_dict[file_name]["comments"] = get_comment_values(row, header, node_range[0], node_range[1])
                    data_dict[file_name]["assay_ref"].append(assay_name)
                    data_dict[file_name]["extract_ref"].append(extract_name)
                    data_dict[file_name]["sample_ref"].append(sample_name)
                    data_dict[file_name]["protocol_ref"].extend([row[i] for i in node_range[2]])
                    if le_name:
                        data_dict[file_name]["le_ref"].append(le_name)

                # For microarray we can have several labeled extracts/assays per file
                else:
                    # This may be redundant but recording this for now
                    if le_name and le_name not in data_dict[file_name]["le_ref"]:
                        data_dict[file_name]["le_ref"].append(le_name)
                    if assay_name not in data_dict[file_name]["assay_ref"]:
                        data_dict[file_name]["assay_ref"].append(assay_name)
                    if extract_name not in data_dict[file_name]["extract_ref"]:
                        data_dict[file_name]["extract_ref"].append(extract_name)
                    if sample_name not in data_dict[file_name]["sample_ref"]:
                        data_dict[file_name]["sample_ref"].append(sample_name)


def parse_idf(idf_file):
    """This parses the raw IDF dictionary and puts the fields and values into a sub-category dictionary."""

    idf_dict = read_idf_file(idf_file)

    study_info = {
        "title": "",
        "accession": "",
        "description": "",
        "experimental_design": [],
        "experimental_factor": [],
        "releaseDate": "",
        "date_of_experiment": "",
        "contacts": [],
        "publications": [],
        "comments": {},
        "protocolRefs": [],
        "idf_filename": idf_file,
        "sdrf_filename": ""
    }

    protocols = list()

    # Contacts
    parse_multi_column_fields(idf_dict, study_info["contacts"], "contact_terms")

    # Publications
    parse_multi_column_fields(idf_dict, study_info["publications"], "publication_terms")

    # Factors
    parse_multi_column_fields(idf_dict, study_info["experimental_factor"], "factor_terms")

    # Experimental Design
    parse_multi_column_fields(idf_dict, study_info["experimental_design"], "design_terms")

    # Protocols
    parse_multi_column_fields(idf_dict, protocols, "protocol_terms")

    study_info["protocolRefs"] = [p.get("title") for p in protocols]

    # Comments
    # Here we allow >1 value but they are not related to the other comment values in the same column
    comment_terms = get_controlled_vocabulary("idf_comment_terms")
    for idf_ct, usi_ct in comment_terms.items():
        if idf_ct in idf_dict:
            # Remove empty list entries, e.g. 'related experiment': ['E-MTAB-7236', '', '', '']
            comment_values = [x for x in idf_dict[idf_ct] if x]
            # A new dict entry for the term with the list of values
            study_info["comments"][usi_ct] = comment_values

    # General Info
    general_terms = {get_name(t): val for t, val in get_controlled_vocabulary("investigation_terms").items()}
    for idf_ct, usi_ct in general_terms.items():
        if idf_ct in idf_dict and idf_dict[idf_ct]:
            # for these terms we only expect/allow 1 value (the first item in the list)
            study_info[usi_ct] = idf_dict[idf_ct][0]

    accession = study_info["comments"].get("accession", None)
    if accession:
        study_info["comments"]["accession"] = accession[0]

    return study_info, protocols


def parse_multi_column_fields(idf_dict, category_list, lookup_term):
    """Helper function to connect the related entries in IDF columns,
    e.g. for protocols, contacts and factor values.
    It translates the IDF terms into terms used in the USI data model
    and saves the values as a list of individual dictionaries."""

    # Get the field names and translation from IDF to USI
    controlled_terms = {get_name(t): val for t, val in get_controlled_vocabulary(lookup_term).items()}
    # Go through terms
    for idf_ct, usi_ct in controlled_terms.items():
        if idf_ct in idf_dict:
            # Go through IDF value lists
            for i, contact_value in enumerate(idf_dict[idf_ct]):
                # Not for empty values
                if contact_value:
                    # Add new list entry for the current position
                    if len(category_list) <= i:
                        category_list.append(dict())
                    # Fill in value for this term and position
                    category_list[i][usi_ct] = contact_value


