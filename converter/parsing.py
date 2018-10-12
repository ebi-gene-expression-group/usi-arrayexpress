import json

import codecs
import re
from collections import defaultdict


def read_json_file(filename):
    try:
        with open(filename) as fh:
            data = json.load(fh, encoding="utf-8")
            return data
    except IOError as err:
        print("Cannot import file: %s" % err)
        # TODO: system exit
    except ValueError as j_err:
        print("Cannot read JSON file: %s" % j_err)
        raise


def read_sdrf_file(sdrf_file):
    with codecs.open(sdrf_file, encoding='utf-8') as sf:
        header = sf.readline().rstrip().split('\t')
        sdrf_raw = sf.readlines()

    sdrf_list = [x.rstrip('\n').split('\t') for x in sdrf_raw]
    header_dict = defaultdict(list)
    for i, field in enumerate(header):
        short_name = get_name(field)
        header_dict[short_name].append(i)

    return sdrf_list, header, header_dict


def get_name(header_string):
    no_spaces = header_string.replace(' ', '')
    field_name = no_spaces.split('[')[0]
    return field_name.lower()


def get_value(header_string):
    field_value = header_string.split('[')[-1]
    return field_value.strip(']')


def get_prefix(filename):
    """This takes a filename string as input and strips off the file extension and any patterns to be ignored"""
    extensions = ('\.fastq\.gz$', '\.fq\.gz$', '\.txt\.gz$', '\.fastq\.bz2$', '\.[a-zA-Z0-9]+$', )
    ignore = ('_?001$', )
    filebase = None

    for ext in extensions:
        if re.search(ext, filename):
            filebase = re.sub(ext, '', filename)
            break
    for ip in ignore:
        if re.search(ip, filebase):
            filebase = re.sub(ip, '', filebase)
            return filebase
    if filebase:
        return filebase
    else:
        return filename


def get_protocol_refs(sdrf_row, header_dict, node_map, node2_index):
    """This collects the protocol refs between the given node and the previous one"""
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
    """This collects all comment columns between the given nodes
    and stores the values in the given row as a dictionary."""
    comments_dict = {}
    for i in range(node1_index, node2_index+1):
        if get_name(header[i]) == 'comment':
            comments_dict[get_value(header[i])] = sdrf_row[i]
    return comments_dict


def get_node_positions(nodes, header):
    """This function takes an SDRF header row as list
    and returns the start and end indexes of the nodes
    and the indexes of the protocol refs together with
    the following node."""

    # Go through the header and find the start and end indexes for all nodes
    # We can have more than one node of the same type, e.g. derived array data file
    # Duplicated nodes will create additional entries in the list value of the dict
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
            # In case there is no protocol between them
        elif h in nodes and node_open:
            node_breakpoints[last_node][-1][1] = i - 1
            last_node = h
            node_breakpoints[last_node].append([i, None, protocols])
        elif h == 'factorvalue' and node_open:
            node_breakpoints[last_node][-1][1] = i - 1
            node_open = False
    # To catch the last one in case there are no factor values
    if node_breakpoints[last_node][-1][1] is None:
        node_breakpoints[last_node][-1][1] = i

    return node_breakpoints


def parse_sdrf(sdrf_file):

    sample_nodes = ("sourcename",)
    le_nodes = ("labeledextractname",)
    extract_nodes = ("extractname",)
    assay_nodes = ("assayname",)
    raw_data_nodes = ("arraydatafile", "arraydatamatrixfile", "scanname",)
    processed_data_nodes = ("derivedarraydatafile", "derivedarraydatamatrixfile",)

    nodes = sample_nodes + extract_nodes + le_nodes + assay_nodes + raw_data_nodes + processed_data_nodes

    # Read in the file and get data (separated from header), the header row, and a breakdown of header nodes/attributes
    sdrf_data, header, header_dict = read_sdrf_file(sdrf_file)

    node_map = get_node_positions(nodes, header)

    samples = dict()
    extracts = dict()
    le = dict()
    assays = dict()
    raw_data = dict()
    processed_data = dict()


    # Guessing experiment type (SEQ or MA) from SDRF header
    is_microarray = False
    if 'arraydesignref' in header_dict or 'labeledextractname' in header_dict:
        is_microarray = True


    # Recording the node names, to skip duplicate rows
    sample_names = []
    extract_names = []
    le_names = []
    assay_names = []

    for row in sdrf_data:
        # new_sample = Sample.from_sdrf()

        #if is_microarray:
        #    sdr.assays.append(MicroarrayAssay.from_sdrf(row, header_dict, header))
        #else:
        #    sdr.assays.append(SeqAssay.from_sdrf(row, header_dict, header))

        #sdr.rawdata.append(RawData.from_sdrf(row, header_dict, header))

    #sdr.group_pairedend_files()

        sample_name = ""
        extract_name = ""
        le_name = ""
        assay_name = ""


        #Samples

        last_attribute = None
        last_unit = None
        last_termsource = None

        sample_attributes = {"name": "",
                             "characteristics": defaultdict(list),
                             "material type": "",
                             "description": "",
                             "comments": {}}

        # Only one Source Name column is allowed
        if "sourcename" in node_map and len(node_map["sourcename"]) == 1:
            node_range = node_map["sourcename"][0]
            sample_name = row[node_range[0]]
            # Skipping the samples we have already seen
            if sample_name not in sample_names:
                sample_attributes["name"] = sample_name
                sample_names.append(sample_name)
                # Go through the header values in between the nodes
                for i in range(node_range[0]+1, node_range[1]+1):
                    # Get Characteristics
                    if "characteristics" in header_dict and i in header_dict["characteristics"]:
                        last_attribute = get_value(header[i])
                        sample_attributes["characteristics"][last_attribute] = {"value": row[i]}
                    # Add units
                    elif "unit" in header_dict and i in header_dict["unit"]:
                        if i - 1 in header_dict["characteristics"] and last_attribute:
                            last_unit = get_value(header[i])
                            sample_attributes["characteristics"][last_attribute]["unit"] = {"value": row[i]}
                            sample_attributes["characteristics"][last_attribute]["unit"]["unit type"] = last_unit
                        else:
                            logger.log(40, "PARSER ERROR: [column {}] Unit found without Characteristics.".format(i+1))
                    # Add Term Source REFs
                    elif "termsourceref" in header_dict and i in header_dict["termsourceref"]:
                        if "characteristics" in header_dict and i-1 in header_dict["characteristics"] and last_attribute:
                            last_termsource = get_name(header[i])
                            sample_attributes["characteristics"][last_attribute]["term source"] = row[i]
                        elif "unit" in header_dict and i-1 in header_dict["unit"] and last_attribute and last_unit:
                            last_termsource = get_name(header[i])
                            sample_attributes["characteristics"][last_attribute]["unit"]["term source"] = row[i]
                        else:
                            logger.log(40, "PARSER ERROR: [colum {}] \"Term source REF\" found without Characteristics or Unit".format(i+1))
                    # Add Term Accession Number
                    elif "termaccessionnumber" in header_dict and i in header_dict["termaccessionnumber"]:
                        if "termsourceref" in header_dict and i-1 in header_dict["termsourceref"] and last_termsource:
                            if "characteristics" in header_dict and i-2 in header_dict["characteristics"] and last_attribute:
                                sample_attributes["characteristics"][last_attribute]["term accession"] = row[i]
                            elif "unit" in header_dict and i-2 in header_dict["unit"] and last_attribute and last_unit:
                                sample_attributes["characteristics"][last_attribute]["unit"]["term accession"] = row[i]
                        else:
                            logger.log(40, "PARSER ERROR: [column {}] \"Term Accession Number\" found without Term Source REF".format(i+1))
                    # Get Material Type
                    elif "materialtype" in header_dict and i in header_dict["materialtype"]:
                        sample_attributes["material type"] = row[i]
                    # Get Description
                    elif "description" in header_dict and i in header_dict["description"]:
                        sample_attributes["description"] = row[i]
                    # Get Comments
                    sample_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1])

                samples[sample_name] = sample_attributes

        else:
            # We don't want to log this for each row. Maybe make an error dict and count the occurances and print at the end of the script.
            # Or simpler test this first before going through the rows
            logger.log(40, "PARSER ERROR: No \"Source Name\" column found, or more than one.")



        # Extract

        extract_attributes = {"name": "",
                              "comments": {},
                              "sample ref": "",
                              "protocol ref": []}

        if "extractname" in header_dict:
            # This creates separate extract dicts if there is more than one column. Should we allow that?
            node_range = node_map["extractname"][0]

            # Expecting the first value to the be node name
            extract_name = row[node_range[0]]

            if extract_name not in extract_names:
                extract_names.append(extract_name)
                extract_attributes["name"] = extract_name
                # Get comments
                extract_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1] + 1)
                # Keep reference to sample in that row
                extract_attributes["sample ref"] = sample_name

                extract_attributes["protocol ref"] = [row[i] for i in node_range[2]]
                extracts[extract_name] = extract_attributes



        # Labeled Extract

        le_attributes = {"name": "",
                         "label": "",
                         "comments": {},
                         "extract ref": "",
                         "protocol ref": []}

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
                le_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1] + 1)
                # Keep reference to sample in that row
                le_attributes["extract ref"] = extract_name

                le_attributes["protocol ref"] = [row[i] for i in node_range[2]]
                le[le_name] = le_attributes


        # Assay

        assay_attributes = {"name": "",
                            "technology type": "",
                            "comments": {},
                            "extract ref": set(),
                            "protocol ref": []}

        if "assayname" in header_dict:
            node_range = node_map["assayname"][0]
            assay_name = row[node_range[0]]

            if assay_name not in assay_names:
                assay_attributes["name"] = assay_name
                assay_names.append(assay_name)

                for i in range(node_range[0]+1, node_range[1]+1):
                    if get_value(header[i]) == "technologytype":
                        assay_attributes["technology type"] = row[i]
                assay_attributes["comments"] = get_comment_values(row, header, node_range[0], node_range[1])
                assay_attributes["protocol ref"] = [row[i] for i in node_range[2]]
                assays[assay_name] = assay_attributes

                if le_name:
                    assay_attributes["extract ref"].add(le_name)
                    continue
                elif extract_name:
                    assay_attributes["extract ref"].add(extract_name)
            # We need to allow more than 1 extract per assay for the two-colour case
            else:
                # Add the extract ref to the previous assay, using set to avoid duplicates
                if le_name:
                    assays[assay_name]["extract ref"].add(le_name)
                    continue
                elif extract_name:
                    assays[assay_name]["extract ref"].add(extract_name)


        # Datafiles

        raw_data = defaultdict()
        datatype = None
        files = []

        for rdn in raw_data_nodes:
            if rdn in header:
                files = [row[x] for x in header[rdn]]
                datatype = rdn

        for f in files:
            # Make entry in the raw data dict under the file prefix
            raw_data[get_prefix(f)].append({
                'filename': f,
                'datatype': datatype,
                'assayref': '',
                'comments': get_comment_values(row, header, header_dict[datatype], header_dict.get('derivedarraydatafile', len(header)-1))
            })

    return samples, extracts, le, assays, raw_data, processed_data