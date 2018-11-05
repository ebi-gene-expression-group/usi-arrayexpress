from parsing import read_json_file, get_controlled_vocabulary, remove_duplicates


class Assay:

    def __init__(self, alias, accession, techtype, protocolrefs, sampleref):
        # Attributes common to all assay types
        self.alias = alias
        self.accession = accession
        self.techtype = techtype
        self.protocolrefs = protocolrefs
        self.sampleref = sampleref

    def get_attributes(self):
        """This returns a list of all additional attributes that have values"""
        all_attributes = self.__dict__.keys()
        all_values = self.__dict__
        fixed_attributes = ("alias", "accession", "techtype", "protocolrefs", "sampleref")
        other_attributes = list()
        for a in all_attributes:
            if a not in fixed_attributes and all_values[a]:
                other_attributes.append(a)
        return other_attributes


class MicroarrayAssay(Assay):
    def __init__(self, assay_obj):
        Assay.__init__(self, assay_obj)
        self.arrayref = self.attributes.get("array ref", None)


class SeqAssay(Assay):
    def __init__(self, alias, accession, techtype, protocolrefs, sampleref, lib_attribs):
        Assay.__init__(self, alias, accession, techtype, protocolrefs, sampleref)

        self.library_layout = lib_attribs.get("library_layout")
        self.library_selection = lib_attribs.get("library_selection")
        self.library_strategy = lib_attribs.get("library_strategy")
        self.library_strand = lib_attribs.get("library_strand")
        self.library_source = lib_attribs.get("library_source")
        self.orientation = lib_attribs.get("orientation")
        self.nominal_length = lib_attribs.get("nominal_length")
        self.nominal_sdev = lib_attribs.get("nominal_sdev")


    @classmethod
    def from_sdrf(cls, extract_attributes, assay_attributes):
        """Intialise assay attributes from sdrf data dicts,
        Translating parsed SDRF attributes into USI format"""

        alias = extract_attributes.get("name")

        # Get library attributes from extract comments
        comments = extract_attributes.get("comments")
        lib_attrib_cv = get_controlled_vocabulary("sdrf_comments_ena").keys()
        lib_attribs = {a.lower(): comments[a] for a in lib_attrib_cv if comments.get(a)}

        # Get technology type(s) from assay attributes
        techtype = remove_duplicates([a.get("technology type") for a in assay_attributes])

        # Get accession from ENA Experiment in assay comments
        accession = remove_duplicates([a.get('comments', {}).get('ENA_EXPERIMENT', "") for a in assay_attributes])
        if len(accession) == 1:
            accession = accession[0]
        #else: report ERROR

        # Get all protocol refs
        protocolrefs = []
        for a in assay_attributes:
            protocolrefs.extend(a.get("protocol ref"))
        protocolrefs.extend(extract_attributes.get("protocol ref", []))
        protocolrefs = remove_duplicates(protocolrefs)

        # Get sample refs
        samplerefs = extract_attributes.get("sample ref")

        return cls(alias, accession, techtype, protocolrefs, samplerefs, lib_attribs)

    @classmethod
    def from_json(cls, assay_object):

        pass


class Protocol:
    def __init__(self, alias, title, description, protocol_type, hardware, software, parameters):
        self.alias = alias
        self.title = title
        self.description = description
        self.protocol_type = protocol_type
        self.hardware = hardware
        self.software = software
        self.parameters = parameters

    @classmethod
    def from_idf(cls, protocol_dict):
        title = protocol_dict.get("title")
        alias = title
        description = protocol_dict.get("description")
        protocol_type = protocol_dict.get("protocol type")
        hardware = protocol_dict.get("hardware")
        software = protocol_dict.get("software")
        parameters = protocol_dict.get("parameters")

        return cls(alias, title, description, protocol_type, hardware, software, parameters)


class Study:
    def __init__(self, accession, alias, title, description, protocolrefs, projectref,
                 experimental_factor, experimental_design, experiment_type, comments):
        self.accession = accession
        self.alias = alias
        self.title = title
        self.description = description
        self.protocolrefs = protocolrefs
        self.projectref = projectref
        self.experimental_factor = experimental_factor
        self.experimental_design = experimental_design
        self.experiment_type = experiment_type
        self.comments = comments

    @classmethod
    def from_idf(cls, study_info):
        accession = study_info.get("accession")
        alias = study_info.get("title")
        title = study_info.get("title")
        description = study_info.get("description")
        experiment_type = study_info.get("experiment type")
        experimental_factor = study_info.get("experimental factor")
        experimental_design = study_info.get("experimental design")
        protocolrefs = study_info.get("protocolRefs")
        projectref = "project_" + alias
        comments = study_info.get("comments")

        return cls(accession, alias, title, description, protocolrefs, projectref,
                   experimental_factor, experimental_design, experiment_type, comments)
