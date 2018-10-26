from parsing import read_json_file, get_controlled_vocabulary


class Assay:

    def __init__(self, alias, accession, techtype, protocolrefs, samplerefs):
        # Attributes common to all assay types
        self.alias = alias
        self.accession = accession
        self.techtype = techtype
        self.protocolrefs = protocolrefs
        self.samplerefs = samplerefs

    def get_attributes(self):
        all_attributes = self.__dict__.keys()
        return all_attributes

    def get_techtype(self):
        try:
            techtype_value = self.techtype[0]["value"]
            return techtype_value
        except KeyError:
            return None


    def generate_assay_json(self, study_info):

        assay_json = {
            "alias": self.alias,
            "attributes": {
                "technology_type": self.techtype
            },
            "studyRef": "",
            "sampleUses": "",
            "protocolUses": "",
            "team": study_info["team"]
        }
        return assay_json


class MicroarrayAssay(Assay):
    def __init__(self, assay_obj):
        Assay.__init__(self, assay_obj)
        self.arrayref = self.attributes.get("array ref", None)


class SeqAssay(Assay):
    def __init__(self, alias, accession, techtype, protocolrefs, samplerefs, lib_attribs):
        Assay.__init__(self, alias, accession, techtype, protocolrefs, samplerefs)

        self.library_layout = lib_attribs.get("library_layout")
        self.library_selection = lib_attribs.get("library_selection")
        self.library_strategy = lib_attribs.get("library_strategy")
        self.library_strand = lib_attribs.get("library_strand")
        self.library_source = lib_attribs.get("library_source")
        self.library_orientation = lib_attribs.get("orientation")
        self.library_nomlen = lib_attribs.get("nominal_length")
        self.library_nomsd = lib_attribs.get("nominal_sdev")


    @classmethod
    def from_sdrf(cls, extract_attributes, assay_attributes):
        """Intialise assay attributes from sdrf data dicts,
        Translating parsed SDRF attributes into USI format"""

        alias = extract_attributes.get("name")
        techtype = assay_attributes.get("technology type")

        # Get library attributes from extract comments
        comments = extract_attributes.get("comments")
        cv = read_json_file("../converter/expected_terms.json")
        lib_attrib_cv = cv["library_terms"]
        # Transform into USI layout
        lib_attribs = {a.lower(): comments[a] for a in lib_attrib_cv if comments.get(a)}

        # Get accession from ENA Experiment in assay comments
        assay_comments = assay_attributes.get('comments')
        if assay_comments:
            accession = assay_comments.get('ENA_EXPERIMENT')
        else:
            accession = None

        # Get all protocol refs
        protocolrefs = extract_attributes.get("protocol ref", [])
        # Also need to add assay protocol refs from assay node
        for p in assay_attributes.get("protocol ref", []):
            protocolrefs.append(p)

        # Get sample refs
        samplerefs = extract_attributes.get("sample ref", None)

        return cls(alias, accession, techtype, protocolrefs, samplerefs, lib_attribs)


    @classmethod
    def from_json(cls, assay_object):

        pass


