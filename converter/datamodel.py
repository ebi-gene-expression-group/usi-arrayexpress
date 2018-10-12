from parsing import read_json_file


class Assay:

    def __init__(self, alias, techtype, protocolrefs, samplerefs):
        # Attributes common to all assay types
        self.alias = alias
        self.techtype = techtype
        self.protocolrefs = protocolrefs
        self.samplerefs = samplerefs


    def get_techtype(self):
        try:
            techtype_value = self.techtype[0]["value"]
            return techtype_value
        except KeyError:
            return None


class MicroarrayAssay(Assay):
    def __init__(self, assay_obj):
        Assay.__init__(self, assay_obj)
        self.arrayref = self.attributes.get("array ref", None)


class SeqAssay(Assay):
    def __init__(self, alias, techtype, protocolrefs, samplerefs, lib_attribs):
        Assay.__init__(self, alias, techtype, protocolrefs, samplerefs)
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
        lib_attribs = {a.lower(): {"value": comments[a]} for a in lib_attrib_cv if comments.get(a)}

        protocolUses = extract_attributes.get("protocol refs", None)
        #protocolUses2 = assay_attributes.get("protocol refs")
        #protocolrefs = protocolUses + protocolUses2

        sampleUses = extract_attributes.get("sample refs", None)

        return cls(alias, techtype, protocolUses, sampleUses, lib_attribs)

