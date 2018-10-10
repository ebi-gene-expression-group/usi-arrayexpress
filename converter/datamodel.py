
class Assay:

    def __init__(self, assay_obj):
        # Attributes common to all assay types
        self.name = assay_obj.get("alias", None)
        self.attributes = assay_obj.get("attributes", None)
        self.techtype = self.attributes.get("Technology Type", None)
        self.protocolrefs = assay_obj.get("protocolUses")
        self.sampleUses = assay_obj.get("sampleUses", None)
        self.comments = self.attributes.get("comments", None)
        #self.factors = factors


    def get_techtype(self):
        try:
            techtype_value = self.techtype[0]["value"]
            return techtype_value
        except KeyError:
            return None


class MicroarrayAssay(Assay):
    def __init__(self, assay_obj):
        Assay.__init__(assay_obj)
        self.arrayref = self.attributes.get("array ref", None)


class SeqAssay(Assay):
    def __init__(self, title, techtype, factors, protocolrefs, samplesref, lib_attribs):
        super().__init__(techtype, factors, protocolrefs, samplesref)
        self.title = title
        self.lib_layout = lib_attribs.get('library_layout', None)
        self.lib_selection = lib_attribs.get('library_selection', None)
        self.lib_strategy = lib_attribs.get('library_strategy', None)
        self.lib_strand = lib_attribs.get('library_strand', None)
        self.lib_source = lib_attribs.get('library_source', None)
        self.lib_orientation = lib_attribs.get('orientation', None)
        self.lib_nomlen = lib_attribs.get('nominal_length', None)
        self.lib_nomsd = lib_attribs.get('nominal_sdev', None)


    @classmethod
    def from_sdrf(cls, extract_attributes, assay_attributes):
        """Intialise assay attributes from sdrf data dicts"""