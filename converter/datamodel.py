from parsing import read_json_file, get_controlled_vocabulary


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
        techtype = _remove_duplicates([a.get("technology type") for a in assay_attributes])

        # Get accession from ENA Experiment in assay comments
        accession = _remove_duplicates([a.get('comments', {}).get('ENA_EXPERIMENT', "") for a in assay_attributes])
        if len(accession) == 1:
            accession = accession[0]
        #else: report ERROR

        # Get all protocol refs
        protocolrefs = []
        for a in assay_attributes:
            protocolrefs.extend(a.get("protocol ref"))
        protocolrefs.extend(extract_attributes.get("protocol ref", []))
        protocolrefs = _remove_duplicates(protocolrefs)

        # Get sample refs
        samplerefs = extract_attributes.get("sample ref")

        return cls(alias, accession, techtype, protocolrefs, samplerefs, lib_attribs)

    @classmethod
    def from_json(cls, assay_object):

        pass


def _remove_duplicates(ref_list):
    return list(set(ref_list))


