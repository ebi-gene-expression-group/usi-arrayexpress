from converter.datamodel.submittable import DependentSubmittable
from utils.converter_utils import remove_duplicates, get_controlled_vocabulary


class Assay(DependentSubmittable):
    """
    Attributes common to all assay types
    :param alias: string, unique name of the assay
    :param technology_type: string, one of: array assay, sequencing assay
    :param protocolrefs: list, protocol accessions/name used to generate assay data
    :param sampleref: string, sample accession/name
    """

    def __init__(self, **kwargs):
        DependentSubmittable.__init__(self, **kwargs)
        self.technology_type = kwargs.get("technology_type")
        self.sampleref = kwargs.get("sampleref")

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, technology_type={self.technology_type}, " \
               "protocolrefs={self.protocolrefs}, sampleref={self.sampleref})".format(self=self)


class MicroarrayAssay(Assay):
    """
    Assay attributes specific to microarray assays
    :param label: string, the label type
    :param array_design: string, ArrayExpress array design format accession
    """

    def __init__(self, **kwargs):
        Assay.__init__(self, **kwargs)
        self.label = kwargs.get("label")
        self.array_design = kwargs.get("array_design")

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, technology_type={self.technology_type}, " \
               "protocolrefs={self.protocolrefs}, sampleref={self.sampleref}, " \
               "label={self.label}, array_design={self.array_design})".format(self=self)

    @classmethod
    def from_magetab(cls, le_attributes, extract_attributes, assay_attributes):
        """Intialise assay attributes from MAGE-TAB data dicts.
        The central node for the microarray assay is the Labeled Extract Name."""

        technology_type = remove_duplicates([a.get("technology_type") for a in assay_attributes])

        # Get all protocol refs
        protocolrefs = []
        for a in assay_attributes:
            protocolrefs.extend(a.get("protocol_ref"))
        protocolrefs.extend(extract_attributes.get("protocol_ref", []))
        protocolrefs.extend(le_attributes.get("protocol_ref", []))
        protocolrefs = remove_duplicates(protocolrefs)

        # Get Array design REF, we are only expecting one unique per extract
        array_design = [a.get("array_design") for a in assay_attributes]

        return cls(alias=le_attributes.get("name"),
                   technology_type=technology_type[0],
                   protocolrefs=protocolrefs,
                   sampleref=extract_attributes.get("sample_ref"),
                   label=le_attributes.get("label"),
                   array_design=array_design[0])


class SeqAssay(Assay):
    """
    Assay attributes specific to sequencing assays
    :param accession: string, (not for microarray) ENA experiment accession
    :param lib_attribs: dictionary of sequencing specific attributes that are required for ENA brokering
                        keys: library attribute
                        values: values from SRA's controlled vocabulary
    """
    def __init__(self, **kwargs):

        Assay.__init__(self, **kwargs)
        self.accession = kwargs.get("accession")
        self.library_layout = kwargs.get("library_layout")
        self.library_selection = kwargs.get("library_selection")
        self.library_strategy = kwargs.get("library_strategy")
        self.library_strand = kwargs.get("library_strand")
        self.library_source = kwargs.get("library_source")
        self.orientation = kwargs.get("orientation")
        self.nominal_length = kwargs.get("nominal_length")
        self.nominal_sdev = kwargs.get("nominal_sdev")
        self.platform_type = kwargs.get("platform_type")
        self.instrument_model = kwargs.get("instrument_model")

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, technology_type={self.technology_type}, " \
               "protocolrefs={self.protocolrefs}, sampleref={self.sampleref} " \
               "accession={self.accession}, library_layout={self.library_layout}, " \
               "library_selection={self.library_selection}, library_strategy={self.library_strategy}, " \
               "library_strand={self.library_strand}, library_source={self.library_source}," \
               "orientation={self.orientation}, nominal_length={self.nominal_length}, " \
               "nominal_sdev={self.nominal_sdev}, platform_type={self.platform_type}, " \
               "instrument_model={self.instrument_model})".format(self=self)

    @classmethod
    def from_magetab(cls, extract_attributes, assay_attributes, protocols):
        """Intialise assay attributes from MAGE-TAB data dicts.
        The central node for the sequencing assay is the Extract Name."""

        # Get library attributes from extract comments
        comments = extract_attributes.get("comments")
        lib_attrib_cv = get_controlled_vocabulary("sdrf_comments_ena").keys()
        lib_attribs = {a.lower(): comments[a] for a in lib_attrib_cv if comments.get(a)}

        # Get technology type(s) from assay attributes
        technology_type = remove_duplicates([a.get("technology_type", "") for a in assay_attributes])
        if len(technology_type) > 0:
            technology_type = technology_type[0]

        # Get accession from ENA Experiment in assay comments
        accession = remove_duplicates([a.get('comments', {}).get('ENA_EXPERIMENT', "") for a in assay_attributes])
        if len(accession) == 1:
            accession = accession[0]
        #else: report ERROR

        # Get all protocol refs
        protocolrefs = []
        for a in assay_attributes:
            protocolrefs.extend(a.get("protocol_ref"))
        protocolrefs.extend(extract_attributes.get("protocol_ref", []))
        protocolrefs = remove_duplicates(protocolrefs)

        # Get platform and instrument from sequencing protocol
        for p in protocols:
            if p.get("title") in protocolrefs and p.get("protocol_type") == "nucleic acid sequencing protocol":
                hardware = p.get("hardware")
                lib_attribs["instrument_model"] = hardware
                # TODO: need to look up the platform type from ENA's controlled vocab
                # Using ILLUMINA as placeholder for now as it fits 90% of cases
                lib_attribs["platform_type"] = "ILLUMINA"
                # Expecting only one sequencing protocol per assay
                break

        return cls(alias=extract_attributes.get("name"),
                   accession=accession,
                   technology_type=technology_type,
                   protocolrefs=protocolrefs,
                   sampleref=extract_attributes.get("sample_ref"),
                   **lib_attribs)


class SingleCellAssay(SeqAssay):
    """Assay attributes specific to single-cell sequencing assays"""

    def __init__(self, **kwargs):
        SeqAssay.__init__(self, **kwargs)

        self.single_cell_isolation = kwargs.get("single_cell_isolation")
        self.library_construction = kwargs.get("library_construction")
        self.input_molecule = kwargs.get("input_molecule")
        self.end_bias = kwargs.get("end_bias")
        self.primer = kwargs.get("primer")
        self.spike_in = kwargs.get("spike_in")
        self.spike_in_dilution = kwargs.get("spike_in_dilution")
        self.umi_barcode_read = kwargs.get("umi_barcode_read")
        self.umi_barcode_size = kwargs.get("umi_barcode_size")
        self.umi_barcode_offset = kwargs.get("umi_barcode_offset")
        self.cell_barcode_read = kwargs.get("cell_barcode_read")
        self.cell_barcode_size = kwargs.get("cell_barcode_size")
        self.cell_barcode_offset = kwargs.get("cell_barcode_offset")
        self.cDNA_read = kwargs.get("cDNA_read")
        self.cDNA_read_size = kwargs.get("cDNA_read_size")
        self.cDNA_read_offset = kwargs.get("cDNA_read_offset")
        self.sample_barcode_read = kwargs.get("sample_barcode_read")
        self.sample_barcode_size = kwargs.get("sample_barcode_size")
        self.sample_barcode_offset = kwargs.get("sample_barcode_offset")

