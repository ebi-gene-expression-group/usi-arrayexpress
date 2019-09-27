"""The classes for the different assay types"""

from converter.datamodel.submittable import DependentSubmittable


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

    def get_assay_attributes(self):
        """A list of all attributes that are specific to the assay object (not in the general assay class)"""
        exclude = ('id', 'alias', 'accession', 'protocolrefs', 'technology_type', 'sampleref')
        return [at for at in self.get_all_attributes() if at not in exclude]

    def get_attributes_with_values(self):
        """A list of all assay attributes that have values"""
        exclude = ('id', 'alias', 'accession', 'protocolrefs', 'sampleref')
        return [at for at in self.get_all_attributes() if getattr(self, at) and at not in exclude]


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
        self.hca_bundle_uuid = kwargs.get("hca_bundle_uuid")
        self.hca_bundle_version = kwargs.get("hca_bundle_version")

    def get_singlecell_attributes(self, invert=False):
        """A list of single-cell specific attributes that are not in the parent sequencing assay class,
        with invert=True returns the sequencing assay attributes."""
        dummy_assay = SeqAssay()
        exclude = dummy_assay.get_all_attributes()
        if invert:
            return dummy_assay.get_assay_attributes()
        return [at for at in self.get_all_attributes() if at not in exclude]
