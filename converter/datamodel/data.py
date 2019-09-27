"""The classes for data objects"""

from converter.datamodel.submittable import DependentSubmittable
from converter.datamodel.components import DataFile


class DataSubmittable(DependentSubmittable):
    """The base class for all data objects"""

    def __init__(self, **kwargs):
        DependentSubmittable.__init__(self, **kwargs)
        self.files = [DataFile(**d) for d in kwargs.get("files", [])]
        self.data_type = kwargs.get("data_type")


class AssayData(DataSubmittable):
    """
    Attributes of the assay data (raw data object)
    :param alias: string, unique name in the experiment (auto-generated from Assay Name in SDRF)
    :param files: list, DataFile class objects
    :param data_type: string, raw or raw matrix
    :param assayrefs: list, assay name/accessions
    :param protocolrefs: list, protocol name/accessions used to generate data file
    :param accession: string, (not for microarray) ENA run accession
    """

    def __init__(self, **kwargs):
        DataSubmittable.__init__(self, **kwargs)
        self.assayrefs = kwargs.get("assayrefs", [])
        self.accession = kwargs.get("accession")

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, files={self.files}, " \
               "data_type={self.data_type}, assayrefs={self.assayrefs}, " \
               "protocolrefs={self.protocolrefs}, accession={self.accession})".format(self=self)


class Analysis(DataSubmittable):
    """
    Attributes for the analysis (processed data object)
    :param alias: string, unique name in the experiment (auto-generated from processed file name in SDRF)
    :param files: list, DataFile class objects
    :param data_type: string, processed or processed matrix
    :param assaydatarefs: list, assay data name/accessions
    :param protocolrefs: list, protocol name/accessions used to generate data file
    """

    def __init__(self, **kwargs):
        DataSubmittable.__init__(self, **kwargs)
        self.samplerefs = kwargs.get("sampleref")
        self.assayrefs = kwargs.get("assayrefs", [])
        self.assaydatarefs = kwargs.get("assaydatarefs", [])

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, files={self.files}, " \
               "data_type={self.data_type}, assaydatarefs={self.assaydatarefs}, " \
               "protocolrefs={self.protocolrefs})".format(self=self)

