
from converter.datamodel.submittable import DependentSubmittable
from converter.datamodel.components import DataFile
from utils.converter_utils import remove_duplicates


class DataSubmittable(DependentSubmittable):
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
        self.assayrefs = kwargs.get("assayrefs")
        self.accession = kwargs.get("accession")

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, files={self.files}, " \
               "data_type={self.data_type}, assayrefs={self.assayrefs}, " \
               "protocolrefs={self.protocolrefs}, accession={self.accession})".format(self=self)

    @classmethod
    def from_magetab(cls, name, datafile_objects, file_attributes):

        alias = name

        # Assuming here that attributes of grouped files are the same
        common_file_attributes = file_attributes[0]

        # Try getting le_ref first for MA
        if common_file_attributes.get("le_ref"):
            assayrefs = remove_duplicates(common_file_attributes.get("le_ref"))
        else:
            assayrefs = remove_duplicates(common_file_attributes.get("extract_ref"))

        mage_data_type = common_file_attributes.get("data_type")
        if mage_data_type == "arraydatafile" or mage_data_type == "scanname":
            data_type = "raw"
        elif mage_data_type == "arraydatamatrixfile":
            data_type = "raw matrix"
        else:
            data_type = None

        comments = common_file_attributes.get("comments")
        accession = comments.get("ENA_RUN")

        protocolrefs = common_file_attributes.get("protocol_ref")

        return cls(alias=alias,
                   files=datafile_objects,
                   data_type=data_type,
                   assayrefs=assayrefs,
                   protocolrefs=protocolrefs,
                   accession=accession)


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
        self.assaydatarefs = kwargs.get("assaydatarefs")

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, files={self.files}, " \
               "data_type={self.data_type}, assaydatarefs={self.assaydatarefs}, " \
               "protocolrefs={self.protocolrefs})".format(self=self)

    @classmethod
    def from_magetab(cls, datafile_objects, file_attributes):

        alias = file_attributes.get("name")

        mage_data_type = file_attributes.get("data_type")
        if mage_data_type == "derivedarraydatamatrixfile":
            data_type = "processed matrix"
        else:
            data_type = "processed"

        protocolrefs = file_attributes.get("protocol_ref")

        # Assay ref is generated based on le_ref (for MA) or extract_ref (for HTS)
        if file_attributes.get("le_ref"):
            assaydatarefs = remove_duplicates(file_attributes.get("le_ref"))
        else:
            assaydatarefs = remove_duplicates(file_attributes.get("extract_ref"))

        return cls(alias=alias,
                   files=datafile_objects,
                   data_type=data_type,
                   assaydatarefs=assaydatarefs,
                   protocolrefs=protocolrefs)


