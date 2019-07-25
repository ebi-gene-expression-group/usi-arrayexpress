from converter.datamodel.components import Attribute
from converter.datamodel.submittable import AccessionedSubmittable
from utils.converter_utils import is_accession


class Protocol(AccessionedSubmittable):
    """
    Attributes of a protocol
    :param alias: string, protocol accession or unique name in experiment
    :param accession: string, ArrayExpress protocol accession
    :param description: string, free-text description of protocol
    :param protocol_type: object, Attribute class object, experiment type ontology term
    :param hardware: string, free-text hardware description, sequencer model for sequencing experiments
    :param software: string, free-text software description
    :param performer: string, protocol performer
    """
    def __init__(self, **kwargs):
        AccessionedSubmittable.__init__(self, **kwargs)
        self.name = kwargs.get("name")
        self.protocol_type = kwargs.get("protocol_type")
        self.hardware = kwargs.get("hardware")
        self.software = kwargs.get("software")
        self.performer = kwargs.get("performer")

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, accession={self.accession}, " \
               "description={self.description}, protocol_type={self.protocol_type}, " \
               "hardware={self.hardware}, software={self.software}, performer={self.performer})".format(self=self)

    def get_ae_attributes(self):
        """Return a list of all AE attributes that have values."""
        all_attributes = self.__dict__.keys()
        all_values = self.__dict__
        ae_attributes = ("protocol_type", "hardware", "software")
        found_attributes = list()
        for a in all_attributes:
            if a in ae_attributes and all_values[a]:
                found_attributes.append(a)
        return found_attributes

    @classmethod
    def from_magetab(cls, protocol_dict):
        alias = protocol_dict.get("title")
        if is_accession(alias):
            accession = alias
        else:
            accession = None
        description = protocol_dict.get("description")
        hardware = protocol_dict.get("hardware")
        software = protocol_dict.get("software")
        protocol_type = Attribute(value=protocol_dict.get("protocol_type"),
                                  unit=None,
                                  term_accession=protocol_dict.get("term_accession"),
                                  term_source=protocol_dict.get("term_source"))

        return cls(alias=alias,
                   accession=accession,
                   description=description,
                   protocol_type=protocol_type,
                   hardware=hardware,
                   software=software)
