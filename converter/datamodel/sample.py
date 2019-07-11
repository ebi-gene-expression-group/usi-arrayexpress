from collections import OrderedDict

from converter.datamodel.submittable import AccessionedSubmittable
from converter.datamodel.components import Attribute, Unit
from utils.converter_utils import get_taxon


class Sample(AccessionedSubmittable):
    """
       :param alias: string, unique sample name in the experiment
       :param accession: string, BioSamples accession
       :param taxon: string, latin species name
       :param taxonId: int, NCBI taxonomy identifier for the species
       :param attributes: dictionary of attribute categories as keys and Attribute class object as value
       :param material_type: string, (optional) one of: whole organism, organism part, cell, RNA, DNA
       :param description: string, (optional) free-text description of sample
       """
    def __init__(self, **kwargs):
        AccessionedSubmittable.__init__(self, **kwargs)
        self.taxon = kwargs.get("taxon")
        self.taxonId = kwargs.get("taxonId")
        self.material_type = kwargs.get("material_type")
        self.attributes = kwargs.get("attributes", {})

        # Set material type if found in attributes
        if not self.material_type and "material_type" in self.attributes:
            self.material_type = self.attributes.get("material_type").value
            del self.attributes["material_type"]

    def __repr__(self):
        return "{self.__class__.__name__}(alias={self.alias}, accession={self.accession}, taxon={self.taxon}, " \
               "taxonId={self.taxonId}, material_type={self.material_type}, description={self.description}, " \
               "attributes={self.attributes})".format(self=self)

    @classmethod
    def from_magetab(cls, sample_attributes):
        alias = sample_attributes.get("name")
        description = sample_attributes.get("description")
        material_type = sample_attributes.get("material_type")

        comments = sample_attributes.get("comments")
        accession = comments.get("BioSD_SAMPLE")

        characteristics = sample_attributes.get("characteristics")
        factors = sample_attributes.get("factors")
        organism = characteristics.get("organism", {})
        taxon = organism.get("value")
        taxonId = get_taxon(taxon)

        # Note this will overwrite the characteristics values if a factor is also a characteristics
        raw_attributes = characteristics.copy()
        raw_attributes.update(factors)

        attributes = OrderedDict()
        for c_name, c_attrib in raw_attributes.items():
            new_unit = None
            if "unit" in c_attrib:
                unit_attrib = c_attrib.get("unit")
                new_unit = Unit(value=unit_attrib.get("value"),
                                unit_type=unit_attrib.get("unit_type"),
                                term_accession=unit_attrib.get("term_accession"),
                                term_source=unit_attrib.get("term_source"))

            attributes[c_name] = Attribute(value=c_attrib.get("value"),
                                           unit=new_unit,
                                           term_accession=c_attrib.get("term_accession"),
                                           term_source=c_attrib.get("term_source"))

        return cls(alias=alias,
                   accession=accession,
                   taxon=taxon,
                   taxonId=taxonId,
                   attributes=attributes,
                   material_type=material_type,
                   description=description)


