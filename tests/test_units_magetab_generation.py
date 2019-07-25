import unittest
from converter.dm2magetab import get_protocol_positions, sort_protocol_refs_to_dict, flatten_sample_attribute
from converter.datamodel.components import Attribute, Unit
from converter.datamodel.protocol import Protocol


class TestReadingProtocolPositions(unittest.TestCase):

    def test_microarray_protocols(self):
        techtype = "microarray"
        protocol_positions = get_protocol_positions(techtype)
        print(protocol_positions)
        assert(protocol_positions[2] == ["nucleic acid labeling protocol"])

    def test_sequencing_protocol(self):
        techtype = "sequencing"
        protocol_positions = get_protocol_positions(techtype)
        print(protocol_positions)
        assert(protocol_positions[4] == ["nucleic acid sequencing protocol"])

    def test_single_cell_protocol(self):
        techtype = "singlecell"
        protocol_positions = get_protocol_positions(techtype)
        print(protocol_positions)
        assert(protocol_positions[4] == ["nucleic acid sequencing protocol"])


class TestTransformingProtocolRefs(unittest.TestCase):

    def setUp(self):
        self.protocol_positions = get_protocol_positions("microarray")
        self.all_protocols = [Protocol(alias="Protocol 1", accession="P-MTAB-1234", description="testing col",
                                       protocol_type=Attribute(value="sample collection protocol")),
                              Protocol(alias="Protocol 2", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="growth protocol")),
                              Protocol(alias="Protocol 3", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="growth protocol")),
                              Protocol(alias="Protocol 4", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="treatment protocol")),
                              Protocol(alias="Protocol 5", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="conversion protocol")),
                              Protocol(alias="Protocol 6", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="nucleic acid extraction protocol")),
                              Protocol(alias="Protocol 7", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="nucleic acid labeling protocol")),
                              Protocol(alias="Protocol 8", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="nucleic acid hybridization to array protocol")),
                              Protocol(alias="Protocol 9", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="array scanning and feature extraction protocol")),
                              Protocol(alias="Protocol 10", accession="P-MTAB-1235", description="testing scan",
                                       protocol_type=Attribute(value="normalization data transformation protocol"))]

    def test_microarray_protocol_refs(self):
        for i in range(5):
            protocol_refs = sort_protocol_refs_to_dict(self.protocol_positions, self.all_protocols)
            print(protocol_refs)
            assert(protocol_refs == {1: {'11~~~Protocol REF': 'Protocol 1',
                                         '12~~~Protocol REF': 'Protocol 2',
                                         '13~~~Protocol REF': 'Protocol 3',
                                         '14~~~Protocol REF': 'Protocol 4',
                                         '15~~~Protocol REF': 'Protocol 5',
                                         '16~~~Protocol REF': 'Protocol 6'},
                                     2: {'21~~~Protocol REF': 'Protocol 7'},
                                     3: {'31~~~Protocol REF': 'Protocol 8'},
                                     5: {'51~~~Protocol REF': 'Protocol 9'},
                                     6: {'61~~~Protocol REF': 'Protocol 10'}})


class TestFlatteningAttributesToList(unittest.TestCase):

    def setUp(self):
        self.test_age_category = "age"
        self.test_age_attribute = Attribute(value="24", unit=Unit(value="year", unit_type="time unit"))
        self.test_dose_category = "dose"
        self.test_dose_attribute = Attribute(value="50", unit=Unit(value="micromolar", unit_type="concentration unit",
                                                                   term_accession="EFO_000123", term_source="EFO"))
        self.test_organism_category = "organism"
        self.test_organism_attribute = Attribute(value="Homo sapiens",
                                                 term_accession="NCBITaxon_9606", term_source="NCBITAXON")
        self.test_sex_category = "sex"
        self.test_sex_attribute = Attribute(value="female")

    def test_flatten_attribute_with_unit(self):
        result = flatten_sample_attribute(self.test_age_category, self.test_age_attribute, "Characteristics")
        assert(result == [("Characteristics[age]", "24"),
                          ("age~~~Unit[time unit]", "year")])

    def test_flatten_factor(self):
        result = flatten_sample_attribute(self.test_dose_category, self.test_dose_attribute, "Factor Value")
        assert(result == [("Factor Value[dose]", "50"),
                          ("dose~~~Unit[concentration unit]", "micromolar"),
                          ("dose-unit~~~Term Source REF", "EFO"),
                          ("dose-unit~~~Term Accession Number", "EFO_000123")])

    def test_flatten_attribute_with_terms(self):
        result = flatten_sample_attribute(self.test_organism_category, self.test_organism_attribute, "Characteristics")
        assert(result == [("Characteristics[organism]", "Homo sapiens"),
                          ("organism~~~Term Source REF", "NCBITAXON"),
                          ("organism~~~Term Accession Number", "NCBITaxon_9606")])

    def test_flatten_simple_attribute(self):
        result = flatten_sample_attribute(self.test_sex_category, self.test_sex_attribute, "Characteristics")
        assert(result == [("Characteristics[sex]", "female")])

    def test_empty_object(self):
        result = flatten_sample_attribute("", Attribute(value=None), "Characteristics")
        assert(result == [])
