import unittest
from converter.datamodel2magetab import get_protocol_positions, sort_protocol_refs_to_dict
from converter.datamodel import Protocol, Attribute


class TestReadingProtocolPositions(unittest.TestCase):

    def test_microarray_protocols(self):
        techtype = "microarray"
        protocol_positions = get_protocol_positions(techtype)
        print(protocol_positions)
        assert protocol_positions[2] == ["nucleic acid labeling protocol"]

    def test_sequencing_protocol(self):
        techtype = "sequencing"
        protocol_positions = get_protocol_positions(techtype)
        print(protocol_positions)
        assert protocol_positions[4] == ["nucleic acid sequencing protocol"]

    def test_single_cell_protocol(self):
        techtype = "singlecell"
        protocol_positions = get_protocol_positions(techtype)
        print(protocol_positions)
        assert protocol_positions[4] == ["nucleic acid sequencing protocol"]


class TestTransformingProtocolRefs(unittest.TestCase):

    def setUp(self):
        self.protocol_positions = get_protocol_positions("microarray")
        self.all_protocols = [Protocol("Protocol 1", "P-MTAB-1234", "testing col",
                                       Attribute("sample collection protocol", None, None), None, None, None),
                              Protocol("Protocol 2", "P-MTAB-1235", "testing scan",
                                       Attribute("growth protocol", None, None), None, None, None),
                              Protocol("Protocol 3", "P-MTAB-1235", "testing scan",
                                       Attribute("growth protocol", None, None), None, None, None),
                              Protocol("Protocol 4", "P-MTAB-1235", "testing scan",
                                       Attribute("treatment protocol", None, None), None, None, None),
                              Protocol("Protocol 5", "P-MTAB-1235", "testing scan",
                                       Attribute("conversion protocol", None, None), None, None, None),
                              Protocol("Protocol 6", "P-MTAB-1235", "testing scan",
                                       Attribute("nucleic acid extraction protocol", None, None), None, None, None),
                              Protocol("Protocol 7", "P-MTAB-1235", "testing scan",
                                       Attribute("nucleic acid labeling protocol", None, None), None, None, None),
                              Protocol("Protocol 8", "P-MTAB-1235", "testing scan",
                                       Attribute("nucleic acid hybridization to array protocol", None, None), None, None, None),
                              Protocol("Protocol 9", "P-MTAB-1235", "testing scan",
                                       Attribute("array scanning and feature extraction protocol", None, None), None, None, None),
                              Protocol("Protocol 10", "P-MTAB-1235", "testing scan",
                                       Attribute("normalization data transformation protocol", None, None), None, None, None)
                              ]

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
