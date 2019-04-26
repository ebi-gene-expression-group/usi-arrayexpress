import unittest
from converter.datamodel2magetab import get_protocol_positions


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
        assert protocol_positions[2] == ["nucleic acid sequencing protocol"]

    def test_single_cell_protocol(self):
        techtype = "singlecell"
        protocol_positions = get_protocol_positions(techtype)
        print(protocol_positions)
        assert protocol_positions[2] == ["nucleic acid sequencing protocol"]
