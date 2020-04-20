import unittest


class UtilitiesTest(unittest.TestCase):
    def test_replace_ambiguous_bases(self):
        from src.probe_qpcr import replace_ambiguous_bases
        sequence = 'AGYTCGY'
        expected = 'AG[CT]TCG[CT]'
        actual = replace_ambiguous_bases(sequence)
        self.assertEqual(expected, actual)

        sequence = 'ACZTG'
        self.assertRaises(KeyError, replace_ambiguous_bases, sequence)

    def test_get_qpcr_hits(self):
        primer_dict = {"Test_Primer": {
            'F': '',
            'R': '',
            'P': ''
        }}
        expected = "something"
        actual = "placeholder"

    def test_reverse_complement(self):
        from src.probe_qpcr import reverse_complement
        sequence_1 = "ACGT"
        expected_1 = 'ACGT'
        actual_1 = reverse_complement(sequence_1)
        self.assertEqual(expected_1, actual_1)

        sequence_2 = 'ATCGTGCTGCTGTCGTCAAGAC'
        expected_2 = 'GTCTTGACGACAGCAGCACGAT'
        actual_2 = reverse_complement(sequence_2)
        self.assertEqual(expected_2, actual_2)

        sequence_3 = 'TGCTAGCATCGAGTCGATCGATATATTTAGCATCAGCATT'
        expected_3 = 'AATGCTGATGCTAAATATATCGATCGACTCGATGCTAGCA'
        actual_3 = reverse_complement(sequence_3)
        self.assertEqual(expected_3, actual_3)

        sequence_4 = 'AG[CT]TCG[CT]'
        expected_4 = '[AG]CGA[AG]CT'
        actual_4 = reverse_complement(sequence_4)
        self.assertEqual(expected_4, actual_4)


if __name__ == '__main__':
    unittest.main(verbosity=2)
