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


if __name__ == '__main__':
    unittest.main(verbosity=2)
