import unittest
from src import probe_qpcr

class UtilitiesTest(unittest.TestCase):
    def test_replace_ambiguous_bases(self):
        sequence = 'AGYTCGY'
        expected = 'AG[CT]TCG[CT]'
        actual = probe_qpcr.replace_ambiguity_codes(sequence)
        self.assertEqual(expected, actual)

if __name__ == '__main__':
    unittest.main(verbosity=3)