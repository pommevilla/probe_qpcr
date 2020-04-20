from io import StringIO
import unittest
from unittest.mock import patch


class UtilitiesTest(unittest.TestCase):
    def test_replace_ambiguous_bases(self):
        from src.probe_qpcr import replace_ambiguous_bases
        sequence = 'AGYTCGY'
        expected = 'AG[CT]TCG[CT]'
        actual = replace_ambiguous_bases(sequence)
        self.assertEqual(expected, actual)

        sequence = 'ACZTG'
        self.assertRaises(KeyError, replace_ambiguous_bases, sequence)

    def test_read_primers(self):
        from src.probe_qpcr import read_primers, reverse_complement, replace_ambiguous_bases
        expected = {
            'China.ORF1ab': {
                'F': 'CCCTGTGGGTTTTACACTTAA',
                'R': 'ACGATTGTGCATCAGCTGA',
                'P': 'CCGTCTGCGGTATGTGGAAAGGTTATGG',
            },
            'USA.2019-nCoV_N3': {
                'F': 'GGGAGCCTTGAATACACCAAAA',
                'R': 'TGTAGCACGATTGCAGCATTG',
                'P': 'A[CT]CACATTGGCACCCGCAATCCTG',
            },
        }

        actual = read_primers('data/test_primers.fa')
        self.assertEqual(expected, actual)

    def test_get_qpcr_hits(self):
        from src.probe_qpcr import get_qpcr_hits

        output_format_string = '{}\t{}\t{}\t{}\t{}\t{}\n'

        primers = {
            'test_primer.1': {
                'F': 'GTACA',
                'R': 'TTGCGG',
                'P': 'TCCTG',
            },
            'test_primer.2': {
                'F': 'TGCTGAAA',
                'R': 'ACGAGG[CGT]',
                'P': 'CCT[CT]TG',
            }
        }

        test_sequences = {
            'test_sequence.1': 'TTTAAGCTGTACAGCCCTCCTGAATGCTGTCCCCGCAATTCGA',
            'test_sequence.2': 'AACTGCTGAAACCCTGTATCGTAGCCCTCTGTTAAGCCTAACCTCGTAGGCTA',
        }

        expected_1 = output_format_string.format('test_primer.1', 'test_sequence.1',
                                                 8, 38, 30,
                                                 'GTACAGCCCTCCTGAATGCTGTCCCCGCAA')

        with patch('sys.stdout', new=StringIO()) as captured_output:
            get_qpcr_hits('test_primer.1', primers['test_primer.1']['F'],
                          primers['test_primer.1']['R'],
                          primers['test_primer.1']['P'],
                          'test_sequence.1',
                          test_sequences['test_sequence.1'],
                          )
            self.assertEqual(expected_1, captured_output.getvalue())

        expected_2 = ''

        with patch('sys.stdout', new=StringIO()) as captured_output:
            get_qpcr_hits('test_primer.2', primers['test_primer.2']['F'],
                          primers['test_primer.2']['R'],
                          primers['test_primer.2']['P'],
                          'test_sequence.1',
                          test_sequences['test_sequence.1'],
                          )

            self.assertEqual(expected_2, captured_output.getvalue())

        expected_3 = output_format_string.format('test_primer.2', 'test_sequence.2',
                                                 3, 47, 44,
                                                 'TGCTGAAACCCTGTATCGTAGCCCTCTGTTAAGCCTAACCTCGT')

        with patch('sys.stdout', new=StringIO()) as captured_output:
            get_qpcr_hits('test_primer.2', primers['test_primer.2']['F'],
                          primers['test_primer.2']['R'],
                          primers['test_primer.2']['P'],
                          'test_sequence.2',
                          test_sequences['test_sequence.2'],
                          )

            self.assertEqual(expected_3, captured_output.getvalue())



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
