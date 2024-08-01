#! /usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
import pandas as pd
from unittest.mock import patch, mock_open

# Assuming all your functions are imported from your script
from Summarize_results_module_improved import (generate_html, quit_program, is_tool, usage, count_files, check_folders,
                         extract_number_and_primers, extract_number_from_filename,
                         extract_primer_sequences, get_amplicon_length_from_seq,
                         run_tests, handle_blasts_and_efetch, get_highest_scoring_accession,
                         print_results, generate_results)

class TestPrimerScript(unittest.TestCase):

    def setUp(self):
        # Setup temporary directories and files for testing
        self.test_dir = tempfile.mkdtemp()
        self.primer_file = Path(self.test_dir) / 'primer_123.txt'
        with self.primer_file.open('w') as f:
            f.write("PRIMER_LEFT\nATCG\nPRIMER_RIGHT\nCGTA\nPRIMER_INTERNAL\nGATC\n")

        self.seqkit_file = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m0.txt'
        with self.seqkit_file.open('w') as f:
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT\n")

    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.test_dir)

    def test_generate_html(self):
        header = "Test Header"
        primer_frwd = "ATCG"
        primer_rev = "CGTA"
        primer_internal = "GATC"
        sensi = {"test": "sensi"}
        speci = {"test": "speci"}
        res_target_str = "Test result string"
        primer_fold = Path("/path/to/primers")
        target_fold = Path("/path/to/targets")
        seqkit_fold = Path("/path/to/seqkit")

        html = generate_html(header, primer_frwd, primer_rev, primer_internal, sensi, speci, res_target_str, primer_fold, target_fold, seqkit_fold)
        self.assertIn("<html", html)
        self.assertIn(header.replace(' ', '&nbsp;').replace('\n', '<br>'), html)

    @patch('sys.exit')
    def test_quit_program(self, mock_exit):
        quit_program("Test message")
        mock_exit.assert_called_once_with(1)    

    def test_is_tool(self):
        self.assertTrue(is_tool("cd"))
        self.assertFalse(is_tool("nonexistenttool"))

    def test_usage(self):
        with patch('builtins.print') as mocked_print:
            usage()
            self.assertTrue(mocked_print.called)

    def test_count_files(self):
        path = Path(self.test_dir)
        self.assertEqual(count_files(path), 2)

    def test_check_folders(self):
        check_folders(Path(self.test_dir))  # Should not raise
        with self.assertRaises(SystemExit):
            check_folders(Path("/nonexistent"))

    def test_extract_number_and_primers(self):
        number, pr_frwd, pr_rev, pr_intern, ampli_len = extract_number_and_primers(self.primer_file, Path(self.test_dir))
        self.assertEqual(number, 123)
        self.assertEqual(pr_frwd, "ATCG")
        self.assertEqual(pr_rev, "CGTA")
        self.assertEqual(pr_intern, "GATC")
        self.assertEqual(ampli_len, 7)

    def test_extract_number_from_filename(self):
        self.assertEqual(extract_number_from_filename("primer_123.txt"), 123)
        with self.assertRaises(ValueError):
            extract_number_from_filename("primer.txt")

    def test_extract_primer_sequences(self):
        pr_frwd, pr_rev, pr_intern = extract_primer_sequences(self.primer_file)
        self.assertEqual(pr_frwd, "ATCG")
        self.assertEqual(pr_rev, "CGTA")
        self.assertEqual(pr_intern, "GATC")

    def test_get_amplicon_length_from_seq(self):
        ampli_len = get_amplicon_length_from_seq(self.seqkit_file)
        self.assertEqual(ampli_len, 7)

    @patch('Summarize_results_module_improved.pd.read_csv')
    def test_get_highest_scoring_accession(self, mock_read_csv):
        mock_read_csv.return_value = pd.DataFrame({
            'sseqid': ['acc1', 'acc2'],
            'evalue': [0.01, 0.02],
            'bitscore': [50.0, 40.0]
        })
        sseqid, evalue, bitscore = get_highest_scoring_accession(self.seqkit_file)
        self.assertEqual(sseqid, 'acc1')
        self.assertEqual(evalue, 0.01)
        self.assertEqual(bitscore, 50.0)

    # Additional tests can be added here for other functions like run_tests, handle_blasts_and_efetch, etc.

if __name__ == '__main__':
    unittest.main()
