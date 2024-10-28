#!/usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
import subprocess
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock

# Assuming all your functions are imported from your script
from Summarize_results_module_improved import (quit_program, is_tool, usage, count_files, check_folders,
                         extract_number_and_primers, extract_number_from_filename,
                         extract_primer_sequences, get_amplicon_length_from_seq,
                         run_tests, handle_blasts_and_efetch, get_highest_scoring_accession,
                         interpret_and_reformat_sensi_speci_tests)

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
        
        self.seqkit_file1 = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m1.txt'
        with self.seqkit_file1.open('w') as f:
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT\n")
        
        self.seqkit_file2 = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m3.txt'
        with self.seqkit_file2.open('w') as f:
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTGT\n")
        
        self.seqkit_file3 = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_neighbour_m0.txt'
        with self.seqkit_file3.open('w') as f:
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTAC\n")

    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.test_dir)
        

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
        self.assertEqual(count_files(path), 5)

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

class TestHandleBlastsAndEfetch(unittest.TestCase):

    
    @patch('Summarize_results_module_improved.Path.glob')
    def test_no_files_found(self, mock_glob):
        # Setup mock for glob to return an empty list, simulating no files found
        mock_glob.return_value = []

        # Call the function with mocked behavior
        result = handle_blasts_and_efetch(Path('/some/folder'), 1)

        # Assert the expected output
        self.assertEqual(result, "No files for *Target_1.txt_blastx_1e-5.txt found.")

    @patch('Summarize_results_module_improved.get_highest_scoring_accession')
    @patch('Summarize_results_module_improved.Path.glob')
    @patch('Summarize_results_module_improved.Path.stat')
    def test_empty_file(self, mock_stat, mock_glob, mock_get_highest_scoring_accession):
        # Setup mocks for glob and stat to simulate finding an empty file
        mock_file = MagicMock()
        mock_glob.return_value = [mock_file]
        mock_stat.return_value.st_size = 0

        # Setup mock for get_highest_scoring_accession to return "Blast did not find hits"
        mock_get_highest_scoring_accession.return_value = ("Blast did not find hits", float('inf'), 0)

        # Call the function with mocked behavior
        result = handle_blasts_and_efetch(Path('/some/folder'), 1)

        # Assert the expected output
        self.assertEqual(result, "The primer 1's BLASTX search did not find any hits.\n")

    @patch('Summarize_results_module_improved.get_highest_scoring_accession')
    @patch('Summarize_results_module_improved.Path.glob')
    @patch('Summarize_results_module_improved.Path.stat')
    @patch('Summarize_results_module_improved.is_tool')
    @patch('Summarize_results_module_improved.subprocess.Popen')
    @patch('Summarize_results_module_improved.quit_program')
    def test_efetch_available_success(self, mock_quit_program, mock_popen, mock_is_tool, mock_stat, mock_glob, mock_get_highest_scoring_accession):
        # Setup mocks
        mock_is_tool.return_value = True  # Simulate that efetch is available
        
        # Setup mock for glob and stat to simulate finding a file and file being non-empty
        mock_file = MagicMock()
        mock_glob.return_value = [mock_file]
        mock_stat.return_value.st_size = 100
        
        # Setup mock for get_highest_scoring_accession
        mock_get_highest_scoring_accession.return_value = ("acc123", 1e-5, 50)
        
        # Setup mock for subprocess.Popen
        mock_efetch_popen = MagicMock()
        mock_xtract_popen = MagicMock()
        mock_popen.side_effect = [mock_efetch_popen, mock_xtract_popen]
        mock_xtract_popen.communicate.return_value = ("Mock Title", "")
        
        # Call the function
        result = handle_blasts_and_efetch(Path('/some/folder'), 1)
        
        # Assert the expected output
        expected_result = (
            "The primer 1's target with the highest bitscore 50 & evalue 1.00e-05 has the accession acc123 and codes for Mock Title. \nFurther information on the accession: https://www.ncbi.nlm.nih.gov/protein/acc123/ \n"
        )
        self.assertEqual(result, expected_result)
        mock_popen.assert_called()  # Check that subprocess.Popen was called
        mock_quit_program.assert_not_called()  # Check that quit_program was not called

    @patch('Summarize_results_module_improved.get_highest_scoring_accession')
    @patch('Summarize_results_module_improved.Path.glob')
    @patch('Summarize_results_module_improved.Path.stat')
    @patch('Summarize_results_module_improved.is_tool')
    @patch('Summarize_results_module_improved.subprocess.Popen')
    @patch('Summarize_results_module_improved.quit_program')
    def test_efetch_command_failure_xtract(self, mock_quit_program, mock_popen, mock_is_tool, mock_stat, mock_glob, mock_get_highest_scoring_accession):
        # Setup mocks
        mock_is_tool.return_value = True  # Simulate that efetch is available
        
        # Setup mock for glob and stat to simulate finding a file and file being non-empty
        mock_file = MagicMock()
        mock_glob.return_value = [mock_file]
        mock_stat.return_value.st_size = 100
        
        # Setup mock for get_highest_scoring_accession
        mock_get_highest_scoring_accession.return_value = ("acc123", 1e-5, 50)
        
        # Setup mock for subprocess.Popen and its communicate method
        mock_xtract_popen = MagicMock()
        mock_popen.return_value = mock_xtract_popen
        mock_xtract_popen.returncode = 1 
        mock_xtract_popen.communicate.return_value = ("", 'error message')
        
        # After communicate, we raise the CalledProcessError exception
        mock_xtract_popen.wait.side_effect = subprocess.CalledProcessError(1, 'xtract')
        mock_quit_program.side_effect = SystemExit
        # Call the function and assert that quit_program is called correctly
        with self.assertRaises(SystemExit):
            handle_blasts_and_efetch(Path('/some/folder'), 1)
        
        # Verify that quit_program was called once
        self.assertEqual(mock_quit_program.call_count, 1)

    @patch('Summarize_results_module_improved.get_highest_scoring_accession')
    @patch('Summarize_results_module_improved.Path.glob')
    @patch('Summarize_results_module_improved.Path.stat')
    @patch('Summarize_results_module_improved.is_tool')
    @patch('Summarize_results_module_improved.subprocess.Popen')
    @patch('Summarize_results_module_improved.quit_program')
    def test_efetch_command_failure_efetch(self, mock_quit_program, mock_popen, mock_is_tool, mock_stat, mock_glob, mock_get_highest_scoring_accession):
        # Setup mocks
        mock_is_tool.return_value = True  # Simulate that efetch is available
        
        # Setup mock for glob and stat to simulate finding a file and file being non-empty
        mock_file = MagicMock()
        mock_glob.return_value = [mock_file]
        mock_stat.return_value.st_size = 100
        
        # Setup mock for get_highest_scoring_accession
        mock_get_highest_scoring_accession.return_value = ("acc123", 1e-5, 50)
        
        # Setup mock for subprocess.Popen and its communicate method
        mock_efetch_popen = MagicMock()
        mock_popen.return_value = mock_efetch_popen
        mock_efetch_popen.returncode = 1 
        mock_efetch_popen.communicate.return_value = ("", 'error message')
        
        # After communicate, we raise the CalledProcessError exception
        mock_efetch_popen.wait.side_effect = subprocess.CalledProcessError(1, 'xtract')
        mock_quit_program.side_effect = SystemExit
        # Call the function and assert that quit_program is called correctly
        with self.assertRaises(SystemExit):
            handle_blasts_and_efetch(Path('/some/folder'), 1)
        
        # Verify that quit_program was called once
        self.assertEqual(mock_quit_program.call_count, 1)



    @patch('Summarize_results_module_improved.get_highest_scoring_accession')
    @patch('Summarize_results_module_improved.Path.glob')
    @patch('Summarize_results_module_improved.Path.stat')
    @patch('Summarize_results_module_improved.is_tool')
    @patch('Summarize_results_module_improved.subprocess.Popen')
    @patch('Summarize_results_module_improved.quit_program')
    def test_efetch_not_available(self, mock_quit_program, mock_popen, mock_is_tool, mock_stat, mock_glob, mock_get_highest_scoring_accession):
        # Setup mocks
        mock_is_tool.return_value = False  # Simulate that efetch is not available
        
        # Setup mock for glob and stat to simulate finding a file and file being non-empty
        mock_file = MagicMock()
        mock_glob.return_value = [mock_file]
        mock_stat.return_value.st_size = 100
        
        # Setup mock for get_highest_scoring_accession
        mock_get_highest_scoring_accession.return_value = ("acc123", 1e-5, 50)
        
        # Call the function
        result = handle_blasts_and_efetch(Path('/some/folder'), 1)
        
        # Assert the expected output
        expected_result = (
            "The primer 1's target with the highest bitscore 50 & evalue 1.00e-05 has the accession acc123. \nFurther information on the accession: https://www.ncbi.nlm.nih.gov/protein/acc123/ \n"
        )
        self.assertEqual(result, expected_result)
        mock_popen.assert_not_called()  # Check that subprocess.Popen was not called
        mock_quit_program.assert_not_called()  # Check that quit_program was not called

class TestRunTests(unittest.TestCase):

    def setUp(self):
        # Setup temporary directories and files for testing
        self.test_dir = tempfile.mkdtemp()
        self.primer_file = Path(self.test_dir) / 'primer_123.txt'
        with self.primer_file.open('w') as f:
            f.write("PRIMER_LEFT\nATCG\nPRIMER_RIGHT\nCGTA\nPRIMER_INTERNAL\nGATC\n")

        #test files for success with a count of 2 for target and count of 4 for neighbour
        self.seqkit_file = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m0.txt'
        with self.seqkit_file.open('w') as f:
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT\n")
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT")
        
        self.seqkit_file = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m1.txt'
        with self.seqkit_file.open('w') as f:
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT\n")
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT")
        
        self.seqkit_file = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m3.txt'
        with self.seqkit_file.open('w') as f:
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTGT\n")
            f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTGT")
        
        self.seqkit_file = Path(self.test_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_neighbour_m0.txt'
        with self.seqkit_file.open('w') as f:
          f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTAC\n")
          f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTAC\n")
          f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTAC\n")
          f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTAC")
        # # set up fail directory with two few amplicons for count 1, amplicons of wrong length for target and amplicons of correct length for neighbour
        # self.fail_dir = tempfile.mkdtemp()
        # #m2 must have same amplicon length: FAIL
        # self.seqkit_file = Path(self.fail_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m0.txt'
        # with self.seqkit_file.open('w') as f:
        #     f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT\n")
        # self.seqkit_file = Path(self.fail_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m1.txt'
        # with self.seqkit_file.open('w') as f:
        #     f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTA\n")
        
        # self.seqkit_file = Path(self.fail_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_target_m3.txt'
        # with self.seqkit_file.open('w') as f:
        #     f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCTGT\n")
        
        # self.seqkit_file = Path(self.fail_dir) / 'Test_Primer_123.txt_seqkit_amplicon_against_neighbour_m0.txt'
        # with self.seqkit_file.open('w') as f:
        #     f.write("some\ttab\tseparated\tvalues\tincluding\tamplicon\tATCGGCT\n")
    
    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.test_dir)
        #shutil.rmtree(self.fail_dir)

    def test_count_files(self):
        path = Path(self.test_dir)
        self.assertEqual(count_files(path), 5)

    @patch('Summarize_results_module_improved.Path.glob')
    def test_runtest_no_files(self, mock_glob):
         # Setup mock for glob to return an empty list, simulating no files found
        mock_glob.return_value = []

        #set up expected result
        expected_result = {
            "NA": {
                "Mismatches tested": "NA",
                "Did the test pass?": "NA",
                "Number of assemblies, in silico PCR was performed on": "NA",
                "Number of assemblies with correct size amplicon": "NA",
                "Number of assemblies with wrong size amplicon": "NA"
            },
            "Number of files that passed:": "NA",
            "Number of files that failed:": "NA"
        }

        # Call the function with mocked behavior
        result = run_tests(Path('/some/folder'),"Test_Primer_123.txt", 7,6, "neighbour")

        # Assert the expected output
        self.assertEqual(result, expected_result)

    @patch('Summarize_results_module_improved.Path.glob')
    def test_no_file_number(self, mock_glob):
        mock_file = MagicMock()
        mock_file.name = "test_file_without_pattern.txt"
        mock_glob.return_value = [mock_file]
        with self.assertRaises(ValueError) as cm:
            run_tests(Path('/fake/folder'), 'test', 0, 0, 'target')
        
        # Check the exception message
        self.assertEqual(str(cm.exception), "No match found in the document name.")

    @patch('Summarize_results_module_improved.Path.open')
    @patch('sys.exit')
    def test_os_error(self, mock_exit,mock_open):
        # Simulate OSError when trying to open a file
        mock_open.side_effect = OSError("Unable to open file")
        
        with self.assertRaises(OSError):
            quit_program("could not open or read file: Unable to open file", OSError("Unable to open file"))

    # # #ASK HENRY!
    # @patch('Summarize_results_module_improved.Path.glob')
    # @patch('Summarize_results_module_improved.Path.open')
    # @patch('Summarize_results_module_improved.quit_program')
    # def test_os_error_script(self, mock_quit_program, mock_open, mock_glob):

    #     # Setting up mock for files
    #     mock_file = MagicMock()
    #     mock_file.name = "test_file_m1.txt"
    #     mock_glob.return_value = [mock_file]

    #     # Simulate OSError when trying to open a file
    #     mock_open.side_effect = OSError("Unable to open file")
        
    #     with self.assertRaises(OSError):
    #         run_tests(Path('/fake/folder'), "test_file_m1", 0, 0, 'target')
        
    #     # Ensure quit_program is called with the correct message
    #     mock_quit_program.assert_called_once_with("could not open or read file: Unable to open file", OSError("Unable to open file"))

   
    def test_sensitvity_successfully(self):

        expected_dict = {"Test_Primer_123.txt_seqkit_amplicon_against_target_m0.txt":{
                "Mismatches tested": "m0",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 2,
                "Number of assemblies with wrong size amplicon": 0},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m1.txt":{
                "Mismatches tested": "m1",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 2,
                "Number of assemblies with wrong size amplicon": 0},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m3.txt":{
                "Mismatches tested": "m3",
                "Did the test pass?": "failed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 0,
                "Number of assemblies with wrong size amplicon": 2}, "Number of files that passed:": 2
                ,"Number of files that failed:": 1, "Number of files (in silico PCR result files for different number of primer mismatches) tested:":3
            }
        #mock run
        result =run_tests(Path(self.test_dir), "Test_Primer_123.txt", 7,2,"target")

        # Assert the expected output
        self.assertEqual(result, expected_dict)
    def test_specificity_successfully(self):

        expected_dict = {
                "Test_Primer_123.txt_seqkit_amplicon_against_neighbour_m0.txt":{
                "Mismatches tested": "m0",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 4,
                "Number of assemblies with correct size amplicon": 0}, "Number of files that passed:": 1
                ,"Number of files that failed:": 0, "Number of files (in silico PCR result files for different number of primer mismatches) tested:":1
            }
        #mock run
        result =run_tests(Path(self.test_dir), "Test_Primer_123.txt", 7,4,"neighbour")

        # Assert the expected output
        self.assertEqual(result, expected_dict)
    def test_specificity_failed(self):
        expected_dict = {
                "Test_Primer_123.txt_seqkit_amplicon_against_neighbour_m0.txt":{
                "Mismatches tested": "m0",
                "Did the test pass?": "failed",
                "Number of assemblies, in silico PCR was performed on": 10,
                "Number of assemblies with correct size amplicon": 4}, "Number of files that passed:": 0
                ,"Number of files that failed:": 1, "Number of files (in silico PCR result files for different number of primer mismatches) tested:":1
            }
        #mock run
        result =run_tests(Path(self.test_dir), "Test_Primer_123.txt", 9,10,"neighbour")

        # Assert the expected output
        self.assertEqual(result, expected_dict)
    
class TestInterpretAndReformatSensiSpeciTests(unittest.TestCase):


    def test_not_a_dictionary_str(self):
        """See if TypeError is raised if string and not dictionary is provided"""
        test_str="This is a string"
        
        with self.assertRaises(TypeError):
            interpret_and_reformat_sensi_speci_tests(test_str,"target")
    
    def test_not_a_dictionary_tuple(self):
        """Test is TypeError is thrown for Tuple"""
        test_tuple=(1, 3, 5, 7, 9)

        with self.assertRaises(TypeError):
            interpret_and_reformat_sensi_speci_tests(test_tuple,"target")
    
    def test_interpret_and_reformat_sensi_speci_fail(self):
        """Test it successfully interprets a fail as a fail"""
        input_dict = {"Test_Primer_123.txt_seqkit_amplicon_against_target_m0.txt":{
                "Mismatches tested": "m0",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 2,
                "Number of assemblies with wrong size amplicon": 0},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m1.txt":{
                "Mismatches tested": "m1",
                "Did the test pass?": "failed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 1,
                "Number of assemblies with wrong size amplicon": 1},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m3.txt":{
                "Mismatches tested": "m3",
                "Did the test pass?": "failed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 0,
                "Number of assemblies with wrong size amplicon": 2}, "Number of files that passed:": 2
                ,"Number of files that failed:": 1, "Number of files (in silico PCR result files for different number of primer mismatches) tested:":3
            }
        expected_output=("FAILED",2,"m1, m3")

        result=interpret_and_reformat_sensi_speci_tests(input_dict, "target")

        self.assertEqual(result, expected_output)

    def test_interpret_and_reformat_sensi_speci_passed(self):
        """Test it successfully interprets a fail as a fail"""
        input_dict = {"Test_Primer_123.txt_seqkit_amplicon_against_target_m0.txt":{
                "Mismatches tested": "m0",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 2,
                "Number of assemblies with wrong size amplicon": 0},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m1.txt":{
                "Mismatches tested": "m1",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 2,
                "Number of assemblies with wrong size amplicon": 0},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m3.txt":{
                "Mismatches tested": "m3",
                "Did the test pass?": "failed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 0,
                "Number of assemblies with wrong size amplicon": 2}, "Number of files that passed:": 2
                ,"Number of files that failed:": 1, "Number of files (in silico PCR result files for different number of primer mismatches) tested:":3
            }
        expected_output=("PASSED",2,"m3")

        result=interpret_and_reformat_sensi_speci_tests(input_dict, "target")

        self.assertEqual(result, expected_output)
    
    def test_correct_output_type(self):
        """Test whether the output is of the correct type"""
        input_dict = {"Test_Primer_123.txt_seqkit_amplicon_against_target_m0.txt":{
                "Mismatches tested": "m0",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 2,
                "Number of assemblies with wrong size amplicon": 0},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m1.txt":{
                "Mismatches tested": "m1",
                "Did the test pass?": "passed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 2,
                "Number of assemblies with wrong size amplicon": 0},
                "Test_Primer_123.txt_seqkit_amplicon_against_target_m3.txt":{
                "Mismatches tested": "m3",
                "Did the test pass?": "failed",
                "Number of assemblies, in silico PCR was performed on": 2,
                "Number of assemblies with correct size amplicon": 0,
                "Number of assemblies with wrong size amplicon": 2}, "Number of files that passed:": 2
                ,"Number of files that failed:": 1, "Number of files (in silico PCR result files for different number of primer mismatches) tested:":3
            }
    
        test, numb, m =interpret_and_reformat_sensi_speci_tests(input_dict, "target")

        #assert the output is of the correct type
        self.assertIsInstance(test, str)
        self.assertIsInstance(numb,int)
        self.assertIsInstance(m, str)



if __name__ == '__main__':
    unittest.main()
