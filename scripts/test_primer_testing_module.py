#! /usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from Primer_Testing_module_optimized import (
    check_folders,
    check_program_installed,
    concat_files,
    delete_concats,
    extract_primer_sequences,
    move_files_with_pattern,
    get_amplicon,
    get_longest_target,
    run_seqkit_amplicon_with_optional_timeout,
    run_seqkit_locate,
)

class test_util_functions (unittest.TestCase):

    def setUp(self):
    # Setup temporary directories and files for testing
        self.test_dir = tempfile.mkdtemp()
        self.neighbour_dir=tempfile.mkdtemp()
        self.assembly_mock=Path(self.neighbour_dir)/"assembly.fa"
        with self.assembly_mock.open('w') as f:
            f.write(">MockAssembly GCA123\nATGCGCTAA")
        self.primer_file = Path(self.test_dir) / 'primer_123.txt'
        with self.primer_file.open('w') as f:
            f.write("PRIMER_LEFT\nATCG\nPRIMER_RIGHT\nCGTA\nPRIMER_INTERNAL\nGATC\n")
        self.fur_output = Path(self.test_dir) / 'FUR.db'
        with self.fur_output.open('w') as f:
            f.write("a database")

    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.test_dir)
        shutil.rmtree(self.neighbour_dir)
    
    @patch('Primer_Testing_module_optimized.Logger')
    def test_check_folders(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        #first one should not raise
        check_folders(Path(self.test_dir), Path(self.neighbour_dir), logger=mock_logger_instance)

        #this should sys.exit
        with self.assertRaises(SystemExit):
            check_folders(Path("/Nonesensefolder"), Path("/AsNoneExistentAsMyIQ"), logger=mock_logger_instance)

        #check if error is correctly logged
        mock_logger_instance.error.assert_called_once()
        self.assertTrue("The folder /Nonesensefolder does not exist" in mock_logger_instance.error.call_args[0][0])
    
    @patch('Primer_Testing_module_optimized.Logger')
    @patch('Primer_Testing_module_optimized.Path.iterdir', return_value=[])
    @patch('Primer_Testing_module_optimized.Path.is_dir', return_value=True)
    @patch('sys.exit')
    def test_check_folders_empty(self, mock_exit, mock_is_dir, mock_iterdir, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
        # Test with empty directory
        target_folder = Path('/fake/target')
        neighbour_folder = Path('/fake/neighbour')
        check_folders(target_folder, neighbour_folder, logger=mock_logger_instance)
        #check if error is correctly logged
        #Verify logger error calls for non-existant files
        expected_error_calls = [
            unittest.mock.call('The folder /fake/target is empty'),
            unittest.mock.call('The folder /fake/neighbour is empty')
        ]
        mock_logger_instance.error.assert_has_calls(expected_error_calls, any_order=True)



    def test_check_program_installed_found(self):
        
        # test with program that should always be installed
        result=check_program_installed("cd")
        self.assertTrue(result)

    def test_check_program_installed_None(self):
        # test with program that should always be installed
        result=check_program_installed("Superkalifragilistikexpialegetisch")
        self.assertFalse(result)
    
    def test_check_program_installed_typeerror(self):
        #test that it raises typeerror when confronted with wrong type
        boolean_var=True
        with self.assertRaises(TypeError):
            check_program_installed(boolean_var)

#class test_file_operations(unittest.TestCase):

        
if __name__ == "__main__":
    unittest.main()