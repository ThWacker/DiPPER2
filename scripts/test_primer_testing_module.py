#! /usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
from unittest.mock import patch, MagicMock, Mock
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

class test_file_ops(unittest.TestCase):

    def setUp(self):
        # Setup temporary directories and files for testing
        self.folder = tempfile.mkdtemp()
        self.source=tempfile.mkdtemp()
        self.file1=Path(self.folder)/"file1_because_I_am_creative"
        with self.file1.open('w') as f:
            f.write("Not sure I hate testing \n")
        self.file2 = Path(self.folder) / 'file2_oh_the_creativity'
        with self.file2.open('w') as f:
            f.write("or I love it \n")
        # make an empty folder
        self.empty=tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.folder)
        shutil.rmtree(self.source)
        shutil.rmtree(self.empty)
    
    @patch('Primer_Testing_module_optimized.Logger')
    def test_concat_files_success(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        #run concat files
        concat_files(Path(self.folder), "target", Path(self.source), mock_logger_instance)

        #assert that what we want is what we get
        outfile=Path(self.source)/"target_concatenated.fasta"
        content_outfile=Path(outfile).read_text()
        self.assertEqual(content_outfile, "Not sure I hate testing \nor I love it \n")

    @patch('Primer_Testing_module_optimized.Path.glob')
    @patch('Primer_Testing_module_optimized.Logger')
    def test_no_files_found_exception_bycount(self, mock_logger_class, mock_glob):
        #magic mock our logger
        mock_logger_instance = MagicMock()#Who even came up with the name MagicMock? 
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Mock glob to return no files
        mock_glob.return_value = []

        with self.assertRaises(FileNotFoundError):
            concat_files(Path(self.empty), "target", Path(self.source), mock_logger_instance)

        mock_logger_instance.exception.assert_called_once_with(
            f"Could not open or read files in {Path(self.empty)}. Concatenation failed.",
            exc_info=1,
        )
    @patch('Primer_Testing_module_optimized.Path.glob')
    @patch('Primer_Testing_module_optimized.Logger')
    def test_concat_files_exception(self, mock_logger_class, mock_glob):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Mock glob to return a list of one file
        mock_file = MagicMock()
        mock_file.name = "file1.fasta"
        mock_glob.return_value = [mock_file]

        # Mock file.open() to raise FileNotFoundError
        mock_file.open.side_effect = FileNotFoundError("Mocked FileNotFoundError")

        # Run concat_files and check for exception
        with self.assertRaises(FileNotFoundError):
            concat_files(Path(self.empty), "target", Path(self.source), mock_logger_instance)

        # Ensure the logger was called with the correct exception
        mock_logger_instance.exception.assert_called_once_with(
            f"Could not open or read files in {self.empty}. Concatenation failed.", exc_info=1,
        )

    @patch('Primer_Testing_module_optimized.Logger')
    def test_delete_concats(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Create two temporary files
        self.temp1 = tempfile.mkdtemp()
        self.temp2 = tempfile.mkdtemp()
        self.target = Path(self.temp1) / "target.tmp"
        self.neighbour = Path(self.temp2) / "neighbour.tmp"
        self.target.touch()  # Create the target file
        self.neighbour.touch()  # Create the neighbour file

        # Run delete_concats
        delete_concats(self.target, self.neighbour, mock_logger_instance)

        # Assert the logger info message was logged
        mock_logger_instance.info.assert_called_once()
        self.assertIn("Concatenated files deleted.", mock_logger_instance.info.call_args[0][0])

        # Assert the files no longer exist
        self.assertFalse(self.target.exists())
        self.assertFalse(self.neighbour.exists())

    @patch('Primer_Testing_module_optimized.Logger')
    def test_delete_concats_files_missing(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Create paths that do not exist
        self.target = Path("non_existent_target.tmp")
        self.neighbour = Path("non_existent_neighbour.tmp")

        # Run delete_concats
        delete_concats(self.target, self.neighbour, mock_logger_instance)

        # Assert logger.error was called
        mock_logger_instance.error.assert_called_once()
        self.assertIn("Error deleting concatenated files:", mock_logger_instance.error.call_args[0][0])

class test_extract_primer_sequences(unittest.TestCase):
    def setUp(self):
        # Setup temporary directories and files for testing
        self.source_dir = tempfile.mkdtemp()
        self.file1 = Path(self.source_dir) / "file1.fasta"
        data= '''
                PRIMER_LEFT_0_SEQUENCE
                TAGTGTCAGACCCTAGGGCC
                PRIMER_RIGHT_0_SEQUENCE
                CTAGCAATCTGGGCAGCTGT
                PRIMER_INTERNAL_0_SEQUENCE
                ACAAGCCTCCCATGCCAGGGCG
            '''
        with self.file1.open("w") as f:
            f.write(data)
        self.file2 = Path(self.source_dir) / "file2.fasta"
        data= '''
                PRIMER_LEFT_0_SEQUENCE
                TAGTGTCAGACCCTAGGGCC
                PRIMER_RIGHT_0_SEQUENCE
                CTAGCAATCTGGGCAGCTGT
                PRIMER_INTERNAL_0_SEQUENCE
            '''
        with self.file2.open("w") as f:
            f.write(data)
    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.source_dir)

    @patch("Primer_Testing_module_optimized.Logger")
    def test_extract_primer_sequence_success(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance 
        
        #build dictionary
        dictionary=( 'TAGTGTCAGACCCTAGGGCC','CTAGCAATCTGGGCAGCTGT','ACAAGCCTCCCATGCCAGGGCG')

        #run and get results
        results=extract_primer_sequences(Path(self.file1),mock_logger_instance)
        #assert it is correct
        self.assertEqual(dictionary, results)

    @patch("sys.exit")
    @patch("Primer_Testing_module_optimized.Logger")
    def test_extract_primer_sequence_fail(self, mock_logger_class, mock_exit):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance 

        #run and get results
        results=extract_primer_sequences(Path(self.file2),mock_logger_instance)

        #assert sys.exit was called
        mock_exit.assert_called_once()

        #check if logger was correctly called
        mock_logger_instance.error.assert_called_once()
        self.assertIn("The primer headers are not correctly formatted and cannot be processed. Please change the headers accordingly.", mock_logger_instance.error.call_args[0][0])

class test_move_files_with_pattern(unittest.TestCase):

    def setUp(self):
        # Setup temporary directories and files for testing
        self.source_dir = tempfile.mkdtemp()
        self.destination_dir = tempfile.mkdtemp()
        self.file1 = Path(self.source_dir) / "file1_pattern.fasta"
        data= ""
        with self.file1.open("w") as f:
            f.write(data)

    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.source_dir)
        shutil.rmtree(self.destination_dir)

    @patch("shutil.move")
    #@patch("pathlib.Path.mkdir")
    @patch("pathlib.Path.iterdir")
    def test_no_matching_files(self, mock_iterdir,  mock_shutil_move):
        # Mock source and destination directories
        source_dir = MagicMock(spec=Path)
        destination_dir = MagicMock(spec=Path)

        # Mock no matching files in the source directory
        mock_file1 = MagicMock(spec=Path)
        mock_file1.name = "file1.txt"
        mock_file1.is_file.return_value = True

        mock_iterdir.return_value = [mock_file1]

        # Call the function
        move_files_with_pattern(source_dir, "pattern", destination_dir)

        # Assert mkdir was called
        destination_dir.mkdir.assert_called_once_with(exist_ok=True) 

        # Assert shutil.move was not called since no files matched
        mock_shutil_move.assert_not_called()

    @patch("pathlib.Path.mkdir")
    @patch("Primer_Testing_module_optimized.shutil.move")
    def test_matching_files(self,  mock_shutil_move, mock_mkdir):
      
        # Call the function
        move_files_with_pattern(Path(self.source_dir), "pattern", Path(self.destination_dir))

        # Assert mkdir was called
        mock_mkdir.assert_called_once_with(exist_ok=True) 

        # Assert shutil.move was not called since no files matched
        mock_shutil_move.assert_called_once()
    
if __name__ == "__main__":
    unittest.main()