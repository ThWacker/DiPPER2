#!/usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
import subprocess
#from logging_handler import Logger
from unittest.mock import patch, MagicMock

# Assuming all your functions are imported from your script
from target_move_module_optimized import (parse_list, process_target_files,move_remaining_files_to_neighbour)

class TestParseList(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.target1_file = Path(self.test_dir) / 'target1.fasta'
        self.target2_file=Path(self.test_dir) / 'target2.fasta'
        self.target3_file=Path(self.test_dir) / 'target3.fa'
        self.target4_file=Path(self.test_dir) / 'target4.fna'
        self.target5_file=Path(self.test_dir)/ 'target5.txt'
        self.neighbour1_file=Path(self.test_dir) / 'neighbour1.fasta'
        self.neighbour2_file=Path(self.test_dir) / 'neighbour2.fna'
        self.listfile=Path(self.test_dir) / "list.txt"
        with  self.listfile.open('w')  as f:
            f.write("target1\n")
            f.write("target2\n")
            f.write("target3\n")
            f.write("target4\n")
            f.write("target5")
    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.test_dir)
    
    @patch('target_move_module_optimized.Logger')
    def test_parse_list_logger(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Simulate FileNotFoundError by providing a non-existent file path
        non_existent_file = Path("non_existent_file.txt")
        
        # Assert that FileNotFoundError is raised
        with self.assertRaises(FileNotFoundError):
            logger = mock_logger_instance
            parse_list(non_existent_file, logger)
        
        # Check if logger.exception was called with the expected message
        mock_logger_instance.exception.assert_called_once()
        self.assertTrue("Error: Target file 'non_existent_file.txt' not found!" in mock_logger_instance.exception.call_args[0][0])

    @patch('target_move_module_optimized.Logger')
    def test_parse_list_return_is_list(self, mock_logger_class):
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
        result=parse_list(Path(self.listfile), mock_logger_instance)
        self.assertIsInstance(result,list)



class testProcessTargetFiles(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.target1_file = Path(self.test_dir) / 'target1.fasta'
        self.target2_file=Path(self.test_dir) / 'target2.fasta'
        self.target3_file=Path(self.test_dir) / 'target3.fa'
        self.target4_file=Path(self.test_dir) / 'target4.fna'
        self.target5_file=Path(self.test_dir)/ 'target5.txt'
        self.neighbour1_file=Path(self.test_dir) / 'neighbour1.fasta'
        self.neighbour2_file=Path(self.test_dir) / 'neighbour2.fna'
        self.listfile=Path(self.test_dir) / "list.txt"
        with  self.listfile.open('w')  as f:
            f.write("target1\n")
            f.write("target2\n")
            f.write("target3\n")
            f.write("target4\n")
            f.write("target5")
    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.test_dir)
    
    #@patch('target_move_module_optimized.parse_list')
    @patch('target_move_module_optimized.Logger')
    @patch('target_move_module_optimized.subprocess.run')
    @patch('target_move_module_optimized.shutil.copy')
    @patch('target_move_module_optimized.Path.cwd')
    def test_fasta_skip(self,mock_cwd, mock_copy, mock_run, mock_logger_class):
        # mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
    
        # Mock the current working directory
        mock_cwd.return_value = Path(self.test_dir)

        # Mock the subprocess.run to return a non-fasta file
        mock_run.return_value.stdout = self.target5_file.name
        mock_run.return_value.returncode = 0
        
        #mock folders & list
        source_folder = Path(self.test_dir)
        fur_target_folder = Path('/mock/dist/FUR.target')
        list=["target5"]

        #run 
        process_target_files(list,source_folder,fur_target_folder, logger=mock_logger_instance)
        # Assert the file was not copied
        mock_copy.assert_not_called()
    
    @patch('target_move_module_optimized.Logger')
    @patch('target_move_module_optimized.subprocess.run')
    @patch('target_move_module_optimized.shutil.copy')
    @patch('target_move_module_optimized.Path.cwd')
    def test_process_target_files_success(self,mock_cwd, mock_copy, mock_run, mock_logger_class):
        # mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
        # Mock the subprocess.run to return different files based on the pattern
        def side_effect_find_command(args, **kwargs):
            accession_pattern = args[-1].strip("*")
            if accession_pattern == "target1":
                return MagicMock(stdout=str(self.target1_file), returncode=0)
            elif accession_pattern == "target2":
                return MagicMock(stdout=str(self.target2_file), returncode=0)
            elif accession_pattern == "target3":
                return MagicMock(stdout=str(self.target3_file), returncode=0)
            elif accession_pattern == "target4":
                return MagicMock(stdout=str(self.target4_file), returncode=0)
            elif accession_pattern == "target5":
                return MagicMock(stdout=str(self.target5_file), returncode=0)
            else:
                return MagicMock(stdout="", returncode=0)

        mock_run.side_effect = side_effect_find_command

        source_folder = Path(self.test_dir)
        fur_target_folder = Path('/mock/dist/FUR.target')
        target_list = ["target1", "target2", "target3", "target4", "target5"]

        process_target_files(target_list, source_folder, fur_target_folder, logger=mock_logger_instance)

        # Assert that there have been 4 files copied (target1, target2, target3, target4)
        self.assertEqual(mock_copy.call_count, 4)

        # check if correct files were copied
        mock_copy.assert_any_call(self.target1_file, fur_target_folder)
        mock_copy.assert_any_call(self.target2_file, fur_target_folder)
        mock_copy.assert_any_call(self.target3_file, fur_target_folder)
        mock_copy.assert_any_call(self.target4_file, fur_target_folder)

    @patch('target_move_module_optimized.Logger')
    @patch('target_move_module_optimized.subprocess.run')
    @patch('target_move_module_optimized.shutil.copy')
    @patch('target_move_module_optimized.Path.cwd')
    def test_file_not_found(self, mock_cwd, mock_copy, mock_run, mock_logger_class):
        # mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
        # Mock the current working directory
        mock_cwd.return_value = Path(self.test_dir)
        
        # Mock the subprocess.run to return no file
        mock_run.return_value.stdout = ''
        mock_run.return_value.returncode = 0

        # Mock list_targets and the dist folder
        list_targets = ["target1", "target2", "target3", "target4", "target5"]
        source_folder = Path(self.test_dir)
        fur_target_folder = Path('/mock/dist/FUR.target')

        # Call the process_targets function directly
        process_target_files(list_targets, source_folder, fur_target_folder, logger=mock_logger_instance)

        # Assert the file was not copied
        mock_copy.assert_not_called()
    
    @patch('target_move_module_optimized.Logger')
    @patch('target_move_module_optimized.subprocess.run')
    @patch('target_move_module_optimized.Path.cwd')
    def test_Called_Process_Error(self,mock_cwd, mock_run, mock_logger_class):
        # mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
        # Mock the current working directory
        mock_cwd.return_value = Path(self.test_dir)
        #mock sideeffect
        mock_run.side_effect = subprocess.CalledProcessError(1, 'find')

        #run all
        source_folder = Path(self.test_dir)
        fur_target_folder = Path('/mock/dist/FUR.target')
        target_list = ["target1", "target2", "target3", "target4", "target5"]

        # assert the error is called
        with self.assertRaises(RuntimeError):
            process_target_files(target_list, source_folder, fur_target_folder, logger=mock_logger_instance)
    
    
class TestMoveRemainingFilesToNeighbour(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for testing
        self.source_folder = tempfile.mkdtemp()
        self.fur_target_folder = tempfile.mkdtemp()
        self.fur_neighbour_folder = tempfile.mkdtemp()

        # Create test files in the source folder
        self.files = {
            'file1.fasta': 'ATGTCGAAA',
            'file2.fasta': '>This is a fasta file\n ATGCTCCAGTC',
            'file3.txt': 'Superkalifragilistikexpialegetisch',
            'file4.fa':"This will be copied"
        }

        for filename in self.files:
            file_path = Path(self.source_folder) / filename
            with file_path.open('w') as f:
                f.write(self.files[filename])

        # Create target files to simulate existing files
        for filename in ['file1.fasta', 'file2.fasta']:
            file_path = Path(self.fur_target_folder) / filename
            with file_path.open('w') as f:
                f.write('')

    def tearDown(self):
        # Remove temporary directories and files after testing
        shutil.rmtree(self.source_folder)
        shutil.rmtree(self.fur_target_folder)
        shutil.rmtree(self.fur_neighbour_folder)

    @patch('target_move_module_optimized.Logger')
    def test_file_exists_in_target_folder(self, mock_logger_class):
        # Create mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Call the function
        move_remaining_files_to_neighbour(
            Path(self.source_folder), 
            Path(self.fur_target_folder), 
            Path(self.fur_neighbour_folder), 
            mock_logger_instance
        )

        # Verify logger debug calls for skipped files
        expected_debug_calls = [
            unittest.mock.call(f'File file1.fasta already exists in {Path(self.fur_target_folder)}, skipping.'),
            unittest.mock.call(f'File file2.fasta already exists in {Path(self.fur_target_folder)}, skipping.')
        ]
        mock_logger_instance.debug.assert_has_calls(expected_debug_calls, any_order=True)

        # Verify warning calls for non-fasta files
        mock_logger_instance.warning.assert_called_once_with("file3.txt is not a fasta file. Skipping.")

    @patch('target_move_module_optimized.Logger')
    @patch('target_move_module_optimized.shutil.copy')
    def test_successful_move(self, mock_copy, mock_logger_class):
        # Mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Call the function
        move_remaining_files_to_neighbour(
            Path(self.source_folder), 
            Path(self.fur_target_folder), 
            Path(self.fur_neighbour_folder), 
            mock_logger_instance
        )

        # Verify that shutil.copy was called with the correct arguments
        mock_copy.assert_called_once_with(
            Path(self.source_folder) / 'file4.fa', 
            Path(self.fur_neighbour_folder)
        )

        # Verify debug call for copied file
        mock_logger_instance.debug.assert_any_call(f'Copied to {Path(self.fur_neighbour_folder)}: file4.fa')


if __name__ == '__main__':
    unittest.main()
        

