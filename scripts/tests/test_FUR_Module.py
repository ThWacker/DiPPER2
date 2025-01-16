#!/usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
import subprocess
from unittest.mock import patch, MagicMock
from FUR_module_optimized import (make_fur_db, run_fur, check_folders, clean_up)
from logging_handler import Logger

class TestUtilFunctions (unittest.TestCase):

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
    
    @patch('FUR_module_optimized.Logger')
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
    
    @patch('FUR_module_optimized.Logger')
    @patch('FUR_module_optimized.Path.iterdir', return_value=[])
    @patch('FUR_module_optimized.Path.is_dir', return_value=True)
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


class TestFUR(unittest.TestCase):
        
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.neighbour_dir=tempfile.mkdtemp()
        self.assembly_mock=Path(self.neighbour_dir)/"assembly.fa"
        with self.assembly_mock.open('w') as f:
            f.write(">MockAssembly GCA123\nATGCGCTAA")
        self.target_dir=tempfile.mkdtemp()
        self.assembly_target_mock=Path(self.neighbour_dir)/"assembly_target.fa"
        with self.assembly_target_mock.open('w') as f:
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
        shutil.rmtree(self.target_dir)
    
    @patch('FUR_module_optimized.Logger')
    @patch("FUR_module_optimized.subprocess.run")
    def test_makeFurDB_error(self, mock_run, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        #mock side effect, assert the correct exception is raised
        mock_run.side_effect=RuntimeError(1, 'makeFurDb')
        with self.assertRaises(RuntimeError):
            make_fur_db(Path(self.target_dir), Path(self.neighbour_dir), Path(self.test_dir), mock_logger_instance)
            # Check if logger.exception was called with the expected message
            mock_logger_instance.exception.assert_called_once()
            self.assertTrue("Error when executing makeFurDb:" in mock_logger_instance.exception.call_args[0][0])

    @patch('FUR_module_optimized.Logger')
    @patch("FUR_module_optimized.subprocess.run")
    @patch("FUR_module_optimized.Path.write_text")
    @patch("FUR_module_optimized.Path.stat")
    @patch("FUR_module_optimized.Path.exists")
    @patch('builtins.print')
    @patch('sys.exit')
    def test_furout_empty(self, mock_exit, mock_print, mock_exists, mock_stat, mock_write_text, mock_run, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
        # Arrange
        mock_exists.return_value = True  # Pretend file exists
        mock_stat.return_value.st_size = 0  # Pretend file is empty
        mock_run.return_value.returncode = 0  # Pretend subprocess.run is successful

        # Test directory and file name
        source = Path('/some/path')
        outfile_prefix = 'test_prefix'
        option = 'some_option'
      
        # run the actual function
        run_fur(option, outfile_prefix, source, mock_logger_instance)

        # Assert
        # Check if logger.exception was called with the expected message
        mock_logger_instance.error.assert_called_once()
        self.assertTrue(f"Fur could not find unique regions or did not run successfully. {source / f'{outfile_prefix}_FUR.db.out.txt'} is empty or does not exist." in mock_logger_instance.error.call_args[0][0])

    @patch('FUR_module_optimized.Logger')
    @patch("FUR_module_optimized.subprocess.run")
    @patch("FUR_module_optimized.Path.write_text")
    @patch("FUR_module_optimized.Path.stat")
    @patch("FUR_module_optimized.Path.exists")
    @patch('builtins.print')
    @patch('sys.exit')
    def test_furout_not_existing(self, mock_exit, mock_print, mock_exists, mock_stat, mock_write_text, mock_run, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Arrange
        mock_exists.return_value = False  # Pretend file doesn't exist
        mock_stat.return_value.st_size = 1  # Pretend file is not empty
        mock_run.return_value.returncode = 0  # Pretend subprocess.run is successful

        # Test directory and file name
        source = Path('/some/path')
        outfile_prefix = 'test_prefix'
        option = 'some_option'
       
        # run the actual function
        run_fur(option, outfile_prefix, source, mock_logger_instance)

        # Assert
        # Check if logger.exception was called with the expected message
        mock_logger_instance.error.assert_called_once()
        self.assertTrue(f"Fur could not find unique regions or did not run successfully. {source / f'{outfile_prefix}_FUR.db.out.txt'} is empty or does not exist." in mock_logger_instance.error.call_args[0][0])

    @patch('FUR_module_optimized.Logger')
    @patch("FUR_module_optimized.subprocess.run")
    @patch('sys.exit')
    def test_FUR_error(self,mock_exit, mock_run, mock_logger_class):

        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance
        mock_run.side_effect=subprocess.CalledProcessError(1, 'fur')
        with self.assertRaises(RuntimeError):
            make_fur_db(Path(self.target_dir), Path(self.neighbour_dir), Path(self.test_dir), mock_logger_instance)

    @patch('FUR_module_optimized.Logger')
    @patch('FUR_module_optimized.subprocess.run')
    def test_successful_makeFurDb(self, mock_run, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        target_folder = Path('/fake/target')
        neighbour_folder = Path('/fake/neighbour')
        source = Path('/fake/source')

        make_fur_db(target_folder, neighbour_folder, source, mock_logger_instance)
        mock_run.assert_called_once_with(['makeFurDb', '-t', str(target_folder), '-n', str(neighbour_folder), '-d', str(source / 'FUR.db')], check=True)

    
    @patch('FUR_module_optimized.Logger')
    @patch('FUR_module_optimized.subprocess.run')
    def test_successful_makeFurDb_reference(self, mock_run, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        target_folder = Path('/fake/target')
        neighbour_folder = Path('/fake/neighbour')
        source = Path('/fake/source')
        reference= Path ('/fake/reference.fasta')

        make_fur_db(target_folder, neighbour_folder, source, mock_logger_instance, reference)
        mock_run.assert_called_once_with(['makeFurDb', '-t', str(target_folder), '-n', str(neighbour_folder),'-d', str(source / 'FUR.db'),'-r', str(reference)], check=True)
    
    @patch('FUR_module_optimized.Logger')
    @patch('FUR_module_optimized.Path.write_text')
    @patch('FUR_module_optimized.Path.stat')
    @patch('FUR_module_optimized.Path.exists')
    @patch('FUR_module_optimized.subprocess.run')
    def test_run_fur(self, mock_subprocess_run, mock_exists, mock_stat, mock_write_text, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # make sure it thinks the file exists and is not empty
        mock_exists.return_value = True
        mock_stat.return_value.st_size = 100
        mock_subprocess_run.return_value.stdout = "FUR output"
        mock_subprocess_run.return_value.stderr = ""
        
        #parameters
        option = "u"
        outfile_prefix = "test_prefix"
        source = Path('/fake/source')

        # run and assert
        run_fur(option, outfile_prefix, source, mock_logger_instance)
        mock_subprocess_run.assert_called_once_with(['fur', '-d', str(source / 'FUR.db'), '-u'], check=True, capture_output=True, text=True)
        mock_write_text.assert_called_once_with("FUR output", encoding="utf-8")

    @patch('FUR_module_optimized.Logger')
    @patch('FUR_module_optimized.Path.write_text')
    @patch('FUR_module_optimized.Path.stat')
    @patch('FUR_module_optimized.Path.exists')
    @patch('FUR_module_optimized.subprocess.run')
    def test_run_fur_error(self, mock_subprocess_run, mock_exists, mock_stat, mock_write_text, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # mock run outcome as processerror
        mock_subprocess_run.side_effect=subprocess.CalledProcessError(1, 'fur')

        # make sure it thinks the file exists and is not empty
        mock_exists.return_value = True
        mock_stat.return_value.st_size = 100
    
        #parameters
        option = "u"
        outfile_prefix = "test_prefix"
        source = Path('/fake/source')

        #assert
        with self.assertRaises(RuntimeError):
            run_fur(option, outfile_prefix, source, mock_logger_instance)

        # Check if logger.exception was called with the info message
        mock_logger_instance.exception.assert_called_once()
        self.assertTrue(
            "FUR did not run to completion and a Runtime Error occured: Command 'fur' returned non-zero exit status 1" in mock_logger_instance.exception.call_args[0][0]
        )
        
class TestCleanUp(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.fur_output = Path(self.test_dir) / 'FUR.db'
        with self.fur_output.open('w') as f:
            f.write("a database")

    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.test_dir)
    
    @patch('FUR_module_optimized.shutil.rmtree')
    @patch('FUR_module_optimized.Logger')
    def test_clean_up(self, mock_logger_class, mock_shutil_rmtree):

         # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        clean_up(Path(self.fur_output), mock_logger_instance)
        mock_shutil_rmtree.assert_called_once()
        mock_logger_instance.info.assert_called_once()
        self.assertTrue(
            f"The directory {self.fur_output} has been deleted." in mock_logger_instance.info.call_args[0][0]
        )
    
    @patch('FUR_module_optimized.shutil.rmtree')
    @patch('FUR_module_optimized.Logger')
    def test_clean_up_nodb(self, mock_logger_class, mock_shutil_rmtree):

         # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        clean_up(Path("notaDB"), mock_logger_instance)
        mock_shutil_rmtree.assert_not_called()
        mock_logger_instance.info.assert_not_called()
       

if __name__ == '__main__':
    unittest.main()
        

