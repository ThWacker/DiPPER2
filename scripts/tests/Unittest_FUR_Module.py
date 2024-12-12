#!/usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
import subprocess
from unittest.mock import patch

# Assuming all your functions are imported from your script
from FUR_module_optimized import (quit, make_fur_db, run_fur, check_folders, usage)

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

    def test_quit_program(self):
        with self.assertRaises(SystemExit):
            quit("More tests")
    
    def test_usage(self):
        with patch('builtins.print') as mocked_print:
            usage()
            self.assertTrue(mocked_print.called)    
    
    def test_check_folders(self):
        #first one should not raise
        check_folders(Path(self.test_dir), Path(self.neighbour_dir))
        #this should sys.exit
        with self.assertRaises(SystemExit):
            check_folders(Path("/Nonesensefolder"), Path("/AsNoneExistentAsMyIQ"))
    
    @patch('FUR_module_optimized.Path.iterdir', return_value=[])
    @patch('FUR_module_optimized.Path.is_dir', return_value=True)
    @patch('FUR_module_optimized.quit')
    def test_check_folders_empty(self, mock_quit, mock_is_dir, mock_iterdir):
        # Test with empty directory
        target_folder = Path('/fake/target')
        neighbour_folder = Path('/fake/neighbour')
        check_folders(target_folder, neighbour_folder)
        mock_quit.assert_called_with(f"The folder {neighbour_folder} is empty")


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
    
    @patch("FUR_module_optimized.subprocess.run")
    @patch('sys.exit')
    def test_makeFurDB_error(self,mock_exit, mock_run):
        mock_run.side_effect=subprocess.CalledProcessError(1, 'makeFurDb')
        with self.assertRaises(subprocess.CalledProcessError):
            make_fur_db(Path(self.target_dir), Path(self.neighbour_dir), Path(self.test_dir))
    
    @patch("FUR_module_optimized.subprocess.run")
    @patch("FUR_module_optimized.Path.write_text")
    @patch("FUR_module_optimized.Path.stat")
    @patch("FUR_module_optimized.Path.exists")
    @patch('builtins.print')
    @patch('sys.exit')
    def test_furout_empty(self, mock_exit, mock_print, mock_exists, mock_stat, mock_write_text, mock_run):
        # Arrange
        mock_exists.return_value = True  # Pretend file exists
        mock_stat.return_value.st_size = 0  # Pretend file is empty
        mock_run.return_value.returncode = 0  # Pretend subprocess.run is successful

        # Test directory and file name
        source = Path('/some/path')
        outfile_prefix = 'test_prefix'
        option = 'some_option'
        
        # Expected result
        expected_result = f"Fur could not find unique regions or did not run successfully. {source / f'{outfile_prefix}_FUR.db.out.txt'} is empty or does not exist."

        # run the actual function
        run_fur(option, outfile_prefix, source)

        # Assert
        mock_print.assert_called_with(expected_result)

    @patch("FUR_module_optimized.subprocess.run")
    @patch("FUR_module_optimized.Path.write_text")
    @patch("FUR_module_optimized.Path.stat")
    @patch("FUR_module_optimized.Path.exists")
    @patch('builtins.print')
    @patch('sys.exit')
    def test_furout_not_existing(self, mock_exit, mock_print, mock_exists, mock_stat, mock_write_text, mock_run):
        # Arrange
        mock_exists.return_value = False  # Pretend file doesn't exist
        mock_stat.return_value.st_size = 1  # Pretend file is not empty
        mock_run.return_value.returncode = 0  # Pretend subprocess.run is successful

        # Test directory and file name
        source = Path('/some/path')
        outfile_prefix = 'test_prefix'
        option = 'some_option'
        
        # Expected result
        expected_result = f"Fur could not find unique regions or did not run successfully. {source / f'{outfile_prefix}_FUR.db.out.txt'} is empty or does not exist."

        # run the actual function
        run_fur(option, outfile_prefix, source)

        # Assert
        mock_print.assert_called_with(expected_result)
    
    @patch("FUR_module_optimized.subprocess.run")
    @patch('sys.exit')
    def test_FUR_error(self,mock_exit, mock_run):
        mock_run.side_effect=subprocess.CalledProcessError(1, 'fur')
        with self.assertRaises(subprocess.CalledProcessError):
            make_fur_db(Path(self.target_dir), Path(self.neighbour_dir), Path(self.test_dir))

    @patch('FUR_module_optimized.subprocess.run')
    def test_successful_makeFurDb(self, mock_run):
        target_folder = Path('/fake/target')
        neighbour_folder = Path('/fake/neighbour')
        source = Path('/fake/source')

        make_fur_db(target_folder, neighbour_folder, source)
        mock_run.assert_called_once_with(['makeFurDb', '-t', str(target_folder), '-n', str(neighbour_folder), '-d', str(source / 'FUR.db')], check=True)
    
    @patch('FUR_module_optimized.Path.write_text')
    @patch('FUR_module_optimized.Path.stat')
    @patch('FUR_module_optimized.Path.exists')
    @patch('FUR_module_optimized.subprocess.run')
    @patch('FUR_module_optimized.quit')
    def test_run_fur(self, mock_quit, mock_subprocess_run, mock_exists, mock_stat, mock_write_text):
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
        run_fur(option, outfile_prefix, source)
        mock_subprocess_run.assert_called_once_with(['fur', '-d', str(source / 'FUR.db'), '-u'], check=True, capture_output=True, text=True)
        mock_write_text.assert_called_once_with("FUR output", encoding="utf-8")

if __name__ == '__main__':
    unittest.main()
        

