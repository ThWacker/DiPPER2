#!/usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
import subprocess
from unittest.mock import patch, MagicMock

# Assuming all your functions are imported from your script
from Primer3_module_optimized import (quit,move_files,parse_primers,find_and_return_following_lines_and_target, main, usage)

class TestHelperFunctions(unittest.TestCase):
        def setUp(self):
        # Setup temporary directories and files for testing
        self.source_dir = tempfile.mkdtemp()
        self.dist_dir
        self.file1 = Path(self.source_dir) / 'file1.fasta'
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
    def test_quit_program(self):
        with self.assertRaises(SystemExit):
            quit("More tests")
    
    def test_usage(self):
        with patch('builtins.print') as mocked_print:
            usage()
            self.assertTrue(mocked_print.called) 
    
    @patch("Primer3_module_optimized.shutil.move")
    @patch("Primer3_module_optimized.Path.glob")
    def test_move_folders_success(self, mock_glob, mock_move):
        # Mock the files returned by glob
        mock_glob.return_value = [Path("file1.fasta"), Path("file2.fasta"), Path("file3.fasta")]
        
        # Create mocked source and destination folders
        mocked_source = MagicMock(spec=Path)
        mock_destination = MagicMock(spec=Path)
        
        # Set up __truediv__ to return actual paths, not more MagicMock objects
        mock_destination.__truediv__.side_effect = lambda x: Path(mock_destination, x)
        mocked_source.__truediv__.side_effect = lambda x: Path(mocked_source, x)

        # Call the function being tested
        move_files(mocked_source, mock_destination, "*.fasta")
        
        # Check if shutil.move is called correctly
        #self.assertEqual(mock_move.call_count, 3)
        mock_move.assert_any_call(Path("file1.fasta"), Path(mock_destination, "file1.fasta"))
        mock_move.assert_any_call(Path("file2.fasta"), Path(mock_destination, "file2.fasta"))
        mock_move.assert_any_call(Path("file3.fasta"), Path(mock_destination, "file3.fasta"))


        
if __name__ == '__main__':
    unittest.main()
        

