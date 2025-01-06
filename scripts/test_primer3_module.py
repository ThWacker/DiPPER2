#!/usr/bin/python

import unittest
from pathlib import Path
import tempfile
import shutil
import subprocess
from unittest.mock import patch, MagicMock

# Assuming all your functions are imported from your script
from Primer3_module_optimized import (
    move_files,
    parse_primers,
    find_and_return_following_lines_and_target,
    main,
)


class TestMoveFolders(unittest.TestCase):
    def setUp(self):
        # Setup temporary directories and files for testing
        self.source_dir = tempfile.mkdtemp()
        self.destination_dir = tempfile.mkdtemp()
        self.file1 = Path(self.source_dir) / "file1.fasta"
        with self.file1.open("w") as f:
            f.write(">MockAssembly GCA123\nATGCGCTAA")

        self.file2 = Path(self.source_dir) / "file2.fasta"
        with self.file2.open("w") as f:
            f.write(">ThisIsATest\n atcgatcgatcg\n")
            self.file3 = Path(self.source_dir) / "file3.fasta"
        with self.file3.open("w") as f:
            f.write(">HelloWorld")

    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.source_dir)
        shutil.rmtree(self.destination_dir)

    @patch("Primer3_module_optimized.Logger")
    @patch("Primer3_module_optimized.shutil.move")
    @patch("Primer3_module_optimized.Path.glob")
    def testMoveFoldersSuccess(self, mock_glob, mock_move, mock_logger_class):
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Mock the files returned by glob
        mock_glob.return_value = [Path(self.file1), Path(self.file2), Path(self.file3)]

        # Call the function being tested
        move_files(
            Path(self.source_dir),
            Path(self.destination_dir),
            "*.fasta",
            mock_logger_instance,
        )

        # Check if shutil.move is called correctly
        self.assertEqual(mock_move.call_count, 3)
        mock_move.assert_any_call(
            Path(self.file1), Path(self.destination_dir, "file1.fasta")
        )
        mock_move.assert_any_call(
            Path(self.file2), Path(self.destination_dir, "file2.fasta")
        )
        mock_move.assert_any_call(
            Path(self.file3), Path(self.destination_dir, "file3.fasta")
        )

    @patch("Primer3_module_optimized.shutil.move")
    @patch("Primer3_module_optimized.Path.glob")
    @patch("Primer3_module_optimized.Logger")
    def testMoveFoldersLogger(self, mock_logger_class, mock_glob, mock_move):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Mock the files returned by glob
        mock_glob.return_value = [Path(self.file1)]

        # Call the function being tested
        move_files(
            Path(self.source_dir),
            Path(self.destination_dir),
            "*.fasta",
            mock_logger_instance,
        )

        # Check if logger.info was called with the info message
        mock_logger_instance.info.assert_called_once()
        self.assertTrue(
            f"Moved {Path(self.file1)} to {Path(self.destination_dir)}"
            in mock_logger_instance.info.call_args[0][0]
        )

    @patch("Primer3_module_optimized.shutil.move")
    @patch("Primer3_module_optimized.Path.glob")
    @patch("Primer3_module_optimized.Logger")
    def testMoveFolderFailureNoGlob(self, mock_logger_class, mock_glob, mock_move):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Mock the files returned by glob
        mock_glob.return_value = []  # No files found

        # Assert that FileNotFoundError is raised
        with self.assertRaises(FileNotFoundError):
            logger = mock_logger_instance
            # Call the function being tested
            move_files(
                Path(self.source_dir), Path(self.destination_dir), "*.fasta", logger
            )

        # Check if logger.exception was called with the info message
        mock_logger_instance.exception.assert_called_once()
        self.assertTrue(
            "No files matching *.fasta found in"
            in mock_logger_instance.exception.call_args[0][0]
        )

    @patch("Primer3_module_optimized.shutil.move")
    @patch("Primer3_module_optimized.Path.glob")
    @patch("Primer3_module_optimized.Logger")
    def testMoveFolderFailureFileNotFound(
        self, mock_logger_class, mock_glob, mock_move
    ):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Configure shutil.move to raise FileNotFoundError
        mock_move.side_effect = FileNotFoundError("Mocked FileNotFoundError")

        # Mock source and destination with explicit names and paths. Slight redundance.
        source_dir = Path("/mock/DieQuelle")
        destination_dir = Path("/mock/FinalDestination")
        pattern = "*.fasta"

        # Mock source folder containing files
        with patch.object(
            Path, "glob", return_value=[Path("/mock/DieQuelle/file1.fasta")]
        ):
            # Assert that FileNotFoundError is raised
            with self.assertRaises(FileNotFoundError):
                move_files(source_dir, destination_dir, pattern, mock_logger_instance)

        # Check if the logger.exception was called
        mock_logger_instance.exception.assert_called_once_with(
            "Error moving files", exc_info=1
        )

        # Check that shutil.move was called
        mock_move.assert_called_once_with(
            Path("/mock/DieQuelle/file1.fasta"),
            Path("/mock/FinalDestination/file1.fasta"),
        )

    @patch("Primer3_module_optimized.shutil.move")
    @patch("Primer3_module_optimized.Path.glob")
    @patch("Primer3_module_optimized.Logger")
    def testMoveFolderFailureOSError(self, mock_logger_class, mock_glob, mock_move):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # Configure shutil.move to raise OSError
        mock_move.side_effect = OSError("Oh No an OSError")

        # Mock source and destination with explicit names and paths. Slight redundance.
        source_dir = Path("/mock/DieQuelle")
        destination_dir = Path("/mock/FinalDestination")
        pattern = "*.fasta"

        # Mock source folder containing files
        with patch.object(
            Path, "glob", return_value=[Path("/mock/DieQuelle/file1.fasta")]
        ):
            # Assert that FileNotFoundError is raised
            with self.assertRaises(OSError):
                move_files(source_dir, destination_dir, pattern, mock_logger_instance)

        # Check if the logger.exception was called
        mock_logger_instance.exception.assert_called_once_with(
            "OS error during moving files with shutil.move: Oh No an OSError"
        )

        # Check that shutil.move was called
        mock_move.assert_called_once_with(
            Path("/mock/DieQuelle/file1.fasta"),
            Path("/mock/FinalDestination/file1.fasta"),
        )


class testParsePrimers(unittest.TestCase):

    def setUp(self):
        self.source_dir = tempfile.mkdtemp()
        self.filePrimer3 = Path(self.source_dir) / "primer3.txt"
        with self.filePrimer3.open("w") as f:
            f.write(
                "PRIMER_PAIR_0_PENALTY=0.456789\nPRIMER_PAIR_0_PENALTY=0.56789\nPRIMER_PAIR_0_PENALTY=0.12345\nPRIMER_PAIR_0_PENALTY=0.891011"
            )
        self.filePrimer_VE = Path(self.source_dir) / "primer_novalues.txt"
        with self.filePrimer_VE.open("w") as f:
            f.write(
                "PRIMER_PAIR_0_PENALTY=NOT\nPRIMER_PAIR_0_PENALTY=A\nPRIMER_PAIR_0_PENALTY=NUMERICAL\nPRIMER_PAIR_0_PENALTY=VALUE"
            )
        
        self.filePrimer_empty = Path(self.source_dir) / "primer_empty.txt"
        with self.filePrimer_empty.open("w") as f:
            f.write(
                "\s"
            )

    def tearDown(self):
        # Remove temporary directory and files after testing
        shutil.rmtree(self.source_dir)

    @patch("Primer3_module_optimized.Logger")
    def testParsePrimersSuccess(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # run function on mocked file, compare return
        returned = parse_primers(Path(self.filePrimer3), mock_logger_instance)

        # expected result
        expect = [
            ["PRIMER_PAIR_0_PENALTY", "0.12345"],
            ["PRIMER_PAIR_0_PENALTY", "0.456789"],
            ["PRIMER_PAIR_0_PENALTY", "0.56789"],
            ["PRIMER_PAIR_0_PENALTY", "0.891011"],
        ]
        # Assert the expected output
        self.assertEqual(returned, expect)

    @patch("Primer3_module_optimized.Logger")
    def testParsePrimersValueError(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # run function on mocked file, compare return
        with self.assertRaises(ValueError):
            parse_primers(Path(self.filePrimer_VE), mock_logger_instance)

        # Check if logger.exception was called with the info message
        mock_logger_instance.exception.assert_called_once()
        self.assertTrue(
            "Error in sorting lines" in mock_logger_instance.exception.call_args[0][0]
        )
    
    @patch("Primer3_module_optimized.Logger")
    def testParsePrimersException(self, mock_logger_class):
        # Create a mock logger instance
        mock_logger_instance = MagicMock()
        mock_logger_class.return_value.get_logger.return_value = mock_logger_instance

        # run function on mocked file, compare return
        with self.assertRaises(Exception):
            parse_primers(Path(self.filePrimer_empty), mock_logger_instance)

        # Check if logger.error was called with the info message
        mock_logger_instance.error.assert_called_once()
        self.assertTrue(
            "Error: No primers found in file" in mock_logger_instance.error.call_args[0][0]
        )
class test_find_and_return_following_lines_and_target(unittest.TestCase):
    def test(self):
        print("placeholder")

if __name__ == "__main__":
    unittest.main()
