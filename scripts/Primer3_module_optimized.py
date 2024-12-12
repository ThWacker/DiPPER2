#!/usr/bin/python

import os
import sys
import argparse
import subprocess
import re
from datetime import datetime
import shutil
from pathlib import Path
from fur2primer3 import remap_keys, write_result, args_to_dict
from logging_handler import Logger


def parse_primers(file_name: str, logger: Logger) -> list:
    """
    Parse the primer penalties from the file and return the top lines.

    Args:
        file_name (str): the file name of the primer3 output

    Returns:
        sorted_lines(list): the first 4 elements of a list of penalty sorted incrementally
    """
    with open(file_name, "r", encoding="utf-8") as file:
        lines = file.readlines()

    # if line with PRIMER_PAIR_0_PENALTY is found, put in list
    filtered_lines = [line for line in lines if "PRIMER_PAIR_0_PENALTY" in line]

    # This line processes each line in 'filtered_lines' by stripping whitespace and splitting 
    # it into two parts at the first occurrence of the '=' symbol. The result is a list of lists,
    # where each sublist contains the parts of the line before and after the '='. The regular 
    # expression 'r"\s*=\s*"' ensures that any surrounding whitespace around the '=' is ignored.
    split_lines = [re.split(r"\s*=\s*", line.strip()) for line in filtered_lines]

    try:
        # Attempt to sort 'split_lines' by the second element in each sublist (x[1]), converting it to a float
        sorted_lines = sorted(split_lines, key=lambda x: float(x[1]))
    except ValueError as e:
        # If a ValueError occurs during sorting (e.g., non-numeric data), log the error with exception info
        logger.exception("Error in sorting lines", exc_info=1)
        # Raise a new ValueError, preserving the original exception details
        raise ValueError(f"Error in sorting lines: {e}") from e

    # If the sorted list is empty, log an error indicating no primers were found
    if not sorted_lines:
        logger.error("Error: No primers found in file")
        # Raise an exception because no valid data was found to return
        raise Exception("Error: No primers found in file")

    # Return the first 4 items from the sorted list
    return sorted_lines[:4]



def find_and_return_following_lines_and_target(
    file_name: str, top_lines: list, qpcr: str
) -> dict:
    """
    Find the primer penalty and return following lines aka primers and targets for the top primer penalties.

    Args:
        file_name(str): the primer3 output file
        top_lines(list): list with the 4 lowest primer penalties
        qPCR(str): yes/no toggle to decide if internal probe is also returned or not

    Returns:
        found_data(dict): a dictionary containing primers and targets, as well as data (Tm etc.) of primers
    """
    #initialize dictionary
    found_data = {}

    with open(file_name, "r", encoding="utf-8") as file:
        lines = file.readlines()

    # Iterate over each 'top_line' in 'top_lines'
    for top_line in top_lines:
        # Join the elements of 'top_line' with an equals sign ('=') between them
        joined_top_line = "=".join(top_line)
        
        # Check if 'qpcr' is equal to 'y' (yes) to determine the flow of processing
        if qpcr == "y":
            # Iterate over the 'lines' list and enumerate them to get both index and content
            for i, line in enumerate(lines):
                # If the 'joined_top_line' is found in the current line
                if joined_top_line in line:
                    # Retrieve the target sequence from 5 lines above the current one
                    target_sequence = lines[i - 5].strip()
                    # Extract primer sequences from lines 4 to 6 after the current line
                    primer_sequences = [
                        lines[i + j].strip() for j in range(4, 7) if i + j < len(lines)
                    ]
                    # Extract additional primer data from lines 7 to 29 after the current line
                    primer_data = [
                        lines[i + j].strip() for j in range(7, 30) if i + j < len(lines)
                    ]
                    # Store the found data in the 'found_data' dictionary with dynamic keys
                    found_data[f"Target_{i}"] = target_sequence
                    found_data[f"Primer_{i}"] = primer_sequences
                    found_data[f"Data_primer_{i}"] = primer_data
        else:
            # If 'qpcr' is not 'y', follow a different data extraction process
            for i, line in enumerate(lines):
                # If the 'joined_top_line' is found in the current line
                if joined_top_line in line:
                    # Retrieve the target sequence from 5 lines above the current one
                    target_sequence = lines[i - 5].strip()
                    # Extract primer sequences from lines 3 to 4 after the current line
                    primer_sequences = [
                        lines[i + j].strip() for j in range(3, 5) if i + j < len(lines)
                    ]
                    # Append a placeholder for a missing primer sequence when 'qpcr' is not 'y'
                    primer_sequences.append("PRIMER_INTERNAL_0_SEQUENCE=NA")
                    # Extract additional primer data from lines 5 to 21 after the current line
                    primer_data = [
                        lines[i + j].strip() for j in range(5, 22) if i + j < len(lines)
                    ]
                    # Store the found data in the 'found_data' dictionary with dynamic keys
                    found_data[f"Target_{i}"] = target_sequence
                    found_data[f"Primer_{i}"] = primer_sequences
                    found_data[f"Data_primer_{i}"] = primer_data

    # Return the 'found_data' dictionary containing all the extracted information
    return found_data


def move_files(source_folder: Path, destination_folder: Path, pattern: str, logger: Logger) -> None:
    """
    Move files matching pattern from source_folder to destination_folder.

    Args:
        source_folder(Path): Path object containing the files that are supposed to be moved into the destination folder
        destination_folder(Path): Path object the files are supposed to be shifted to
        pattern(str): pattern to be matched to identify which files from the source folder go into the destination folder
        
    Returns:
        None
    """
    # glob returns an iterator of Path objects that match the pattern
    try:
        for file in source_folder.glob(pattern):
            shutil.move(file, destination_folder / file.name)
            logger.info(f"Moved {file} to {destination_folder}")
    # throw error if no files are found
    except FileNotFoundError as e:
        logger.exception("Error moving files", exc_info=1)
        raise FileNotFoundError(f"Error moving files: {e}") from e


def main():
    """
    main function doing too much at the moment.
    tbc
    """
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")

    parser = argparse.ArgumentParser(
        prog="primer3core module of DiPPER2",
        description="Please ensure that primer3 is found in path",
        epilog="Theresa Wacker T.Wacker2@exeter.ac.uk",
    )
    parser.add_argument(
        "-p",
        "--parameter",
        type=str,
        default="primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200 Oligo=1",
        help="string from config file that defines the primer3_core parameters for fur2prim. Default: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200",
    )
    parser.add_argument("-V", "--version", action="version", version="%(prog)s 0.1.0-beta")
    parser.add_argument(
        "-o",
        "--outfile_prefix",
        default=dt_string,
        type=str,
        help="Outfile prefix. Default is date and time in d-m-y-h-m-s-tz format",
    )
    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        required=True,
        help="Results folder name, which includes results folders from previous steps",
    )
    parser.add_argument(
        "-q",
        "--qpcr",
        type=str,
        default="n",
        help="Toggle whether to expect an internal probe sequence (qPCR) or not. Default is no(n).",
    )
    parser.add_argument('-v', '--verbose', action="store_true", help="increase logging verbosity" )

    args = parser.parse_args()

    # set source folder
    source_folder = Path(args.folder)

    # configures the logger
    module_name = Path(__file__).name
    logger_instance = Logger(module_name, source_folder, args.verbose)
    logger = logger_instance.get_logger()

    # info
    logger.info("Starting the Primer3 module...")

    # change into the folder
    os.chdir(source_folder)

    # looks for all FUR.db.out.txt files, but should only find one. However, it should not find none.
    try:
        target = next(source_folder.glob("*FUR.db.out.txt"))
    except StopIteration as e:
        logger.error("Could not find {target}", exc_info=1)
        raise StopIteration(
            "Exception while trying to find the outfile of the FUR_module.py FUR.db.out.txt."
        ) from e

    # check if there is the right number of parameters that are supposed to be fed into fur2prim
    param_list = args.parameter.split()
    if len(param_list) != 9:
        logger.error(
            f"Error: {args.parameter} does not have 9 elements. You must define all changed values: primMinTm, primOptTm, primMaxTm, inMinTm, inOptTm, inMaxTm, prodMinSize=100, prodMaxSize=200 and Oligo=1"
        )
        sys.exit()

    # convert fur output into primer3 compatible output using the functions from the fur2primer3 script
    logger.info(
        "Convert the FUR output to Primer3 compatible output using fur2primer3 functions..."
    )
    try:
        # parse parameters (dictionary)
        params = args_to_dict(args.parameter)
        # remap/ convert parameter keys to Primer3 conventions
        form_param = remap_keys(params)
        # write Primer3 compatible file.
        write_result(Path(target), form_param)

    # Error handling
    except FileNotFoundError as e:
        logger.error("Could not find {target}", exc_info=1)
        raise FileNotFoundError(f"Could not find {target}") from e

    except ValueError as e:
        logger.exception(f"Input or output values are not as expected: {e}", exc_info=1)
        raise ValueError(f"Input or output values are not as expected: {e}") from e

    except TypeError as e:
        logger.exception("Wrong type", exc_info=1)
        raise TypeError(f"Wrong type: {e}") from e
    except OSError as e:
        logger.exception(f"Could not write or open file: {e}", exc_info=1)
        raise OSError(f"Could not write or open file: {e}") from e
    except Exception as e:
        logger.exception(
            f"Unknown exception occured while running fur2primer3 functions: {e}",
            exc_info=1,
        )
        raise Exception(
            f"Unknown exception occured while running fur2primer3 functions: {e}"
        ) from e

    # define results file
    resultf2p = target.with_suffix(".primers.txt")

    # check if the results file is empty, if so sys.exit 1 with message (might not raise correct exception)
    if resultf2p.stat().st_size == 0:
        logger.error(f"{resultf2p} is empty.")
        sys.exit()

    # try running primer3_core, write it to file, check if the file exists and/or is empty. If so,
    try:
        logger.info("Running primer3_core.")
        primer3 = subprocess.run(
            ["primer3_core", str(resultf2p)], check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        logger.error("Primer 3 did not run successfully.", exc_info=1)
        raise subprocess.CalledProcessError(
            e.returncode, e.cmd, output=e.output, stderr=e.stderr
        ) from e

    # write the resulyts of the primer3 run to a file
    resultp3 = resultf2p.with_suffix(".primer3_out.txt")
    resultp3.write_text(primer3.stdout.strip(), encoding="utf-8")

    # complain about no results and raise exception
    if not resultp3.exists() or resultp3.stat().st_size == 0:
        logger.error(
            "Primer3 did not manage to generate primers or did not run successfully. {resultp3} does not exist or is empty.",
            exc_info=1,
        )
        raise FileExistsError(
            f"Primer3 did not manage to generate primers or did not run successfully. {resultp3} does not exist or is empty."
        )

    # parse primer3 output and find the 4 lowest penalty scores.
    top_4_lines = parse_primers(resultp3, logger)
    logger.info("The top 4 lowest primer penalties are:")
    for line in top_4_lines:
        logger.info("\t".join(line))

    # Inform the user that we are not extracting primers and targets for the lowest penalties
    logger.info("Retrieving primers and targets for the lowest penalties...")

    # get a dictionary that contains the primers and the data. Primers and targets are identified by the line in the primer3 output their corresponding lowest penalty score is found at
    found_data = find_and_return_following_lines_and_target(
        resultp3, top_4_lines, args.qpcr
    )

    logger.info("The resulting targets and primers are:")

    # Iterate over the items (key-value pairs) in the 'found_data' dictionary
    for key, value in found_data.items():
        # Log the current key-value pair
        logger.info(f"{key}: {value}")
        
        # Check if the current key is related to "Primer" by checking if the string "Primer" is in the key
        is_primer = "Primer" in key
        
        # Define the name of the temporary file where the data will be written
        temp_file = f"{args.outfile_prefix}{key}.txt"
        
        # Open the temporary file for writing with UTF-8 encoding
        with open(temp_file, "w", encoding="utf-8") as tf:
            # If the current key is related to primer data
            if is_primer:
                # Iterate through each element in the value (which should be a list of primer data)
                for element in value:
                    # Split the element by '=' and remove unnecessary whitespace around it
                    split_lines = re.split(r"\s*=\s*", element.strip())
                    
                    # If the split results in at least two parts (a name and a sequence)
                    if len(split_lines) >= 2:
                        # Write the primer data to the file in the format ">name\nsequence\n"
                        tf.write(f">{split_lines[0]}\n{split_lines[1]}\n")
            
            # If the key contains "Data" (indicating it's data-related)
            elif "Data" in key:
                # Write all the data (joined with newlines) to the file
                tf.write("\n".join(value))
            
            # For other cases (non-primer, non-data)
            else:
                # Split the value by '=' and remove whitespace
                split_line = re.split(r"\s*=\s*", value.strip())
                
                # If the split results in at least two parts (a name and a sequence)
                if len(split_line) >= 2:
                    # Write the name and sequence to the file in the format ">name\nsequence\n"
                    tf.write(f">{split_line[0]}\n{split_line[1]}\n")

    # Define the paths for the destination folders where files will be moved
    destination_folder_pr = source_folder / "FUR.P3.PRIMERS"
    destination_folder_tar = source_folder / "FUR.P3.TARGETS"
    destination_folder_data = destination_folder_pr / "primer_data"

    # Create the destination folders if they don't already exist
    destination_folder_pr.mkdir(parents=True, exist_ok=True)
    destination_folder_tar.mkdir(parents=True, exist_ok=True)
    destination_folder_data.mkdir(parents=True, exist_ok=True)

    # Move files with specific patterns from the source folder to the appropriate destination folders
    move_files(source_folder, destination_folder_tar, "*Target*.txt", logger)
    move_files(source_folder, destination_folder_pr, "*Primer*.txt", logger)
    move_files(source_folder, destination_folder_data, "*Data*.txt", logger)

    # Log that the script has completed successfully
    logger.info("Primer3_module.py ran to completion: exit status 0")


if __name__ == "__main__":
    main()
