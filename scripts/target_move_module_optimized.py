#! /usr/bin/python

import sys
import argparse
import shutil
import subprocess
from pathlib import Path
from logging_handler import Logger

def parse_list(targets_file: Path, logger: Logger)-> list:
    """Parse the list of targets from a file.
    Args:
        targets_file(Path): a file with a list of targets
    Returns:
        A list of Accessions with the version number stripped off
    Raises:
        FileNotFoundError
    """
    try:
        with open(targets_file, 'r',encoding="utf-8") as file:
            return [line.split('\t')[0].split('.')[0] for line in file]
    except FileNotFoundError as e:
        logger.exception(f"Error: Target file '{targets_file}' not found!: {e}")
        raise FileNotFoundError(f"Error: Target file '{targets_file}' not found!: {e}") from e
        
def process_target_files(list_targets: list, source_folder:Path, fur_target_folder:Path, logger: Logger):
    """Process each target file based on the list of targets.
    Args:
        list_targets(list): The list with the target accessions
        source_folder(Path): a path object of the folder to rule them all (all results of this DiPPER2 run are in this folder)
        fur_target_folder(Path): the target folder the accessions that are targets will be copiedn in. 
        logger(Logger): the logger
    Returns:
        None
    Raises:
        RunTimeError
    """
    for accession in list_targets:
        accession_pattern = accession.strip() + "*"
        try:
            result = subprocess.run(['find', str(source_folder), '-maxdepth', '1', '-type', 'f', '-iname', accession_pattern], capture_output=True, text=True, check=True)
            filename = result.stdout.strip()
            test_list=filename.split(sep="\n")
            if len(test_list) > 1:
                logger.error("Make sure that there is only one assembly per accession. Please delete the other assemblies and only keep the highest quality one.")
                sys.exit()
            if filename:
                if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
                    target_file = Path(filename)
                    logger.debug(f"{target_file} will be copied to {fur_target_folder}")
                    shutil.copy(target_file, fur_target_folder)
                else:
                    logger.debug(f"{filename} is not a fasta file. Skipping.")
                    continue
            else:
                logger.warning(f"{accession_pattern} not found! Are you in the folder with the assemblies? Please double-check!")
        except subprocess.CalledProcessError as e:
            logger.exception(f"Error finding file {accession_pattern}: {e}")
            raise RuntimeError(f"Error finding {accession_pattern} or copying the file: {e}") from e

def move_remaining_files_to_neighbour(source_folder: Path, fur_target_folder: Path, fur_neighbour_folder: Path, logger:Logger):
    """
    Move remaining files to the neighbour folder.
    Args:
        source_folder(Path): the folder with all results of this run
        fur_target_folder(Path): folder with the target accession assemblies
        fur_neighbour_folder(Path): the folder with the neighbour accession assemblies
        logger(Logger): the logger
    Returns:
        None
    """
    for file_path in source_folder.iterdir():
        if file_path.is_file() and not (fur_target_folder / file_path.name).exists():
            if file_path.name.endswith(".fasta") or file_path.name.endswith(".fa") or file_path.name.endswith (".fna"):
                shutil.copy(file_path, fur_neighbour_folder)
                logger.debug(f'Copied to {fur_neighbour_folder}: {file_path.name}')
            else:
                logger.warning(f"{file_path.name} is not a fasta file. Skipping.")
                continue
        else:
            logger.debug(f'File {file_path.name} already exists in {fur_target_folder}, skipping.')

def main():
    parser = argparse.ArgumentParser(prog='Dipper2',description='Target and Neighbour Folder Sorting', epilog="Bugs, suggestions, criticism, cake, capybaras and praise to t.wacker2@exeter.ac.uk")
    parser.add_argument('-t', '--target', type=str, required=True, help='Please add target list')
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which will include the FUR.target and FUR.neighbour subfolders')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.0-beta')

    args = parser.parse_args()

    # Define source folder
    source_folder = Path.cwd()
    
    # configures the logger
    module_name = Path(__file__).name
    logger_instance = Logger(module_name, source_folder, args.verbose)
    logger = logger_instance.get_logger()

    #parse file with accessions to list
    list_targets = parse_list(args.target, logger)
    dist_folder = Path(args.folder)
    
    # Make directories in the parent folder
    fur_target_folder = dist_folder / "FUR.target"
    fur_neighbour_folder = dist_folder / "FUR.neighbour"
    fur_target_folder.mkdir(parents=True, exist_ok=True)
    fur_neighbour_folder.mkdir(parents=True, exist_ok=True)

    # Process each target file
    process_target_files(list_targets, source_folder, fur_target_folder, logger)

    # Move remaining files to the neighbour folder
    move_remaining_files_to_neighbour(source_folder, fur_target_folder, fur_neighbour_folder, logger)

if __name__ == '__main__':
    main()
