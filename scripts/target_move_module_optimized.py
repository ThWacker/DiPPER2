#! /usr/bin/python

import sys
import argparse
import shutil
import subprocess
from pathlib import Path

def quit(message=None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    usage()
    sys.exit(1)

def usage():
    """Display usage information."""
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Theresa Wacker T.Wacker2@exeter.ac.uk')
    print('')
    print('Module/ Script that moves targets and neighbours into folders called FUR.target and FUR.neighbour')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-f/--folder - results folder name, which includes results folders from previous steps')
    print('$0 --target/-t <list of targets in tab delimited file> -f/--folder <the folder name (required, target and neighbour folders will be located in there)>')
    print('optional: -v/--version for the version')
    print('CAUTION: needs to be in the folder with the assemblies to be sorted in target and neighbour')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')

def parse_list(targets_file):
    """Parse the list of targets from a file."""
    try:
        with open(targets_file, 'r',encoding="utf-8") as file:
            return [line.split('\t')[0].split('.')[0] for line in file]
    except FileNotFoundError as e:
        quit(f"Error: Target file '{targets_file}' not found!")
        raise e
        
def process_target_files(list_targets, source_folder, fur_target_folder):
    """Process each target file based on the list of targets."""
    for accession in list_targets:
        accession_pattern = accession.strip() + "*"
        try:
            result = subprocess.run(['find', str(source_folder), '-maxdepth', '1', '-type', 'f', '-iname', accession_pattern], capture_output=True, text=True, check=True)
            filename = result.stdout.strip()
            test_list=filename.split(sep="\n")
            if len(test_list) > 1:
                quit("Make sure that there is only one assembly per accession. Please delete the other assemblies and only keep the highest quality one.")
            if filename:
                if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
                    target_file = Path(filename)
                    print(f"{target_file} will be copied to {fur_target_folder}")
                    shutil.copy(target_file, fur_target_folder)
                else:
                    print(f"{filename} is not a fasta file. Skipping.")
                    continue
            else:
                print(f"{accession_pattern} not found! Are you in the folder with the assemblies? Please double-check!")
        except subprocess.CalledProcessError as e:
            print(f"Error finding file {accession_pattern}: {e}")
            raise e

def move_remaining_files_to_neighbour(source_folder, fur_target_folder, fur_neighbour_folder):
    """Move remaining files to the neighbour folder."""
    for file_path in source_folder.iterdir():
        if file_path.is_file() and not (fur_target_folder / file_path.name).exists():
            if file_path.name.endswith(".fasta") or file_path.name.endswith(".fa") or file_path.name.endswith (".fna"):
                shutil.copy(file_path, fur_neighbour_folder)
                print(f'Copied to {fur_neighbour_folder}: {file_path.name}')
            else:
                print(f"{file_path.name} is not a fasta file. Skipping.")
                continue
        else:
            print(f'File {file_path.name} already exists in {fur_target_folder}, skipping.')

def main():
    parser = argparse.ArgumentParser(description='Target and Neighbour Folder Sorting')
    parser.add_argument('-t', '--target', type=str, required=True, help='Please add target list')
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which will include the FUR.target and FUR.neighbour subfolders')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')

    args = parser.parse_args()

    list_targets = parse_list(args.target)
    dist_folder = Path(args.folder)

    # Make directories in the parent folder
    fur_target_folder = dist_folder / "FUR.target"
    fur_neighbour_folder = dist_folder / "FUR.neighbour"
    fur_target_folder.mkdir(parents=True, exist_ok=True)
    fur_neighbour_folder.mkdir(parents=True, exist_ok=True)

    # Define source folder
    source_folder = Path.cwd()
    
    # Process each target file
    process_target_files(list_targets, source_folder, fur_target_folder)

    # Move remaining files to the neighbour folder
    move_remaining_files_to_neighbour(source_folder, fur_target_folder, fur_neighbour_folder)

if __name__ == '__main__':
    main()
