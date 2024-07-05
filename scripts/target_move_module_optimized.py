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
    print('$0 --target/-t <list of targets in tab delimited file> -f/--folder <the folder name (required, target and neighbour folders will be located in there)>')
    print('optional: -v/--version for the version')
    print('CAUTION: needs to be in the folder with the assemblies to be sorted in target and neighbour')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')

def parse_list(targets_file):
    """Parse the list of targets from a file."""
    try:
        with open(targets_file, 'r') as file:
            return [line.split('\t')[0].split('.')[0] for line in file]
    except FileNotFoundError:
        quit(f"Error: Target file '{targets_file}' not found!")

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
    for accession in list_targets:
        accession_pattern = accession.strip() + "*"
        try:
            result = subprocess.run(['find', str(source_folder), '-maxdepth', '1', '-type', 'f', '-iname', accession_pattern], capture_output=True, text=True, check=True)
            filename = result.stdout.strip()
            if filename:
                target_file = Path(filename)
                print(f"{target_file} will be copied to {fur_target_folder}")
                shutil.copy(target_file, fur_target_folder)
            else:
                print(f"{accession_pattern} not found! Are you in the folder with the assemblies? Please double-check!")
        except subprocess.CalledProcessError as e:
            print(f"Error finding file {accession_pattern}: {e}")

    # Move remaining files to the neighbour folder
    for file_path in source_folder.iterdir():
        if file_path.is_file() and not (fur_target_folder / file_path.name).exists():
            shutil.copy(file_path, fur_neighbour_folder)
            print(f'Copied to {fur_neighbour_folder}: {file_path.name}')
        else:
            print(f'File {file_path.name} already exists in {fur_neighbour_folder}, skipping.')

if __name__ == '__main__':
    main()
