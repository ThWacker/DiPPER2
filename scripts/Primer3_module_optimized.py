#!/usr/bin/python

import os
import sys
import argparse
import subprocess
import re
from datetime import datetime
import shutil
import glob
from pathlib import Path

def quit(message=None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    usage()
    sys.exit(1)

def parse_primers(file_name: str) -> list:
    """Parse the primer penalties from the file and return the top lines."""
    with open(file_name, 'r') as file:
        lines = file.readlines()

    filtered_lines = [line for line in lines if 'PRIMER_PAIR_0_PENALTY' in line]
    split_lines = [re.split(r'\s*=\s*', line.strip()) for line in filtered_lines]

    try:
        sorted_lines = sorted(split_lines, key=lambda x: float(x[1]))
    except ValueError as e:
        quit(f"Error in sorting lines: {e}")

    if not sorted_lines:
        quit("Error: No primers found in file")

    return sorted_lines[:4]

def find_and_return_following_lines_and_target(file_name: str, top_lines: list) -> dict:
    """Find and return following lines and targets for the top primer penalties."""
    found_data = {}

    with open(file_name, 'r') as file:
        lines = file.readlines()

    for top_line in top_lines:
        joined_top_line = '='.join(top_line)
        for i, line in enumerate(lines):
            if joined_top_line in line:
                target_sequence = lines[i-5].strip()
                primer_sequences = [lines[i + j].strip() for j in range(4, 7) if i + j < len(lines)]
                found_data[f"Target_{i}"] = target_sequence
                found_data[f"Primer_{i}"] = primer_sequences

    return found_data

def move_files(source_folder: Path, destination_folder: Path, pattern: str) -> None:
    """Move files matching pattern from source_folder to destination_folder."""
    try:
        for file in source_folder.glob(pattern):
            shutil.move(file, destination_folder / file.name)
            print(f"Moved {file} to {destination_folder}")
    except Exception as e:
        print(f"Error moving files: {e}")

def usage():
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Theresa Wacker T.Wacker2@exeter.ac.uk')
    print('')
    print('Module runs seqkit amplicon, blastx and RRW-Primer-Blast')
    print('')
    print('Usage:')
    print('optional:')
    print('-v/--version for the version')
    print('-o/--outfile_prefix - for outfile prefix [default: date and time in %d-%m-%Y_%Hh%Mmin%Ss_%z format]')
    print('-p/--parameter - string from config file that defines the primer3_core parameters that for fur2prim. The default is: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67')
    print('')
    print('CAUTION: both FUR and primer3_core need to be in path.')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')

def main():
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")

    parser = argparse.ArgumentParser(description='primer3core module')
    parser.add_argument('-p', '--parameter', type=str, default="primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67", 
                        help='string from config file that defines the primer3_core parameters for fur2prim. Default: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67')
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument('-o', '--outfile_prefix', default=dt_string, type=str, help='Outfile prefix. Default is date and time in %%d-%%m-%%Y_%%Hh%%Mmin%%Ss_%%z format')
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which includes results folders from previous steps')

    args = parser.parse_args()

    source_folder = Path(args.folder)

    os.chdir(source_folder)

    try:
        target = next(source_folder.glob('*FUR.db.out.txt'))
    except StopIteration:
        quit("Exception while trying to find the outfile of the FUR_module.py FUR.db.out.txt.")

    param_list = args.parameter.split()
    if len(param_list) != 6:
        quit(f"Error: {args.parameter} does not have 6 elements. You must define all changed values: primMinTm, primOptTm, primMaxTm, inMinTm, inOptTm, and inMaxTm")

    try:
        print("Running fur2prim.")
        out = subprocess.run(['fur2prim'] + param_list + [str(target)], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        quit(f"fur2prim failed: {e}")

    resultf2p = target.with_suffix('.primers.txt')
    resultf2p.write_text(out.stdout.strip())

    if resultf2p.stat().st_size == 0:
        quit(f"{resultf2p} is empty.")

    try:
        print("Running primer3_core.")
        primer3 = subprocess.run(['primer3_core', str(resultf2p)], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        quit(f"primer3_core failed: {e}")

    resultp3 = resultf2p.with_suffix('.primer3_out.txt')
    resultp3.write_text(primer3.stdout.strip())

    if not resultp3.exists() or resultp3.stat().st_size == 0:
        quit(f"Fur could not find unique regions or did not run successfully. {resultp3} does not exist or is empty.")

    top_4_lines = parse_primers(resultp3)
    print("The top 4 lowest primer penalties are:")
    for line in top_4_lines:
        print('\t'.join(line))

    print("Retrieving primers and targets for the lowest penalties...")
    found_data = find_and_return_following_lines_and_target(resultp3, top_4_lines)

    print("The resulting targets and primers are:")
    for key, value in found_data.items():
        print(f"{key}: {value}")
        is_primer = "Primer" in key
        temp_file = f"{args.outfile_prefix}{key}.txt"
        with open(temp_file, 'w') as tf:
            if is_primer:
                for element in value:
                    split_lines = re.split(r'\s*=\s*', element.strip())
                    if len(split_lines) >= 2:
                        tf.write(f">{split_lines[0]}\n{split_lines[1]}\n")
            else:
                split_line = re.split(r'\s*=\s*', value.strip())
                if len(split_line) >= 2:
                    tf.write(f">{split_line[0]}\n{split_line[1]}\n")

    destination_folder_pr = source_folder / "FUR.P3.PRIMERS"
    destination_folder_tar = source_folder / "FUR.P3.TARGETS"
    destination_folder_pr.mkdir(parents=True, exist_ok=True)
    destination_folder_tar.mkdir(parents=True, exist_ok=True)

    move_files(source_folder, destination_folder_tar, "*_target.txt")
    move_files(source_folder, destination_folder_pr, "*_primer.txt")

    print('Primer3_module.py ran to completion: exit status 0')

if __name__ == '__main__':
    main()
