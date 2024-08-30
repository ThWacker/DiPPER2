#!/usr/bin/python

import os
import sys
import argparse
import subprocess
import re
from datetime import datetime
import shutil
from pathlib import Path
from fur2primer3 import (remap_keys,write_result, args_to_dict)

def quit(message=None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    usage()
    sys.exit(1)

def parse_primers(file_name: str) -> list:
    """Parse the primer penalties from the file and return the top lines."""
    with open(file_name, 'r',encoding="utf-8") as file:
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

    with open(file_name, 'r',encoding="utf-8") as file:
        lines = file.readlines()

    for top_line in top_lines:
        joined_top_line = '='.join(top_line)
        for i, line in enumerate(lines):
            if joined_top_line in line:
                target_sequence = lines[i-5].strip()
                primer_sequences = [lines[i + j].strip() for j in range(4, 7) if i + j < len(lines)]
                primer_data = [lines[i + j].strip() for j in range(8, 30) if i + j < len(lines)]
                found_data[f"Target_{i}"] = target_sequence
                found_data[f"Primer_{i}"] = primer_sequences
                found_data[f"Data_primer_{i}"]=primer_data

    return found_data

def move_files(source_folder: Path, destination_folder: Path, pattern: str) -> None:
    """Move files matching pattern from source_folder to destination_folder."""
    try:
        for file in source_folder.glob(pattern):
            shutil.move(file, destination_folder / file.name)
            #print(f"Moved {file} to {destination_folder}")
    except FileNotFoundError as e:
        quit(f"Error moving files: {e}")
        raise e

def usage():
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Theresa Wacker T.Wacker2@exeter.ac.uk')
    print('')
    print('Module runs primer3_core and sorts targets, primers and data in folders')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-f/--folder - results folder name, which includes results folders from previous steps')
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
    parser.add_argument('-p', '--parameter', type=str, default="primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200 Oligo=1", 
                        help='string from config file that defines the primer3_core parameters for fur2prim. Default: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200')
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument('-o', '--outfile_prefix', default=dt_string, type=str, help='Outfile prefix. Default is date and time in d-m-y-h-m-s-tz format')
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which includes results folders from previous steps')

    args = parser.parse_args()

    source_folder = Path(args.folder)

    #change into the folder
    os.chdir(source_folder)

    # looks for all FUR.db.out.txt files, but should only find one. However, it should not find none.
    try:
        target = next(source_folder.glob('*FUR.db.out.txt'))
    except StopIteration:
        quit("Exception while trying to find the outfile of the FUR_module.py FUR.db.out.txt.")
        raise
    
    # check if there is the right number of parameters that are supposed to be fed into fur2prim
    param_list = args.parameter.split()
    if len(param_list) != 9:
        quit(f"Error: {args.parameter} does not have 9 elements. You must define all changed values: primMinTm, primOptTm, primMaxTm, inMinTm, inOptTm, inMaxTm, prodMinSize=100, prodMaxSize=200 and Oligo=1")

    # convert fur output into primer3 compatible output using the functions from the fur2primer3 script
    print ("Convert the FUR output to Primer3 compatible output using fur2primer3 functions...")
    try:
        # parse parameters (dictionary)
        params=args_to_dict(args.parameter)
        #remap/ convert parameter keys to Primer3 conventions
        form_param=remap_keys(params)
        #write Primer3 compatible file.
        write_result(Path(target),form_param)
    except FileNotFoundError:
        quit(f"Could not find {target}")
        raise
    except ValueError as e:
        quit(f"Input or output values are not as expected: {e}")
    except TypeError as e:
        quit(f"Wrong type: {e}")
    except OSError as e:
        quit(f"Could not write or open file: {e}")
    except Exception as e:
        quit(f"Unknown exception occured while running fur2primer3 functions: {e}")
    
    # define results file
    resultf2p = target.with_suffix('.primers.txt')

    # check if the results file is empty, if so sys.exit 1 with message (might not raise correct exception)
    if resultf2p.stat().st_size == 0:
        quit(f"{resultf2p} is empty.")

    # try running primer3_core, write it to file, check if the file exists and/or is empty. If so,
    try:
        print("Running primer3_core.")
        primer3 = subprocess.run(['primer3_core', str(resultf2p)], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        quit(f"primer3_core failed: {e}")

    resultp3 = resultf2p.with_suffix('.primer3_out.txt')
    resultp3.write_text(primer3.stdout.strip(),encoding="utf-8")

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
        with open(temp_file, 'w', encoding="utf-8") as tf:
            if is_primer:
                for element in value:
                    split_lines = re.split(r'\s*=\s*', element.strip())
                    if len(split_lines) >= 2:
                        tf.write(f">{split_lines[0]}\n{split_lines[1]}\n")
            elif "Data" in key:
                tf.write("\n".join(value))
            else:
                split_line = re.split(r'\s*=\s*', value.strip())
                if len(split_line) >= 2:
                    tf.write(f">{split_line[0]}\n{split_line[1]}\n")

    destination_folder_pr = source_folder / "FUR.P3.PRIMERS"
    destination_folder_tar = source_folder / "FUR.P3.TARGETS"
    destination_folder_data=destination_folder_pr/"primer_data"
    destination_folder_pr.mkdir(parents=True, exist_ok=True)
    destination_folder_tar.mkdir(parents=True, exist_ok=True)
    destination_folder_data.mkdir(parents=True, exist_ok=True)

    move_files(source_folder, destination_folder_tar, "*Target*.txt")
    move_files(source_folder, destination_folder_pr, "*Primer*.txt")
    move_files(source_folder, destination_folder_data, "*Data*.txt")

    print('Primer3_module.py ran to completion: exit status 0')

if __name__ == '__main__':
    main()
