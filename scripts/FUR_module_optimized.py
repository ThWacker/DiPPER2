#!/usr/bin/python

import sys
import argparse
import subprocess
from datetime import datetime
from pathlib import Path
import shutil

def quit(message=None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    sys.exit(1)

def usage():
    """ Usage of the script"""
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Theresa Wacker T.Wacker2@exeter.ac.uk')
    print('')
    print('Module runs FUR')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-f/--folder - results folder name, which includes results folders from previous steps')
    print('optional:')
    print('-v/--version for the version')
    print('-o/--outfile_prefix - for outfile prefix [default: date and time in %%d-%%m-%%Y_%%Hh%%Mmin%%Ss_%%z format]')
    print('-p/--parameter - tells fur how to run. -u -> only the first subtraction step; -U -> runs first Subtraction step and Intersection step. -m -> runs megablast instead of blastn in the last step. [default: " " (runs FUR to completion with all steps using blastn)]')
    print('')
    print('-r/--reference - refence assembly for the makeFurDb step')
    print('')
    print('CAUTION: needs to be in the folder with the assemblies to be sorted in target and neighbour')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')

def check_folders(target_folder, neighbour_folder):
    """Check if folders exist and are non-empty."""
    for folder in [target_folder, neighbour_folder]:
        if not folder.is_dir():
            quit(f"The folder {folder} does not exist")
        if not any(folder.iterdir()):
            quit(f"The folder {folder} is empty")

def make_fur_db(target_folder, neighbour_folder, source, reference=None):
    """Create the FUR database."""
    try:
        print("Making FUR DB.")
        fur_db = source / 'FUR.db'
        if reference:
            subprocess.run(['makeFurDb', '-t', str(target_folder), '-n', str(neighbour_folder), '-d', str(fur_db), '-r', str(reference)], check=True)
        else:
            subprocess.run(['makeFurDb', '-t', str(target_folder), '-n', str(neighbour_folder), '-d', str(fur_db)], check=True)
    except subprocess.CalledProcessError as e:
        quit("makeFurDb failed.")
        raise e

def run_fur(option, outfile_prefix, source):
    """Run the FUR command and handle the output."""
    option_flag = f"-{option.strip()}" if option.strip() else ""
    fur_out = source / f"{outfile_prefix}_FUR.db.out.txt"
    fur_db = source / 'FUR.db'
    try:
        print("Running FUR.")
        result = subprocess.run(['fur', '-d', str(fur_db), option_flag], check=True, capture_output=True, text=True)
        fur_output = result.stdout.strip()
        error=result.stderr.strip()
    except subprocess.CalledProcessError as e:
        quit(f"FUR failed: {e}")
        raise e

    fur_out.write_text(fur_output, encoding="utf-8" )

    if not fur_out.exists() or fur_out.stat().st_size == 0:
        quit(f"Fur could not find unique regions or did not run successfully. {fur_out} is empty or does not exist.")

    return error

def clean_up(db):
    """Clean up temporary files."""
    if db.exists():
        shutil.rmtree(db)
        print(f"The directory {db} has been deleted.", file=sys.stderr)

def main():
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")
    #Parse all arguments
    parser = argparse.ArgumentParser(description='Runs FUR')
    parser.add_argument('-p', '--parameter', type=str, default=" ", help='Tells FUR how to run. -u -> only the first subtraction step; -U -> runs first Subtraction step and Intersection step; -m -> runs megablast instead of blastn in the last step; without argument runs FUR in full.')
    parser.add_argument('-o', '--outfile_prefix', default=dt_string, type=str, help='Outfile prefix. Default is date and time in d-m-y-h-m-s-tz format')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which includes results folders from previous steps')
    parser.add_argument('-r', '--reference', type=str, default=None,help="Optional reference for the makeFurDb step" )

    args = parser.parse_args()
    
    print(f'The parameters are:\n for the FUR parameter {args.parameter}\n for the outfile prefix: {args.outfile_prefix}', file=sys.stderr)

    # find the source folder, define the FUR.target and FUR.neighbour folders
    source_folder = Path(args.folder)
    target_folder = source_folder / "FUR.target"
    neighbour_folder = source_folder / "FUR.neighbour"
    
    #check if they exist
    check_folders(target_folder, neighbour_folder)

    # run makefurdb
    if args.reference:
        try:
            reference= args.reference
            make_fur_db(target_folder, neighbour_folder, source_folder, reference)
        except (FileNotFoundError, OSError):
            print(f"Could not find file {reference}. Defaulting to running without a reference")
            make_fur_db(target_folder, neighbour_folder, source_folder)
    else:
        make_fur_db(target_folder, neighbour_folder, source_folder)

    #run fur
    fur_out = run_fur(args.parameter, args.outfile_prefix, source_folder)

    #print the stats to stderr & save in file(might not want that on server)
    print(f'{fur_out}')
    fur_file= source_folder/ "FUR_Summary_output.txt"
    try:
        with open(fur_file, 'w', encoding="utf-8") as f:
            f.write(fur_out)
    except OSError as e:
        raise OSError(f"Could not open or write the FUR summary output to {fur_file}") from e



    # clean up after yourself
    clean_up(source_folder / "FUR.db")

    #Heureka
    print('FUR_module.py ran to completion: exit status 0')

if __name__ == '__main__':
    main()
