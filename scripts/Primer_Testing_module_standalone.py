#!/usr/bin/python

import os
import sys
import argparse
import subprocess
import shutil
import re
from datetime import datetime
import glob

def exit_program(message=None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    sys.exit(1)

def check_seqkit_installed():
    """Check if seqkit is installed and available in the PATH."""
    try:
        subprocess.run(["seqkit"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        exit_program("seqkit is not installed or not found in the PATH. Please install seqkit and ensure it is in the PATH.")

def extract_primer_sequences(file: str) -> tuple[str, str, str]:
    """Extract primer sequences that are formatted in the primer3 output way"""
    forward = "PRIMER_LEFT"
    reverse = "PRIMER_RIGHT"
    internal = "PRIMER_INTERNAL"
    forward_sequence, reverse_sequence, internal_sequence = None, None, None
    
    with open(file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if forward in line:
                forward_sequence = lines[i+1].strip()
            elif reverse in line:
                reverse_sequence = lines[i+1].strip()
            elif internal in line:
                internal_sequence = lines[i+1].strip()
    
    if not forward_sequence or not reverse_sequence: # or not internal_sequence:
        exit_program("The primer headers are not correctly formatted and cannot be processed. Please change the headers accordingly.")
    
    return forward_sequence, reverse_sequence, internal_sequence

def extract_target_sequence(file: str) -> str:
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if not ">" in line:
                output = line.strip()
    return output

def run_seqkit_amplicon(frwd: str, rev: str, concat: str, number: int) -> str:
    try:
        print(f"Running seqkit amplicon against {concat} with -m {number}.")
        cat = subprocess.Popen(["cat", concat], stdout=subprocess.PIPE, text=True)
        seqkit_out = subprocess.Popen(
            ['seqkit', 'amplicon', '-F', frwd, '-R', rev, '--bed', '-m', str(number)],
            stdin=cat.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        output, error = seqkit_out.communicate()

        if seqkit_out.returncode != 0:
            print(f"Error output from seqkit: {error}")
            raise subprocess.CalledProcessError(seqkit_out.returncode, 'seqkit amplicon')

        return output
    except subprocess.CalledProcessError:
        exit_program("seqkit amplicon failed.")

def usage():
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Theresa Wacker T.Wacker2@exeter.ac.uk')
    print('')
    print('Module runs fur2prim & primer3_core')
    print('')
    print('Usage:')
    print('optional:')
    print('-v/--version for the version')
    print('-o/--outfile_prefix - for outfile prefix [default: date and time in %d-%m-%Y_%Hh%Mmin%Ss_%z format]')
    print('-d/--delete_concat - if included, the concatenated fastas will be deleted.')
    print('')
    print('CAUTION: This module requires seqkit, blast+ and RRW-Primer-Blast (and thus also snakemake, as well as all dependencies) to be in path.')
    print('For more info on how to install and run RRW-Primer-Blast, compare ')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')

def main():
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")

    parser = argparse.ArgumentParser(description='primer testing module')
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument('-o', '--outfile_prefix', default=dt_string, type=str, help='Outfile prefix. Default is date and time in %d-%m-%Y_%Hh%Mmin%Ss_%z format')
    parser.add_argument('-i', '--infile_prefix', required=True, type=str, help="The prefix found on the _primer.fasta, _target.fasta and the targets_concat.fasta and neighbours_concat.fasta files.")
    parser.add_argument('-d', '--delete_concat', action='store_true', help='If option is used, the concatenated fastas will be deleted after module has finished.')

    args = parser.parse_args()

    check_seqkit_installed()
    
    # Define the directory containing the files
    source_folder = os.getcwd()
    prefix = args.infile_prefix

    # File names
    primer_f = f"{prefix}_primer.fasta"
    #target_f = f"{prefix}_target.fasta"
    target_concat = f"{prefix}_targets_concat.fasta"
    neighbour_concat = f"{prefix}_neighbours_concat.fasta"
    files = [primer_f, target_concat, neighbour_concat]#target_f

    # Check if files exist and are not empty
    for file_name in files:
        file_path = os.path.join(source_folder, file_name)
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            exit_program(f"{file_name} does not exist or is empty.")
        if not os.path.isfile(file_path):
            exit_program(f"{file_name} is not a file.")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # PRIMERS - SEQKIT AMPLICON #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    try:
        pr_frwd, pr_rev, pr_intern = extract_primer_sequences(primer_f)
    except Exception as e:
        exit_program(f"Could not extract primer sequences from {primer_f}: {e}")

    for i in range(3):  # Allow up to 2 mismatches and test each case
        try: 
            out_seqk_target = run_seqkit_amplicon(pr_frwd, pr_rev, target_concat, i)
        except Exception as e:
            exit_program(f"Error running seqkit amplicon: {e}")

        # Test that output is not empty, if empty, go to next iteration
        if not out_seqk_target:
            print(f"Seqkit amplicon did not return any matches for the primers in the targets with -m flag at {i}")
            continue

        try:
            filename = f"{prefix}_seqkit_amplicon_against_target_m{i}.txt"
            with open(filename, "w") as file:
                file.write(out_seqk_target)
        except Exception as e:
            exit_program(f"Error writing output of seqkit amplicon to file {filename}: {e}")

    for i in range(5):  # For neighbours, allow up to 4 mismatches and test each case
        try: 
            out_seqk_neighbour = run_seqkit_amplicon(pr_frwd, pr_rev, neighbour_concat, i)
        except Exception as e:
            exit_program(f"Error running seqkit amplicon: {e}")

        # Test that output is not empty, if empty, go to next iteration
        if not out_seqk_neighbour:
            print(f"Seqkit amplicon did not return any matches for the primers in the neighbours with -m flag at {i}")
            continue

        try:
            filename = f"{prefix}_seqkit_amplicon_against_neighbour_m{i}.txt"
            with open(filename, "w") as file:
                file.write(out_seqk_neighbour)
        except Exception as e:
            exit_program(f"Error writing output of seqkit amplicon to file {filename}: {e}")

    # #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # # TARGETS - BLASTX #
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    # try:
    #     target_blast = subprocess.run(
    #         ['blastx', '-query', target_f, '-remote', '-db', 'nr', '-evalue', '0.00001', '-outfmt', '7'],
    #         check=True,
    #         stdout=subprocess.PIPE,
    #         stderr=subprocess.PIPE,
    #         text=True
    #     )
    # except subprocess.CalledProcessError as e:
    #     exit_program(f"Blastx failed: {e}")

    # # Test that output is not empty, if empty, print a message
    # if not target_blast.stdout:
    #     print("Blastx did not return any results. No matches found.")


    print('Primer_Testing_module.py ran to completion: exit status 0')

    # Delete concatenated files if requested
    if args.delete_concat:
        for file_name in [target_concat, neighbour_concat]:
            try:
                os.remove(os.path.join(source_folder, file_name))
                print(f"Deleted {file_name}")
            except OSError as e:
                print(f"Error deleting {file_name}: {e}")

if __name__ == '__main__':
    main()
