#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
import re
from Bio import SeqIO # type: ignore
 
def check_folders(*folders):
    """Check if folders exist and are non-empty."""
    for folder in folders:
        if not folder.is_dir():
            sys.exit(f"The folder {folder} does not exist")
        if not any(folder.iterdir()):
            sys.exit(f"The folder {folder} is empty")

def exit_program(message=None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    usage()
    sys.exit(1)

def check_program_installed(program: str):
    """Check if a program is installed and available in the PATH."""
    try:
        subprocess.run(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    except FileNotFoundError:
        exit_program(f"{program} is not installed or not found in the PATH.Please install {program} and ensure it is in the PATH.")

def extract_primer_sequences(file: Path) -> tuple[str, str, str]:
    """Extract primer sequences from a file."""
    forward, reverse, internal = "PRIMER_LEFT", "PRIMER_RIGHT", "PRIMER_INTERNAL"
    sequences = {forward: None, reverse: None, internal: None}

    with file.open('r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            for key in sequences.keys():
                if key in line:
                    sequences[key] = lines[i + 1].strip()
  
    if not all(sequences.values()):
        exit_program("The primer headers are not correctly formatted and cannot be processed.Please change the headers accordingly.")
    return sequences[forward], sequences[reverse], sequences[internal]

def concat_files(folder: Path, name: str, source:Path) -> str:
    """Concatenate the content of all files found in a folder."""
    outfilename = source/ f"{name}_concatenated.fasta"
    with open(outfilename, 'wb') as outfile:
        for filename in folder.glob('*'):
            if filename.name == outfilename:
                continue
            with filename.open('rb') as readfile:
                shutil.copyfileobj(readfile, outfile)
    return outfilename

def run_seqkit_amplicon_with_optional_timeout(frwd: str, rev: str, concat: str, number: int, timeout: int = None) -> str:
    """Run seqkit amplicon with an optional timeout and return the output."""
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
        if timeout is not None:
            output, error = seqkit_out.communicate(timeout=timeout)
        else:
            output, error = seqkit_out.communicate()

        if seqkit_out.returncode != 0:
            print(f"Error output from seqkit: {error}")
            raise subprocess.CalledProcessError(seqkit_out.returncode, 'seqkit amplicon')

        return output
    except subprocess.TimeoutExpired:
        print(f"seqkit amplicon timed out after {timeout} seconds.")
        seqkit_out.kill()
        cat.kill()
        return None
    except subprocess.CalledProcessError:
        exit_program("seqkit amplicon failed.")

def move_files_with_pattern(source_dir: Path, pattern: str, destination_dir: Path):
    """Move files containing a specific pattern from source to destination."""
    destination_dir.mkdir(exist_ok=True)

    for file in source_dir.iterdir():
        if pattern in file.name and file.is_file():
            destination_file = destination_dir / file.name
            shutil.move(str(file), str(destination_file))

def run_seqkit_locate(amplicon: str, ref_file: Path):
    """Runs seqkit locate on a reference or an assembly with the amplicon, finds its coordinates and returns them in a bed file"""
    try:
        print(f"Running seqkit locate on the following assembly {ref_file} with the amplicon.")
        cat = subprocess.Popen(["cat", ref_file], stdout=subprocess.PIPE, text=True)
        seqkit_out = subprocess.Popen(
            ['seqkit', 'locate', '-p', amplicon,'--bed', '-m','2'],
            stdin=cat.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        output, error = seqkit_out.communicate()

        if seqkit_out.returncode != 0:
            print(f"Error output from seqkit: {error}")
            raise subprocess.CalledProcessError(seqkit_out.returncode, 'seqkit locate')
        
        return output
    
    except subprocess.CalledProcessError:
        exit_program("seqkit locate failed.")

def get_amplicon(file: Path) -> str:
    """Extract the amplicon from the sequence file."""
    with file.open('r') as f:
        first_line = f.readline()
        amplicon = first_line.strip().split('\t')[6]
        return amplicon
    
def get_longest_target(directory: Path) -> Path:
    """Within the target folder find the longest fasta """
    longest_length = 0
    longest_file = None

    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith (".fna"):
            filepath = os.path.join(directory, filename)
            # Read sequences in the file
            sequences = SeqIO.parse(filepath, "fasta")
            total_length = sum(len(seq) for seq in sequences)
            
            if total_length > longest_length:
                longest_length = total_length
                longest_file = filepath
    print(f"the output is a str? {isinstance(longest_file, str)}")
    return longest_file


def usage():
    """Print usage information."""
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Theresa Wacker T.Wacker2@exeter.ac.uk')
    print('')
    print('Module runs in silico PCR for primers and blastx for targets')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-f/--folder - results folder name, which includes results folders from previous steps')
    print('optional:')
    print('-v/--version for the version')
    print('-o/--outfile_prefix - for outfile prefix [default: date and time in %%d-%%m-%%Y_%%Hh%%Mmin%%Ss_%%z format]')
    print('-d/--delete_concat - if included, the concatenated fastas will be deleted.')
    print('')
    print('CAUTION: This module requires seqkit, blast+ and RRW-Primer-Blast (and thus also snakemake, as well as all dependencies) to be in path.')
    print('For more info on how to install and run RRW-Primer-Blast, compare ')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')

def delete_concats():
    """Delete concatenated fasta files."""
    try:
        os.remove("target_concatenated.fasta")
        os.remove("neighbour_concatenated.fasta")
        print("Concatenated files deleted.")
    except OSError as e:
        print(f"Error deleting concatenated files: {e}")

def main():
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")

    parser = argparse.ArgumentParser(description='primer testing module')
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which includes results folders from previous steps')
    parser.add_argument('-o', '--outfile_prefix', default=dt_string, type=str, help='Outfile prefix. Default is date and time in d-m-y-h-m-s-tz format')
    parser.add_argument('-d', '--delete_concat', type=int, default=1, help='If set to 0, the concatenated fastas will not be deleted after module has finished. Default: true.')
    parser.add_argument('-r','--ref', type=str,help='Reference assembly of the targets.' )

    args = parser.parse_args()

    # Check if programs are installed
    check_program_installed("seqkit")
    check_program_installed("blastx")
    
    # Define the folders with primers and targets
    source_folder = Path(args.folder)
    destination_folder_pr = source_folder / "FUR.P3.PRIMERS"
    destination_folder_tar = source_folder / "FUR.P3.TARGETS"

    # Check if they exist
    check_folders(destination_folder_pr, destination_folder_tar)
    
    # Define folders 
    fur_target = source_folder / "FUR.target"
    fur_neighbour = source_folder / "FUR.neighbour"
    
    # Check if they exist and are not empty
    check_folders(fur_target, fur_neighbour)
    
    # Concatenate all target and neighbour fastas
    concat_t = concat_files(fur_target, "target",source_folder)
    concat_n = concat_files(fur_neighbour, "neighbour", source_folder)

    all_files = list(destination_folder_pr.glob('*'))

    for file_path in all_files:
        if file_path.is_file():
            print(f"Testing your primers in {file_path}:\n")
            try:
                pr_frwd, pr_rev, pr_intern = extract_primer_sequences(file_path)
            except Exception as e:
                exit_program(f"Could not extract primer sequences from {file_path}: {e}")

            for i in range(4):
                try: 
                    out_seqk_target = run_seqkit_amplicon_with_optional_timeout(pr_frwd, pr_rev, concat_t, i)
                except subprocess.CalledProcessError as e:
                    exit_program(f"Error running seqkit amplicon: {e}")
                except OSError as e:
                    exit_program(f"Error with the operating system while running seqkit amplicon: {e}")
                except Exception as e:
                    exit_program(f"Unexpected error running seqkit amplicon: {e}")

                if not out_seqk_target:
                    print(f"Seqkit amplicon did not return any matches for the primers in the targets with -m flag at {i}")
                    continue

                try:
                    filename = f"{file_path}_seqkit_amplicon_against_target_m{i}.txt"
                    with open(filename, "w", encoding="utf-8") as file:
                        file.write(out_seqk_target)
                except (IOError, OSError) as e:
                    exit_program(f"Error writing output of seqkit amplicon to file {filename}: {e}")

            for i in range(5):
                try: 
                    out_seqk_neighbour = run_seqkit_amplicon_with_optional_timeout(pr_frwd, pr_rev, concat_n, i, timeout=480)
                    print(f"ran seqkit amplicon for {i} mismatches")
                except Exception as e:
                    exit_program(f"Error running seqkit amplicon: {e}")

                if not out_seqk_neighbour:
                    print(f"Seqkit amplicon did not return any matches for the primers in the neighbours with -m flag at {i}")
                    continue

                try:
                    filename = f"{file_path}_seqkit_amplicon_against_neighbour_m{i}.txt"
                    with open(filename, "w",encoding="utf-8") as file:
                        file.write(out_seqk_neighbour)
                except (IOError, OSError)  as e:
                    exit_program(f"Error writing output of seqkit amplicon to file {filename}: {e}")


    print("Running blastx on the targets...")

    all_files_tar = list(destination_folder_tar.glob('*'))
    for file_path_tar in all_files_tar:
        if file_path_tar.is_file():
            print (f"{file_path_tar}")
            try:
                result = subprocess.run(
                    ['blastx', '-query', str(file_path_tar), '-remote', '-db', 'nr', '-evalue', '0.00001', '-outfmt', '6'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
                )
                output_tar = result.stdout
            except Exception as e:
                exit_program(f"Blastx failed: {e}")

            if not output_tar:
                print("Blastx did not return any results. No matches found.")
           
                # Get the amplicon
                file = destination_folder_pr / f"{file_path_tar.name}_seqkit_amplicon_against_target_m0.txt"
                file = Path(str(file).replace("Target", "Primer"))
                amp = get_amplicon(file)
                print(f"The amplicon is {amp}")

                # See if reference is present, otherwise use longest assembly
                try:
                    if args.ref:
                        ref=Path(args.ref)
                        seqk_loc_out = run_seqkit_locate(amp, ref)
                    else:
                        print("No reference found, using longest target assembly")
                        ref = get_longest_target(fur_target)
                        if not ref:
                            print(f"No valid assembly found in {fur_target} to run seqkit locate. Do the assembly fasta files end on .fa, .fasta, or .fna?")
                            continue
                        print(f"Using longest target assembly: {ref}")
                        seqk_loc_out = run_seqkit_locate(amp, ref)

                except Exception as e:
                    exit_program(f"Error running seqkit locate: {e}")

                if not seqk_loc_out:
                    print(f"Seqkit locate did not return a bed file for the assembly {args.ref if args.ref else ref} with the amplicon \"{amp}\".\n")
                    continue

                try:
                    match_no = re.search(r'_(\d+)\.txt', file_path_tar.name)
                    ref=Path(ref)
                    filename = f"Primer_{int(match_no.group(1))}_amplicon_locate_in_{ref.name}.bed"
                    filename= source_folder / filename
                    print(f"Printing bed file for seqkit locate to {filename}")
                    with open(filename, "w", encoding="utf-8") as file:
                        print("Writing results of seqkit locate to bed file...")
                        file.write(seqk_loc_out)
                except OSError as e:
                    exit_program(f"Error writing output of seqkit locate to file {filename}: {e}")

                

            try:
                filenamed = f"{file_path_tar}_blastx_1e-5.txt"
                with open(filenamed, "w", encoding="utf-8") as file_1:
                    file_1.write(output_tar)
            except OSError as e:
                exit_program(f"Error writing output of blastx to file {filenamed}: {e}")

    # Move everything in the FUR.P3.PRIMER folder that is a seqkit testing file in a subfolder called "in_silico_tests"
    in_silico_folder = destination_folder_pr / "in_silico_tests"
    in_silico_folder.mkdir(parents=True, exist_ok=True)
    pattern_to_match = 'seqkit_amplicon_against'

    print(f"Moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}...")
    try:
        move_files_with_pattern(destination_folder_pr, pattern_to_match, in_silico_folder)
    except Exception as e:
        exit_program(f"Error moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}: {e}")


    print('Primer_Testing_module.py ran to completion: exit status 0')
    
    if args.delete_concat:
        delete_concats()

    sys.exit(0)

if __name__ == '__main__':
    main()
    
