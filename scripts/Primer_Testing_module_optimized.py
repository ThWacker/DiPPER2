#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
import re
from Bio import SeqIO  # type: ignore
from logging_handler import Logger


def check_folders(*folders: Path):
    """
    Check if folders exist and are non-empty.
    Args:
        folders(Path): path object(s) of one or more folders
    Returns:
        None
    """
    for folder in folders:
        if not folder.is_dir():
            logger.error(f"The folder {folder} does not exist")
            sys.exit(f"The folder {folder} does not exist")
        if not any(folder.iterdir()):
            logger.error(f"The folder {folder} is empty")
            sys.exit(f"The folder {folder} is empty")


def check_program_installed(program: str):
    """
    Check if a program is installed and available in the PATH.
    Args:
        programm(str): the program to check
    Returns:
        None
    """
    try:
        subprocess.run(
            program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False
        )
    except FileNotFoundError as e:
        logger.exception(
            f"{program} is not installed or not found in the PATH.Please install {program} and ensure it is in the PATH.",
            exc_info=1,
        )
        raise FileNotFoundError(
            f"{program} is not installed or not found in the PATH.Please install {program} and ensure it is in the PATH."
        ) from e


def extract_primer_sequences(file: Path) -> tuple[str, str, str]:
    """
    Extract primer sequences from a file.
    Args:
        file(Path): path object of file that contains primers
    Returns:
        A tuple of strings with the forward, reverse and internal probe.
    """
    forward, reverse, internal = "PRIMER_LEFT", "PRIMER_RIGHT", "PRIMER_INTERNAL"
    sequences = {forward: None, reverse: None, internal: None}

    with file.open("r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            for key in sequences.keys():
                if key in line:
                    sequences[key] = lines[i + 1].strip()

    if not all(sequences.values()):
        logger.error(
            "In function 'extract_primer_sequences': The primer headers are not correctly formatted and cannot be processed.Please change the headers accordingly."
        )
        raise Exception(
            "The primer headers are not correctly formatted and cannot be processed.Please change the headers accordingly."
        )
    return sequences[forward], sequences[reverse], sequences[internal]


def concat_files(folder: Path, name: str, source: Path) -> str:
    """
    Concatenate the content of all files found in a folder.
    Args:
        folder(Path): a path object of the folder which contains the fastas to be concatenated
        name(str): target or neighbour
        source(Path): Path object to the folder that rules it all (the results folder of this DiPPER2 run)
    Returns:
        the filename of the concatenated fasta file
    Raises:
        FileNotFoundError
    """
    outfilename = source / f"{name}_concatenated.fasta"
    with open(outfilename, "wb") as outfile:
        try:
            for filename in folder.glob("*"):
                if filename.name == outfilename:
                    continue
                with filename.open("rb") as readfile:
                    shutil.copyfileobj(readfile, outfile)
        except FileNotFoundError as e:
            logger.exception(
                f"Could not open or read files {folder}. Concatenation failed",
                exc_info=1,
            )
            raise FileNotFoundError(
                f"Could not open or read files {folder}. Concatenation failed: {e}"
            ) from e
    return outfilename


def run_seqkit_amplicon_with_optional_timeout(
    frwd: str, rev: str, concat: str, number: int, timeout: int = None
):
    """
    Run seqkit amplicon with an optional timeout and return the output.

    Args:
        frwd (str): The forward primer sequence.
        rev (str): The reverse primer sequence.
        concat (str): Path to the concatenated file.
        number (int): The number of allowed mismatches.
        timeout (Optional[int]): Timeout in seconds for the subprocess. If None, no timeout is applied.

    Returns:
        Optional[str]: Output of seqkit amplicon as a string if successful, None if it times out.

    Raises:
        ValueError: If invalid arguments are provided.
        subprocess.CalledProcessError: If seqkit amplicon fails with a non-zero exit code.
    """
    # Validate inputs
    if not all([frwd, rev, concat]):
        raise ValueError(
            "Forward primer, reverse primer, and concatenated file must be provided."
        )
    if number < 0:
        logger.error(
            "Error in function call of run_seqkit_amplicon_with_optional_timeout. Mismatch number must be a non-negative integer."
        )
        raise ValueError("Mismatch number must be a non-negative integer.")

    logger.info(
        f"Running seqkit amplicon on {concat} with primers {frwd} (forward) and {rev} (reverse)."
    )

    try:
        # Step 1: Create subprocess to read the input file using 'cat'
        cat = subprocess.Popen(
            ["cat", concat],
            stdout=subprocess.PIPE,  # Send output to the next process
            text=True,  # Enable text mode for I/O
        )
        logger.debug(f"'cat {concat}' subprocess started successfully.")

        # Step 2: Pipe the output of 'cat' into 'seqkit amplicon'
        logger.debug("'seqkit amplicon' subprocess started successfully.")
        seqkit_out = subprocess.Popen(
            ["seqkit", "amplicon", "-F", frwd, "-R", rev, "--bed", "-m", str(number)],
            stdin=cat.stdout,  # Input from 'cat' command
            stdout=subprocess.PIPE,  # Capture standard output
            stderr=subprocess.PIPE,  # Capture standard error
            text=True,  # Enable text mode for I/O
        )

        # Step 3: Handle timeout if provided
        if timeout is not None:
            output, error = seqkit_out.communicate(timeout=timeout)
        else:
            output, error = seqkit_out.communicate()

        # Step 4: Check for errors
        if seqkit_out.returncode != 0:
            logger.error(f"Seqkit error output: {error.strip()}")
            raise subprocess.CalledProcessError(
                seqkit_out.returncode, "seqkit amplicon", output=error.strip()
            )

        logger.info(
            f"Seqkit amplicon ran successfully. Output size: {len(output)} characters."
        )
        return output

    except subprocess.TimeoutExpired:
        # Handle timeout scenarios
        logger.warning(f"Seqkit amplicon timed out after {timeout} seconds.")
        seqkit_out.kill()  # Ensure seqkit process is terminated
        if cat:
            cat.kill()  # Ensure 'cat' process is terminated
        return None

    except subprocess.CalledProcessError as e:
        # Handle non-zero exit codes from seqkit
        logger.exception(
            f"Seqkit amplicon failed with exit code {e.returncode}: {e.output}"
        )
        raise subprocess.CalledProcessError(
            returncode=-1, cmd="seqkit amplicon", output="", stderr=str(e)
        ) from e

    except Exception as e:
        # Handle any unexpected exceptions
        logger.exception(f"An unexpected error occurred: {str(e)}")
        raise Exception(f"An unexpected error occurred: {str(e)}") from e

    finally:
        # Cleanup: Ensure all subprocesses are terminated
        if cat and cat.poll() is None:
            cat.terminate()
        if seqkit_out and seqkit_out.poll() is None:
            seqkit_out.terminate()


def move_files_with_pattern(source_dir: Path, pattern: str, destination_dir: Path):
    """
    Move files containing a specific pattern from source to destination.
    Args:
        source_dir(Path): path object of folder containing files to be moved
        pattern(str): the pattern that needs to be in the file name for the file to be shifted into the destination folder
        destination_dir(Path): path object of the folder the files with a matching pattern in the file are moved to
    Returns:
        None
    """
    # create destination folder unless it already exists, do not raise FileExistsError
    destination_dir.mkdir(exist_ok=True)

    # iterate through the directory, move files that have a certain pattern in their name
    for file in source_dir.iterdir():
        if pattern in file.name and file.is_file():
            destination_file = destination_dir / file.name
            shutil.move(str(file), str(destination_file))


def run_seqkit_locate(amplicon: str, ref_file: Path):
    """Runs seqkit locate on a reference or an assembly with the amplicon, finds its coordinates and returns them in a bed file
    Args:
        amplicon(str): the amplicon
        ref_file(Path): path object of file used as reference
    Returns:
        seqkit locate output in bed format
    Raises:
        subprocess.CalledProcessError
    """
    try:
        logger.info(
            f"Running seqkit locate on the following assembly {ref_file} with the amplicon."
        )
        # Create subprocess to read the input file using 'cat'
        cat = subprocess.Popen(["cat", ref_file], stdout=subprocess.PIPE, text=True)
        # Pipe that output into subprocess that runs seqkit locate
        logger.debug("started subprocess seqkit locate")
        seqkit_out = subprocess.Popen(
            ["seqkit", "locate", "-p", amplicon, "--bed", "-m", "2"],
            stdin=cat.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        # get both the output and potential errors
        output, error = seqkit_out.communicate()

        # check for errors
        if seqkit_out.returncode != 0:
            logger.error(f"Error output from seqkit: {error}")
            raise subprocess.CalledProcessError(seqkit_out.returncode, "seqkit locate")

        logger.debug("seqkit locate ran successfully, returning the output...")

        return output
    # if subprocess failed and exception was not handled elsewhere
    except subprocess.CalledProcessError as e:
        logger.exception(
            f"seqkit locate failed with errorcode {e.returncode}: {e.output}",
            exc_info=1,
        )
        raise subprocess.CalledProcessError(
            returncode=-1, cmd="seqkit locate", output="", stderr=str(e)
        ) from e


def get_amplicon(file: Path) -> str:
    """
    Extract the amplicon from the sequence file.
    Args:
        file(Path): the bed-file with the amplicon in column 7
    Returns:
        The amplicon as a string.
    """
    # open file, split first line by tab delimiter, return the amplicon in column 7
    with file.open("r") as f:
        first_line = f.readline()
        amplicon = first_line.strip().split("\t")[6]
        return amplicon


def get_longest_target(directory: Path) -> Path:
    """
    Within the target folder find the longest fasta
    Args:
        directory(Path): within the directory with the target assemblies, find the longest assembly.
    Returns:
        The path to the file with the longest assembly
    """
    longest_length = 0
    longest_file = None

    for filename in os.listdir(directory):
        if (
            filename.endswith(".fasta")
            or filename.endswith(".fa")
            or filename.endswith(".fna")
        ):
            filepath = os.path.join(directory, filename)
            # Read sequences in the file
            sequences = SeqIO.parse(filepath, "fasta")
            total_length = sum(len(seq) for seq in sequences)

            if total_length > longest_length:
                longest_length = total_length
                longest_file = filepath
    return longest_file


def delete_concats(target: Path, neighbour: Path):
    """Delete concatenated fasta files."""
    try:
        os.remove(target)
        os.remove(neighbour)
        print("Concatenated files deleted.")
    except OSError as e:
        print(f"Error deleting concatenated files: {e}")


def main():
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")

    parser = argparse.ArgumentParser(
        prog="DiPPER2",
        description="primer testing module. Primers generated by Primer3 are validated in silico.",
        epilog="complaints, bugs, coffee, puppies and suggestions to t.wacker2@exeter.ac.uk",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        required=True,
        help="Results folder name, which includes results folders from previous steps",
    )
    parser.add_argument(
        "-o",
        "--outfile_prefix",
        default=dt_string,
        type=str,
        help="Outfile prefix. Default is date and time in d-m-y-h-m-s-tz format",
    )
    parser.add_argument(
        "-d",
        "--delete_concat",
        type=int,
        default=1,
        help="If set to 0, the concatenated fastas will not be deleted after module has finished. Default: true.",
    )
    parser.add_argument(
        "-r", "--ref", type=str, help="Reference assembly of the targets."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="increase logging verbosity"
    )
    args = parser.parse_args()
    source_folder = Path(args.folder)
    
    # configures the logger
    module_name = Path(__file__).name
    logger_instance = Logger(module_name, source_folder, args.verbose)
    logger = logger_instance.get_logger()

    # Check if programs are installed
    check_program_installed("seqkit")
    check_program_installed("blastx")

    # Define the folders with primers and targets
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
    concat_t = concat_files(fur_target, "target", source_folder)
    concat_n = concat_files(fur_neighbour, "neighbour", source_folder)

    all_files = list(destination_folder_pr.glob("*"))

    for file_path in all_files:
        if file_path.is_file():
            logger.info(f"Testing your primers in {file_path}:\n")
            try:
                pr_frwd, pr_rev, pr_intern = extract_primer_sequences(file_path)
            except Exception as e:
                logger.exception(
                    f"Could not extract primer sequences from {file_path}: {e}",
                    exc_info=1,
                )
                raise Exception(
                    f"Could not extract primer sequences from {file_path}: {e}"
                ) from e

            # runs seqkit amplicon for targets with max mismatches of 4
            for i in range(4):
                try:
                    out_seqk_target = run_seqkit_amplicon_with_optional_timeout(
                        pr_frwd, pr_rev, concat_t, i
                    )
                except subprocess.CalledProcessError as e:
                    logger.exception(f"Error running seqkit amplicon: {e}")
                    raise subprocess.CalledProcessError(
                        returncode=-1, cmd="seqkit amplicon", output="", stderr=str(e)
                    ) from e
                except OSError as e:
                    logger.exception(
                        f"Error with the operating system while running seqkit amplicon: {e}"
                    )
                    raise OSError(
                        f"Error with the operating system while running seqkit amplicon: {e}"
                    ) from e
                except Exception as e:
                    logger.exception(
                        f"Unknown exception/ unexpected error running seqkit amplicon: {e}"
                    )
                    raise Exception(
                        f"Unknown exception/ unexpected error running seqkit amplicon: {e}"
                    ) from e

                if not out_seqk_target:
                    logger.warning(
                        f"Seqkit amplicon did not return any matches for the primers in the targets with -m flag at {i}"
                    )
                    continue

                try:
                    # Construct filename for output
                    filename = f"{file_path}_seqkit_amplicon_against_target_m{i}.txt"

                    # Write output to the file
                    with open(filename, "w", encoding="utf-8") as file:
                        file.write(out_seqk_target)
                except (IOError, OSError, PermissionError) as e:
                    # Log error and raise exception with additional context
                    logger.exception(
                        f"Error writing output of seqkit amplicon to file {filename}: {e}"
                    )
                    raise RuntimeError(
                        f"Error writing output of seqkit amplicon to file {filename}: {e}"
                    ) from e
            # run seqkit amplicon for neighbours with up to 5 mismatches. Time out after 8 min
            for i in range(5):
                try:
                    out_seqk_neighbour = run_seqkit_amplicon_with_optional_timeout(
                        pr_frwd, pr_rev, concat_n, i, timeout=480
                    )
                    logger.info(f"ran seqkit amplicon for {i} mismatches")
                except Exception as e:
                    logger.exception(f"Error running seqkit amplicon: {e}")
                    raise Exception(f"Error running seqkit amplicon: {e}") from e

                if not out_seqk_neighbour:
                    logger.warning(
                        f"Seqkit amplicon did not return any matches for the primers in the neighbours with -m flag at {i}"
                    )
                    continue

                try:
                    # Construct filename for output
                    filename = f"{file_path}_seqkit_amplicon_against_target_m{i}.txt"

                    # Write output to the file
                    with open(filename, "w", encoding="utf-8") as file:
                        file.write(out_seqk_target)
                except (IOError, OSError, PermissionError) as e:
                    # Log error and raise exception with additional context
                    logger.error(
                        f"Error writing output of seqkit amplicon to file {filename}: {e}"
                    )
                    raise RuntimeError(
                        f"Error writing output of seqkit amplicon to file {filename}: {e}"
                    )

    logger.info("Running blastx on the targets...")

    all_files_tar = list(destination_folder_tar.glob("*"))
    for file_path_tar in all_files_tar:
        if file_path_tar.is_file():
            print(f"{file_path_tar}")
            try:
                result = subprocess.run(
                    [
                        "blastx",
                        "-query",
                        str(file_path_tar),
                        "-remote",
                        "-db",
                        "nr",
                        "-evalue",
                        "0.00001",
                        "-outfmt",
                        "6",
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                output_tar = result.stdout
            except Exception as e:
                logger.exception(f"Blastx failed: {e}")
                raise Exception(f"Blastx failed: {e}") from e

            if not output_tar:
                logger.info("Blastx did not return any results. No matches found.")

                # Get the amplicon
                file = (
                    destination_folder_pr
                    / f"{file_path_tar.name}_seqkit_amplicon_against_target_m0.txt"
                )
                file = Path(str(file).replace("Target", "Primer"))
                amp = get_amplicon(file)
                logger.info(f"The amplicon is {amp}")

                # See if reference is present, otherwise use longest assembly
                try:
                    if args.ref:
                        ref = Path(args.ref)
                        seqk_loc_out = run_seqkit_locate(amp, ref)
                    else:
                        logger.warning(
                            "No reference found, using longest target assembly"
                        )
                        ref = get_longest_target(fur_target)
                        if not ref:
                            logger.warning(
                                f"No valid assembly found in {fur_target} to run seqkit locate. Do the assembly fasta files end on .fa, .fasta, or .fna?"
                            )
                            continue
                        logger.info(f"Using longest target assembly: {ref}")
                        seqk_loc_out = run_seqkit_locate(amp, ref)

                except Exception as e:
                    logger.error(f"Error running seqkit locate:{e}")
                    raise Exception(
                        f"Unknown exception running seqkit locate: {e}"
                    ) from e

                if not seqk_loc_out:
                    logger.warning(
                        f'Seqkit locate did not return a bed file for the assembly {args.ref if args.ref else ref} with the amplicon "{amp}".\n'
                    )
                    continue

                try:
                    match_no = re.search(r"_(\d+)\.txt", file_path_tar.name)
                    ref = Path(ref)
                    filename = f"Primer_{int(match_no.group(1))}_amplicon_locate_in_{ref.name}.bed"
                    filename = source_folder / filename
                    logger.info(f"Printing bed file for seqkit locate to {filename}")
                    with open(filename, "w", encoding="utf-8") as file:
                        logger.info("Writing results of seqkit locate to bed file...")
                        file.write(seqk_loc_out)
                except OSError as e:
                    logger.error(
                        f"Error writing output of seqkit locate to file {filename}: {e}"
                    )
                    raise OSError(
                        f"Error writing output of seqkit locate to file {filename}: {e}"
                    ) from e

            try:
                filenamed = f"{file_path_tar}_blastx_1e-5.txt"
                with open(filenamed, "w", encoding="utf-8") as file_1:
                    file_1.write(output_tar)
            except OSError as e:
                logger.error(f"Error writing output of blastx to file {filenamed}: {e}")
                raise OSError(
                    f"Error writing output of blastx to file {filenamed}: {e}"
                ) from e

    # Move everything in the FUR.P3.PRIMER folder that is a seqkit testing file in a subfolder called "in_silico_tests"
    in_silico_folder = destination_folder_pr / "in_silico_tests"
    in_silico_folder.mkdir(parents=True, exist_ok=True)
    pattern_to_match = "seqkit_amplicon_against"

    logger.info(
        f"Moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}..."
    )
    try:
        move_files_with_pattern(
            destination_folder_pr, pattern_to_match, in_silico_folder
        )
    except Exception as e:
        logger.error(
            f"Error moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}: {e}"
        )
        raise Exception(
            f"Error moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}: {e}"
        ) from e

    logger.info("Primer_Testing_module.py ran to completion: exit status 0")

    if args.delete_concat:
        delete_concats(concat_t, concat_n)

    sys.exit(0)


if __name__ == "__main__":
    main()
