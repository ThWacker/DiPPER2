#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
import re
import threading
import time
from multiprocessing import Pool, set_start_method
from typing import Optional
import psutil
from Bio import SeqIO  # type: ignore
from logging_handler import Logger

# ────────────────────────────────────────────────────────────────
# Memory monitoring
# ────────────────────────────────────────────────────────────────
def monitor_memory(pid, threshold, logger):
    """ 
    Memory monitoring to ensure that process does not exceed the limits.
    Will kill process if it exceeds limits.
    Run as a deamon, so that it is terminated when all non-daemon threads finish.
    Args:
        pid (int): the process ID
        threshold_mb (int): maximum memory in bytes
        logger (Logger): the logger
    Returns:
        None
    Raises:
        None
    """
    try:
        #use psutil. Is a C extension, will not be under GIL
        proc = psutil.Process(pid)
        while proc.is_running():
            #get resident set size (rss) in bytes
            mem = proc.memory_info().rss 
            if mem > threshold:
                logger.warning(f"[KILL] PID {pid} exceeded {mem / (1024 ** 2):.2f} MB")
                proc.kill()
                return
            time.sleep(1)
    # process is already killed
    except psutil.NoSuchProcess:
        pass

def wait_for_memory(min_available, logger, CHECK_INTERVAL):
    """
    Will wait till enough memory in bytes is available. Times out after 6h and kills the entire DiPPER2 run.
    Args:
        min_available (int): the minimum memory that needs to be available to exit the function and continue with the script
        logger (Logger): the logger
    Returns:
        None
    Raises:
        TimeoutError
    """
    #hardcoded timeout, needs to be changed to dynamically changeable.
    timeout=21600
    #get time in seconds
    start_time = time.time()

    while True:
        #get available memory
        available = psutil.virtual_memory().available

        # there is enough available. Yay. Let's proceed with the script
        if available >= min_available:
            return
        
        #meh, not enough memory. When if clause becomes True, raise TimeoutError exception
        if time.time() - start_time > timeout:
            logger.warning(f"[TIMEOUT] Memory check timed out after {timeout} seconds.")
            raise TimeoutError(f"Memory did not reach {min_available} bytes within timeout")
        
        #there is not enough memory but also we have not run out of time. Waiting...
        logger.info(f"[WAIT] Available memory {available}  < threshold. Waiting for {CHECK_INTERVAL} seconds.")
        time.sleep(CHECK_INTERVAL)

# ────────────────────────────────────────────────────────────────
# Main processing function
# ────────────────────────────────────────────────────────────────
def process_file_amplicon(

    file_path: Path,
    frwd: str,
    rev: str,
    mismatch: int,
    logger,
    MAX_MEMORY_MB_PER_JOB:int,
    timeout: Optional[int]= None
    ):
    """
    This is where the actual magic happens. At least part of it. Seqkit amplicon is spawned as a separate process.
    A memory monitoring thread is started as a daemon. While the current process waits for seqkit_out.communicate, 
    this daemon monitors memory usage and kills the process if it starts to eat too much memory.

    Args:
        file_path (Path): the assembly file
        frwd (str): the forward primer
        rev (str): the reverse primer
        mismatch (int): mismatches allowed in seqkit
        logger (Logger): the logger
        timeout (int): how many seconds seqkit amplicon can run before it times out

    Returns:
       output (str): text of output of seqkit amplicon (in bed format)

    Raises:
        subprocess.TimeoutExpired
        Exception
    """
    try:
        with open(file_path, 'r') as f:
            file_contents = f.read()
    except OSError as e:
        logger.exception("could not open or read file %s, error %s", file_path, e)
        raise OSError(f"could not open or read file {file_path} with {e}") from e
    try:
        seqkit_out = subprocess.Popen(
            ["seqkit", "amplicon", "-F", frwd, "-R", rev, "--bed", "-m", str(mismatch)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Initializing of memory monitoring daemon
        memory_thread = threading.Thread(target=monitor_memory, args=(seqkit_out.pid, MAX_MEMORY_MB_PER_JOB, logger))
        memory_thread.daemon = True
        memory_thread.start()
        
        # is there a timeout? get the results from seqkit
        if timeout is not None:
            print(f"{timeout} for process {seqkit_out.pid} for file {file_path}. the time is {time.time()}")
            output, error = seqkit_out.communicate(input=file_contents, timeout=timeout)
            seqkit_out.wait() 
            print(f"the output of {seqkit_out.pid} is {len(output)} characters.The time is {time.time()}")
        else:
            output, error = seqkit_out.communicate(input=file_contents)
            seqkit_out.wait() 
        

        # well, this could have been put into an exception but here we are
        if seqkit_out.returncode != 0:
            logger.warning(f"[{file_path.name}] Seqkit failed: {error.strip()}")

        # Eureka
        logger.info(f"[{file_path.name}] Success with mismatch {mismatch}, length: {len(output)}")
        return output

    except subprocess.TimeoutExpired:
        logger.warning(f"[{file_path.name}] Timeout after {timeout}s")
        seqkit_out.kill()
        return None

    except Exception as e:
        logger.exception(f"[{file_path.name}] Unexpected error: {e}")
        return None

# ────────────────────────────────────────────────────────────────
# Worker wrapper
# ────────────────────────────────────────────────────────────────
def process_file_wrapper(args):
    """
    Pool.map manages that each process gets exactly one item (one tuple) from the args list at a time. 
    Each tuple corresponds to a single file f and its associated arguments. 
    Our wrapper passes the arguments of the tuple to the function. 
    Thus the function does not receive a tuple (one var), but each element. 
    """
    return process_file_amplicon(*args)

# ────────────────────────────────────────────────────────────────
# Run multiprocessing pool
# ────────────────────────────────────────────────────────────────
def run_seqkit_amplicon_with_optional_timeout(
    frwd: str,
    rev: str,
    concat: Path,
    mismatch: int,
    logger,
    CHECK_INTERVAL: int, 
    MIN_AVAILABLE_MEMORY: int,
    MAX_MEMORY_MB_PER_JOB: int,
    timeout: Optional[int]=None,
    max_workers: int = 6
    
):
    """
    Parallelizing seqkit using the multiprocessing package and pool.map.
    Args:
        frwd (str): forward primer
        rev (str): reverse primer
        concat (Path): the path to the folder with the assemblies, either FUR.target or FUR.neighbour
        mismatch (int): allowed mismatch for seqkit amplicon
        timeout (int): the time in seconds that can pass before seqkit amplicon runs time out
        logger (Logger): the logger
        max_workers (int): max worker processes spawned by pool.map

    Returns:
        successful (list): list of results
    
    Raises:
        TimeoutError
    """
    files = [f for f in concat.glob("*") if f.is_file()]
    if not files:
        logger.error(f"No files found in {concat}")
        raise FileNotFoundError

    logger.info(f"Processing {len(files)} files with mismatch={mismatch}...")
    
    # Check if enough memory is available
    try:
        wait_for_memory(MIN_AVAILABLE_MEMORY, logger, CHECK_INTERVAL)
    except TimeoutError as e:
        logger.error("[TIMEOUT] The memory timed out when trying to run seqkit amplicon with error %s", str(e))
        sys.exit(1)

    args = [(f, frwd, rev, mismatch, logger, MAX_MEMORY_MB_PER_JOB, timeout) for f in files]

    with Pool(processes=max_workers) as pool:
        results = pool.map(process_file_wrapper, args)

    successful = [res for res in results if res is not None]
    logger.info(f"Finished mismatch {mismatch}: {len(successful)} successful")
    return successful

# ────────────────────────────────────────────────────────────────
# Utility functions 
# ────────────────────────────────────────────────────────────────
def check_folders(*folders: Path,logger: Logger):
    """
    Check if folders exist and are non-empty.

    Args:
        folders (Path): path object(s) of one or more folders

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
        programm (str): the program to check

    Returns:
        None
    """
    # which is a tool in shutil.
    return shutil.which(program) is not None

def move_files_with_pattern(source_dir: Path, pattern: str, destination_dir: Path):
    """
    Move files containing a specific pattern from source to destination.

    Args:
        source_dir (Path): path object of folder containing files to be moved
        pattern (str): the pattern that needs to be in the file name for the file to be shifted into the destination folder
        destination_dir (Path): path object of the folder the files with a matching pattern in the file are moved to

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

# ────────────────────────────────────────────────────────────────
# Primer testing related functions
# ────────────────────────────────────────────────────────────────
def extract_primer_sequences(file: Path, logger: Logger) -> tuple[str, str, str]:
    """
    Extract primer sequences from a file.

    Args:
        file (Path): a path object of the file containing primers (generated by Primer3 module)

    Returns:
        a tuple of the primers and the internal probe.
    """

    # initialize empty dictionary
    sequences = {"PRIMER_LEFT": None, "PRIMER_RIGHT": None, "PRIMER_INTERNAL": None}

    # if you find a line with the keys from sequences, then add the next line as a value to the key in this dictionary
    with file.open("r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            for key in sequences.keys():
                if key in line:
                    sequences[key] = lines[i + 1].strip()

    # did not find all keys (not all values in dictionary are truthy)? Throw error!
    if not all(sequences.values()):
        logger.error(
            "The primer headers are not correctly formatted and cannot be processed. Please change the headers accordingly."
        )
        sys.exit()
    return (
        sequences["PRIMER_LEFT"],
        sequences["PRIMER_RIGHT"],
        sequences["PRIMER_INTERNAL"],
    )

def run_seqkit_locate(amplicon: str, ref_file: Path,logger: Logger):
    """Runs seqkit locate on a reference or an assembly with the amplicon, finds its coordinates and returns them in a bed file

    Args:
        amplicon (str): the amplicon
        ref_file (Path): path object of file used as reference

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

        # check for errors THIS IS REDUNDANT, REMOVE IN NEXT ITERATION OF IMPROVAL
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
        file (Path): the bed-file with the amplicon in column 7

    Returns:
        The amplicon as a string.
    """
    # open file, split first line by tab delimiter, return the amplicon in column 7
    with file.open("r") as f:
        first_line = f.readline()
        amplicon = first_line.strip().split("\t")[6]
        return amplicon


def get_longest_target(directory: Path) -> str:
    """
    Within the target folder find the longest fasta

    Args:
        directory (Path): within the directory with the target assemblies, find the longest assembly.

    Returns:
        The string representation of the path to the file with the longest assembly
    """
    longest_length = 0
    longest_file = None

    # don't read in anything but fasta files
    for filename in os.listdir(directory):
        if (
            filename.endswith(".fasta")
            or filename.endswith(".fa")
            or filename.endswith(".fna")
        ):
            filepath = os.path.join(directory, filename)
            # Read sequences in the file
            sequences = SeqIO.parse(filepath, "fasta")

            # total length of all contigs together
            total_length = sum(len(seq) for seq in sequences)

            # if longer than previous longest assembly, replace with current assembly
            if total_length > longest_length:
                longest_length = total_length
                longest_file = filepath
        
    # make sure that this function returns something or fails gracefully    
    if longest_file is None:
        raise RuntimeError("No valid FASTA files found in the directory.")
    return longest_file

# ────────────────────────────────────────────────────────────────
# main
# ────────────────────────────────────────────────────────────────

def main():
    """
    One function to rule them all.
    """
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")

    parser = argparse.ArgumentParser(
        prog="DiPPER2",
        description="primer testing module. Primers generated by Primer3 are validated in silico.",
        epilog="complaints, bugs, coffee, puppies and suggestions to t.wacker2@exeter.ac.uk",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.1.0")
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
        "-m",
        "--max_mem",
        type=int,
        default=16000000000,
        help="Max memory in bytes a job can use before getting killed. Default:  16000000000."
    )
    parser.add_argument(
       "-c",
       "--check_interval",
       type=int,
       default=2,
       help="If sufficient memory for starting a job is not available, wait <-c> seconds and check again. Default: 2 " 
    )
    parser.add_argument(
        "-d",
        "--min_mem",
        type=int,
        default=5000000000,
        help="Min memory that needs to be available for a job to get started. Integer in bytes.",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int, 
        default=6,
        help="Number of workers in a pool of worker processes to which jobs can be submitted. Default: 6"
    )
    parser.add_argument(
        "-r", "--ref", type=str, help="Reference assembly of the targets."
    )
    parser.add_argument(
        "-V", "--verbose", action="store_true", help="increase logging verbosity"
    )
    args = parser.parse_args()
    source_folder = Path(args.folder)
    
    # configure multiprocessing: spawn and not fork
    set_start_method("spawn")

    # configures the logger
    module_name = Path(__file__).name
    logger_instance = Logger(module_name, source_folder, args.verbose)
    logger = logger_instance.get_logger()

    # Check if programs are installed
    check_program_installed("seqkit")
    check_program_installed("blastx")

    # get args for constants
    MAX_WORKERS = args.workers
    MAX_MEMORY_MB_PER_JOB =  args.max_mem
    MIN_AVAILABLE_MEMORY = args.min_mem
    CHECK_INTERVAL = args.check_interval

    # Define the folders with primers and targets
    destination_folder_pr = source_folder / "FUR.P3.PRIMERS"
    destination_folder_tar = source_folder / "FUR.P3.TARGETS"

    # Check if they exist
    check_folders(destination_folder_pr, destination_folder_tar, logger=logger)

    # Define folders
    fur_target = source_folder / "FUR.target"
    fur_neighbour = source_folder / "FUR.neighbour"

    # Check if they exist and are not empty
    check_folders(fur_target, fur_neighbour, logger=logger)

    all_files = list(destination_folder_pr.glob("*"))

    for file_path in all_files:
        if file_path.is_file():
            logger.info("Testing your primers in: %s ", file_path)
            try:
                pr_frwd, pr_rev, pr_intern = extract_primer_sequences(file_path, logger)
            except Exception as e:
                logger.exception(
                    "Could not extract primer sequences from %s: %s",file_path, e
                )
                raise RuntimeError(
                    f"Could not extract primer sequences from {file_path}: {e}"
                ) from e

            # runs seqkit amplicon for targets with max mismatches of 4
            for i in range(4):
                try:
                    out_seqk_target = run_seqkit_amplicon_with_optional_timeout(
                        frwd=pr_frwd,
                        rev=pr_rev,
                        concat=fur_target,
                        mismatch=i,
                        logger=logger,
                        CHECK_INTERVAL=CHECK_INTERVAL,
                        MIN_AVAILABLE_MEMORY=MIN_AVAILABLE_MEMORY,
                        MAX_MEMORY_MB_PER_JOB=MAX_MEMORY_MB_PER_JOB,
                        timeout=None,
                        max_workers=MAX_WORKERS
                    )
                except subprocess.CalledProcessError as e:
                    logger.exception("Error running seqkit amplicon: %s", e)
                    raise subprocess.CalledProcessError(
                        returncode=-1, cmd="seqkit amplicon", output="", stderr=str(e)
                    ) from e
                except OSError as e:
                    logger.exception(
                        "Error with the operating system while running seqkit amplicon:%s OR no files found in %s", e, fur_target
                    )
                    raise OSError(
                        "Error with the operating system while running seqkit amplicon:%s OR no files found in %s", e, fur_target
                    ) from e
                except Exception as e:
                    logger.exception(
                        "Unknown exception/ unexpected error running seqkit amplicon: %s", e
                    )
                    raise RuntimeError(
                       f"Unknown exception/ unexpected error running seqkit amplicon: {e}"
                    ) from e

                if not out_seqk_target:
                    logger.warning(
                        "Seqkit amplicon did not return any matches for the primers in the targets with -m flag at %d", i
                    )
                    continue

                try:
                    # Construct filename for output
                    filename = f"{file_path}_seqkit_amplicon_against_target_m{i}.txt"
                    # Write output to the file
                    with open(filename, "w", encoding="utf-8") as file:
                        for line in out_seqk_target:
                            file.write(line)
                except (IOError, OSError, PermissionError) as e:
                    # Log error and raise exception with additional context
                    logger.exception(
                        "Error writing output of seqkit amplicon to file %s: %s", filename, e
                    )
                    raise RuntimeError(
                        f"Error writing output of seqkit amplicon to file {filename}: {e}"
                    ) from e
            # run seqkit amplicon for neighbours with up to 5 mismatches. Time out after 30s
            for i in range(5):
                try:
                    out_seqk_neighbour = run_seqkit_amplicon_with_optional_timeout(
                        frwd=pr_frwd,
                        rev=pr_rev,
                        concat=fur_neighbour,
                        mismatch=i,
                        logger=logger,
                        CHECK_INTERVAL=CHECK_INTERVAL,
                        MIN_AVAILABLE_MEMORY=MIN_AVAILABLE_MEMORY,
                        MAX_MEMORY_MB_PER_JOB=MAX_MEMORY_MB_PER_JOB,
                        timeout=30,
                        max_workers=MAX_WORKERS,
                    )
                except Exception as e:
                    logger.exception("Error running seqkit amplicon: %s", e)
                    raise Exception(f"Error running seqkit amplicon: {e}") from e

                if not out_seqk_neighbour:
                    logger.warning(
                        "Seqkit amplicon did not return any matches for the primers in the neighbours with -m flag at %d", i
                    )
                    continue

                try:
                    # Construct filename for output
                    filename = f"{file_path}_seqkit_amplicon_against_neighbour_m{i}.txt"

                    # Write output to the file
                    with open(filename, "w", encoding="utf-8") as file:
                        for line in out_seqk_neighbour:
                            file.write(line)
                except (IOError, OSError, PermissionError) as e:
                    # Log error and raise exception with additional context
                    logger.error(
                        "Error writing output of seqkit amplicon to file %s: %s", filename, e
                    )
                    raise RuntimeError(
                        f"Error writing output of seqkit amplicon to file {filename}: {e}"
                    ) from e

   # Log an informational message to indicate that the blastx command is starting
    logger.info("Running blastx on the targets...")

    # Get a list of all files in the destination folder
    all_files_tar = list(destination_folder_tar.glob("*"))

    # Iterate through each file in the list
    for file_path_tar in all_files_tar:
        if file_path_tar.is_file():  # Check if the current item is a file, not a directory
            print(f"{file_path_tar}")  # Print the file path for tracking purposes
            output_tar = "" #initialize
            try:
                # Run the blastx command with the necessary parameters and capture stdout and stderr
                result = subprocess.run(
                    [
                        "blastx",  # The blastx command for sequence alignment
                        "-query", str(file_path_tar),  # Input file (query)
                        "-remote",  # Use remote database (instead of local)
                        "-db", "nr",  # Database to query against (nr - non-redundant)
                        "-evalue", "0.00001",  # E-value threshold for the alignment
                        "-outfmt", "6",  # Output format (tabular)
                    ],
                    stdout=subprocess.PIPE,  # Capture standard output
                    stderr=subprocess.PIPE,  # Capture standard error
                    text=True,  # Ensure output is returned as text
                    check=True,  # Raise an error if the subprocess fails
                )
                output_tar = result.stdout  # Store the output of the blastx command
            except Exception as e:
                logger.exception("Blastx failed: %s", e)
                logger.warning("Proceeding with seqkit locate instead.")

            # If no output is generated by blastx, log a message and proceed with further steps
            if not output_tar:
                logger.info("Blastx did not return any results. No matches found.")

                # Get the amplicon related to the current target file
                file = (
                    destination_folder_pr
                    / f"{file_path_tar.name}_seqkit_amplicon_against_target_m0.txt"
                )
                file = Path(str(file).replace("Target", "Primer"))  # Adjust the file path for primer
                amp = get_amplicon(file)  # Get the amplicon from the file
                logger.info("The amplicon is %s", amp)

                # Check if a reference is provided, otherwise use the longest target assembly
                try:
                    if args.ref:
                        ref = Path(args.ref)  # Use provided reference
                        seqk_loc_out = run_seqkit_locate(amp, ref, logger)  # Run seqkit locate with the reference
                    else:
                        logger.warning("No reference found, using longest target assembly")
                        ref = get_longest_target(fur_target)  # Get the longest target assembly
                        if not ref:  # If no valid assembly is found, log and continue
                            logger.warning(
                                "No valid assembly found in %s to run seqkit locate. Do the assembly fasta files end on .fa, .fasta, or .fna?", fur_target
                            )
                            continue
                        logger.info("Using longest target assembly: %s", ref)
                        seqk_loc_out = run_seqkit_locate(amp, ref, logger)  # Run seqkit locate with the longest assembly

                except Exception as e:
                    # If an error occurs during seqkit locate, log it and raise an exception
                    logger.error("Error running seqkit locate: %s", e)
                    raise Exception(f"Unknown exception running seqkit locate: {e}") from e

                # If no results are returned from seqkit locate, log a warning and continue
                if not seqk_loc_out:
                    logger.warning(
                        f'Seqkit locate did not return a bed file for the assembly {args.ref if args.ref else ref} with the amplicon "{amp}".\n'
                    )
                    continue
            
                try:
                    # Generate a file name for the bed file and write seqkit locate results to it
                    match_no = re.search(r"_(\d+)\.txt", file_path_tar.name)
                    ref = Path(ref)
                    filename = f"Primer_{int(match_no.group(1))}_amplicon_locate_in_{ref.name}.bed"
                    filename = source_folder / filename
                    logger.info("Printing bed file for seqkit locate to %s", filename)
                    with open(filename, "w", encoding="utf-8") as file:
                        logger.info("Writing results of seqkit locate to bed file...")
                        file.write(seqk_loc_out)  # Write the locate results to the bed file
                except OSError as e:
                    # If an error occurs while writing the output to the file, log it and raise an exception
                    logger.error("Error writing output of seqkit locate to file %s: %s", filename, e)
                    raise OSError(f"Error writing output of seqkit locate to file {filename}: {e}") from e
            else:
                try:
                    # Write the blastx output to a text file
                    filenamed = f"{file_path_tar}_blastx_1e-5.txt"
                    with open(filenamed, "w", encoding="utf-8") as file_1:
                        file_1.write(output_tar)  # Save blastx results to a file
                except Exception as e:
                    # If an error occurs while writing the blastx output, log it and raise an exception
                    logger.error("Error writing output of blastx to file %s: %s", filenamed, e)
                    raise Exception(f"Error writing output of blastx to file {filenamed}: {e}") from e

    # Move files that are related to seqkit testing into a subfolder called "in_silico_tests"
    in_silico_folder = destination_folder_pr / "in_silico_tests"
    in_silico_folder.mkdir(parents=True, exist_ok=True)  # Create the subfolder if it doesn't exist
    pattern_to_match = "seqkit_amplicon_against"  # Pattern to search for in the file names

    logger.info(f"Moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}...")

    try:
        # Move files matching the pattern from the destination folder to the in_silico_tests folder
        move_files_with_pattern(destination_folder_pr, pattern_to_match, in_silico_folder)
    except Exception as e:
        # If an error occurs during the file moving process, log it and raise an exception
        logger.error(f"Error moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}: {e}")
        raise Exception(f"Error moving files with {pattern_to_match} in name from {destination_folder_pr} into {in_silico_folder}: {e}") from e

    # Log that the script has completed successfully
    logger.info("Primer_Testing_module.py ran to completion: exit status 0")


    # Exit the script successfully
    sys.exit(0)

if __name__ == "__main__":
    main()
