#!/usr/bin/python

import sys
import argparse
from pathlib import Path
import subprocess
from shutil import which
import re
import os
from pprint import pformat
from pprint import pprint
import jinja2  # type: ignore
import pandas as pd  # type: ignore
from logging_handler import Logger


def generate_html_jinja(
    header: str,
    primer_frwd: str,
    primer_rev: str,
    primer_internal: str,
    sensi_pass: str,
    sensi_ass_no: int,
    sensi_m: str,
    speci_pass: str,
    speci_ass_no: int,
    speci_m: str,
    sensi: dict,
    speci: dict,
    res_target_str: str,
    primer_fold: Path,
    target_fold: Path,
    seqkit_fold: Path,
    source_folder: Path,
    qPCR: str,
    ampli_len: int,
):
    """
    Generate HTML content for the results.

    Args:
        All args are generated by the generate_results function
        header (str): the header with type of primer and file name of primer
        primer_frwd (str): the forward primer
        primer_rev (str): the reverse primer
        primer_internal (str): the internal probe (qPCR only, otherwise set to NA)
        sensi_pass (str): did it pass the sensitivity tests? (either "passed" or "failed")
        speci_pass (str): did it pass the specificity tests? (either "passed" or "failed")
        speci_ass_no (int): number of assemblies tested for specificity
        speci_m (str): the number of mismatches allowed in the test for seqkit
        sensi (dict): the dictionary with all sensitivity testing results
        speci (dict): the dictionary with all the specificity testing results
        res_target_str (str): results of the blastx testing
        primer_fold (Path): the path object of the folder the primers are in
        target_fold (Path): the Path object of the folder the targets are in
        seqkit_fold (Path): the Path object of the folder with the in silico testing results
        source_fold (Path): the Path object of the folder the entire thing is in (parent folder of all DiPPER2 results)
        qPCR (str): toggle of y or n to change formatting according to whether it is qPCR or not

    Returns:
        None
        
    """
    # hardcoded output file name
    HTML_OUTFILE = source_folder / "Results.html"

    # return the directory of the script the function is called from, making sure that relative paths are converted to absolut paths
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # create jinja2 environments, load templates
    environment = jinja2.Environment(
        loader=jinja2.FileSystemLoader(os.path.join(script_dir, "templates/"))
    )
    if qPCR == "y":
        template = environment.get_template("results_qPCR.html")
    else:
        template = environment.get_template("results.html")

    # render the html file with provided arguments
    content = template.render(
        header=header,
        primer_frwd=primer_frwd,
        primer_rev=primer_rev,
        primer_internal=primer_internal,
        ampli_length=ampli_len,
        sensi=pformat(sensi, depth=4),
        speci=pformat(speci, depth=4),
        res_target_str=res_target_str,
        primer_fold=primer_fold,
        target_fold=target_fold,
        seqkit_fold=seqkit_fold,
        sensi_pass=sensi_pass,
        speci_pass=speci_pass,
        sensi_ass_no=sensi_ass_no,
        sensi_m=sensi_m,
        speci_m=speci_m,
        speci_ass_no=speci_ass_no,
        source_folder=source_folder,
    )

    # appending
    with open(HTML_OUTFILE, "a", encoding="utf-8") as stream:
        stream.write(content)


def is_tool(name: str) -> bool:
    """
    Check whether `name` is on PATH and marked as executable.
    
    Args:
        name (str): the name of the program
    
    Returns:
        True or False (boolean)
        
    """
    # which is a tool in shutil.
    return which(name) is not None


def count_files(folda: Path) -> int:
    """
    Counts the number of files in a folder.

    Args:
        folda (Path): a Path object that is the folder you want to count the files in

    Returns:
        the number of files as an integer
    """
    # kind of an elegant way to count files in folders
    return sum(1 for x in folda.iterdir() if x.is_file())


def check_folders(*folders: Path, logger: Logger):
    """Check if folders exist and are non-empty. Forces sys.exit if either is not true.

    Args:
        folders (Path): one or multiple path objects of folders

    Returns:
        None
    """
    for folder in folders:
        if not folder.is_dir():
            logger.error(
                f" Error in check_folders function: the folder {folder} does not exist"
            )
            sys.exit(f"The folder {folder} does not exist")
        if not any(folder.iterdir()):
            logger.error(f"The folder {folder} does not exist")
            sys.exit(f"Error in check_folders function: the folder {folder} is empty")


def extract_number_and_primers(
    file_path: Path, folder: Path, logger: Logger
) -> tuple[int, str, str, str, int]:
    """
    Extract the primer number, primer sequences, and amplicon length from a file.

    Args:
        file_path (Path): the filepath of the primerset in question (to get the number; generated by Primer3 module)
        folder (Path): the path object with the in silico testing results files

    Returns:
        a tuple with the identification number of the primer, the primers and the amplicon length
    """
    # the number is our identifier for the primer
    number = extract_number_from_filename(file_path.name, logger)

    # extract primer sequence
    pr_frwd, pr_rev, pr_intern = extract_primer_sequences(file_path, logger)

    # get the amplicon length from the amplified sequences in the in silico files
    target_seq = f"*Primer_{number}.txt_seqkit_amplicon_against_target_m0.txt"
    file = next(folder.glob(target_seq), None)
    if file is None:
        logger.error(
            "No file with in silico results for the target at m0 found.", exc_info=1
        )
        raise FileNotFoundError(f"No file matches the pattern {target_seq}")
    ampli_len = get_amplicon_length_from_seq(file)

    return number, pr_frwd, pr_rev, pr_intern, ampli_len


def extract_number_from_filename(filename: str, logger: Logger) -> int:
    """
    Extract the primer number from the file name.

    Args:
        filename (str): the file's name

    Returns:
        an integer that is the identifier of the primer and all its associated results
    """

    # regex match the number
    match = re.search(r"_(\d+)\.txt", filename)
    # if match is found, return the number
    if match:
        return int(match.group(1))
    else:
        logger.error(f"No number found in the filename: {filename}", exc_info=1)
        raise ValueError(f"No number found in the filename: {filename}")


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


def get_amplicon_length_from_seq(file: Path) -> int:
    """
    Extract the amplicon length from the sequence file.

    Args:
        file (Path): a path object of the m0 mismatch target in silico test file (seqkit amplicon results for Sensitivity testing at m0)
        
    Returns:
        the integer that is the amplicon length"""
    # open the file, read the first line, extract amplicon from the 7th column of the bed file
    with file.open("r") as f:
        first_line = f.readline()
        amplicon = first_line.strip().split("\t")[6]
        return len(amplicon)


def run_tests(
    folder: Path, file: str, length: int, count: int, target_type: str, logger: Logger
) -> dict:
    """
    Parse the Specificity and Sensitivity test results based on the in silico PCRs.

    Args:
        folder (Path): a path object corresponding to the folder with the in silico results (in the FUR.P3.PRIMERS) folder
        file (str): the primer file
        length (int): the length of the amplicon
        count (int): the number of targets or neighbours
        target_type (str): "target" or "neighbour"

    Returns:
        a dictionary with the results for targets and neighbours (specificity and sensitivity). Raw numbers not interpreted.
    """
    logger.info(
        "Starting the run_tests function and parsing the information of the in silico tests."
    )
    # find the in silico testing files
    file_pattern = f"{file}_seqkit_amplicon_against_{target_type}_m*"
    files = list(folder.glob(file_pattern))
    # initialize counters
    amp_f_overall, passed, failed, passed_n, failed_n = 0, 0, 0, 0, 0
    results_dict = {}
    # default if for instance Specificity tests do not amplify anything in silico (= no files)
    if not files:
        logger.info(
            f"No files with {file_pattern} found. Setting dictionary entries to NA."
        )
        results_dict["NA"] = {
            "Mismatches tested": "NA",
            "Did the test pass?": "NA",
            "Number of assemblies, in silico PCR was performed on": "NA",
            "Number of assemblies with correct size amplicon": "NA",
            "Number of assemblies with wrong size amplicon": "NA",
        }
        results_dict["Number of files that passed:"] = "NA"
        results_dict["Number of files that failed:"] = "NA"
        return results_dict
    # how many files are there
    no_files = len(files)
    for doc in files:
        # for each file initialize
        try:
            # find out the number of mismatches based on file name
            mismatch = re.search(r"_(m\d)\.txt", doc.name)
            if mismatch:
                m_no = mismatch.group(1)
                # print(f"Found match: {m_no}")
            else:
                logger.error("Could not determine value of mismatches")
                raise ValueError("No match found in the document name.")

        # except Exception as e:
        #     raise
        except ValueError as ve:
            logger.exception(
                f"Value error when trying to find matching document name: {ve}",
                exc_info=1,
            )
            raise ValueError(
                f"Value error when trying to find matching document name: {ve}"
            ) from ve
        passed_amp = 0
        failed_amp = 0
        try:
            with doc.open("r") as fp:
                lines = fp.readlines()
                # have all assemblies yielded one amplicon?
                if len(lines) == count:
                    # Process each line to count correct and incorrect amplicon lengths
                    for line in lines:
                        split_line = line.strip().split("\t")
                        if len(split_line) < 7:
                            logger.warning(f"Skipping line with insufficient columns: {split_line}")
                            continue  # Skip lines that do not have enough columns
                        amp_len = len(split_line[6])
                        if amp_len == length:
                            passed_amp += 1
                        else:
                            failed_amp += 1
                        amp_f_overall += 1
                    # Check if all assemblies have correct-length amplicons
                    if passed_amp == count:
                        passed += 1
                        failed_n += 1  # For specificity, correct counts in neighbor tests are a failure
                    if failed_amp > 0:
                        failed += 1  # for sensitivity, this is a failure
                        passed_n += 1  # for specificity, thisnis what we want
                else:
                    # If count mismatches, adjust pass/fail logic accordingly
                    failed += 1 if target_type == "target" else 0
                    passed_n += 1 if target_type == "neighbour" else 0
                    for line in lines:
                        split_line = line.strip().split("\t")
                        if len(split_line) < 7:
                            logger.warning(f"Skipping line with insufficient columns: {split_line}")
                            continue  # Skip lines that do not have enough columns
                        amp_len = len(split_line[6])
                        if amp_len == length:
                            passed_amp += 1
                        else:
                            failed_amp += 1
                        amp_f_overall += 1
                    if target_type == "neighbour" and passed_amp > 0:
                        passed_n -= 1
                        failed_n += 1
        except OSError as e:
            logger.exception(f"Unable to open or read file {doc}: {e}")
            raise OSError("Unable to open file") from e

        # Define whether it is Specificity or Sensitivity test and fail or pass test accordingly
        if target_type == "target":
            # only passed if all targets have been amplified to the correct length
            note = "passed" if count == passed_amp else "failed"
            results_dict[doc.name] = {
                "Mismatches tested": m_no,
                "Did the test pass?": note,
                "Number of assemblies, in silico PCR was performed on": count,
                "Number of assemblies with correct size amplicon": passed_amp,
                "Number of assemblies with wrong size amplicon": failed_amp,
            }
        else:
            # only pass if no neighbours have been amplified with the correct size amplicon.
            note = "passed" if (passed == 0 and passed_amp == 0) else "failed"
            results_dict[doc.name] = {
                "Mismatches tested": m_no,
                "Did the test pass?": note,
                "Number of assemblies, in silico PCR was performed on": count,
                "Number of assemblies with correct size amplicon": passed_amp,
            }
    results_dict[
        "Number of files (in silico PCR result files for different number of primer mismatches) tested:"
    ] = no_files
    if target_type == "target":
        results_dict["Number of files that passed:"] = passed
        results_dict["Number of files that failed:"] = failed
    else:
        results_dict["Number of files that passed:"] = passed_n
        results_dict["Number of files that failed:"] = failed_n
    logger.info("End of run_tests function. Returning dictionary.")
    return results_dict


def handle_blasts_and_efetch(
    destination_folder_tar: Path, number: int, logger: Logger
) -> str:
    """
    Handle BLAST results and use efetch to get the highest scoring accession's title (description).

    Args:
        destination_folder_tar (Path): the FUR.P3.TARGETS folder
        number (int): the identification number of the primer

    Returns:
        a string with information on the highest scoring target
    """
    logger.info(
        "running handle_blasts_and_efetch function to find highest scoring blast hits and determine what they are or run seqkit locate to generate bed files. "
    )

    # get the blast result
    target_f = f"*Target_{number}.txt_blastx_1e-5.txt"
    found_target = list(destination_folder_tar.glob(target_f))

    # 1: get the file
    if not found_target:
        logger.info(f"No files for {target_f} were found")
        return f"No files for {target_f} found."
    if any(file.stat().st_size == 0 for file in found_target):
        logger.info(
            f"One or more files for {target_f} are empty. Blastx did not return result"
        )
        return (
            f"One or more files for {target_f} are empty. Blastx did not return result"
        )

    # 2: get highest scoring accession
    try:
        id, evalue, bitscore = get_highest_scoring_accession(found_target[0])
        evalue = format(evalue, ".2e")
        if id == "Blast did not find hits":
            return f"The primer {number}'s BLASTX search did not find any hits.\n"

    # Error handling
    except ValueError as ve:
        logger.error(
            f"Value error while finding highest scoring accession: {ve}", exc_info=1
        )
        raise ValueError(
            f"Value error while finding highest scoring accession: {ve}"
        ) from ve
    except TypeError as te:
        logger.error(
            f"Type error while finding highest scoring accession: {te}", exc_info=1
        )
        raise TypeError(
            f"Type error while finding highest scoring accession: {te}"
        ) from te
    except FileNotFoundError as fnf:
        logger.error(f"File not found error: {fnf}", exc_info=1)
        raise FileNotFoundError(f"File not found error: {fnf}") from fnf
    except Exception as e:
        logger.error(
            f"Unknown exception occured.Could not find highest scoring accession for blastx results: {e}",
            exc_info=1,
        )
        raise Exception(
            f"Could not find highest scoring accession for blastx results: {e}"
        ) from e

    # 3: make sure eutilities are in path
    if is_tool("efetch"):

        # if eutilities are in path, find the annotation of the accession
        try:
            efetch_r = subprocess.Popen(
                ["efetch", "-db", "protein", "-id", id, "-format", "docsum"],
                stdout=subprocess.PIPE,
                text=True,
            )
            efetch_r.wait()
        except subprocess.CalledProcessError as e:
            logger.exception(f"efetch command failed: {e}", exc_info=1)
            raise RuntimeError(f"efetch command failed: {e}") from e
        try:
            xtract_r = subprocess.Popen(
                ["xtract", "-pattern", "DocumentSummary", "-element", "Title"],
                stdin=efetch_r.stdout,
                stdout=subprocess.PIPE,
                text=True,
            )
            xtract_r.wait()
        except subprocess.CalledProcessError as e:
            logger.error(f"xtract command failed: {e}", exc_info=1)
            raise RuntimeError(f"xtract command failed: {e}") from e
        title, _ = xtract_r.communicate()
        accession = "https://www.ncbi.nlm.nih.gov/protein/" + id + "/"
        return f"The primer {number}'s target with the highest bitscore {bitscore} & evalue {evalue} has the accession {id} and codes for {title.strip()}. \nFurther information on the accession: {accession} \n"
    else:
        accession = "https://www.ncbi.nlm.nih.gov/protein/" + id + "/"
        return f"The primer {number}'s target with the highest bitscore {bitscore} & evalue {evalue} has the accession {id}. \nFurther information on the accession: {accession} \n"


def get_highest_scoring_accession(blastx_file: Path) -> tuple[str, float, int]:
    """
    Parse a BLASTX output file to find the highest scoring accession and its evalue.

    Args:
        blastx_file (Path): a path object of the blastx results file for the target in question.

    Returns:
        a tuple of the Accession of the highest scoring entry, its evalue and its bitscore
    """
    if not blastx_file.exists() or blastx_file.stat().st_size == 0:
        raise FileNotFoundError(
            f"Blastx did not find hits. No file with the suffix {blastx_file} exists."
        )
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    try:
        df = pd.read_csv(blastx_file, sep="\t", comment="#", names=columns)
        if df.empty:
            raise ValueError(
                f"The BLASTX file {blastx_file} is empty or does not contain expected data."
            )
        highest_scoring_row = df.loc[df["bitscore"].idxmax()]
        return (
            highest_scoring_row["sseqid"],
            highest_scoring_row["evalue"],
            highest_scoring_row["bitscore"],
        )
    except:
        return "Blast did not find hits", float("inf"), 0


def interpret_and_reformat_sensi_speci_tests(
    test_res: dict, flavour: str, logger: Logger
) -> tuple[str, int, str]:
    """
    Interpret the Sensitivity and Specificity Tests, reformat the interpretation to a easy to grasp PASS or FAIL answer, make this \
    printable

    Args:
        test_res (dict): a dictionary of the test results for sensitivity or specificity testing
        flavour (str): target or neighbour

    Returns:
        a tuple with a pass or fail verdict (str), the number of processed assemblies (int) and the mismatches the test failed for (str; a concatenated string separated by colons)  
    """
    passed = "PASSED"
    ass_no = 0
    fail_m = []
    found_dict = False
    try:
        if not isinstance(test_res, dict):
            raise TypeError(
                f"Expected a dictionary, but got {type(test_res).__name__} instead."
            )
    except TypeError as e:
        print(e)
        raise
    # set flag to see if nested dictionary is true
    found_dict = False
    # Filtering keys to only process those that map to dictionaries
    for key in test_res.keys():
        result = test_res[key]

        # Skip keys that do not map to dictionaries
        if not isinstance(result, dict):
            continue
        # flag whether there is any dictionaries for result set to TRUE
        found_dict = True
        if flavour == "target":
            if result["Did the test pass?"] == "passed":
                ass_no = result["Number of assemblies, in silico PCR was performed on"]
                continue
            elif (
                result["Mismatches tested"] == "m3"
            ):  # if it is sensitivity testing, we tolerate that not all amplicons are the correct size for m3 or above
                fail = result["Mismatches tested"]
                fail_m.append(fail)
                ass_no = result["Number of assemblies, in silico PCR was performed on"]
                continue
            else:
                passed = "FAILED"
                fail = result["Mismatches tested"]
                fail_m.append(fail)
                ass_no = result["Number of assemblies, in silico PCR was performed on"]

        else:
            if (
                result["Did the test pass?"] == "passed"
                or result["Did the test pass?"] == "NA"
            ):
                ass_no = result["Number of assemblies, in silico PCR was performed on"]
                fail = result["Mismatches tested"]
                fail_m.append(fail)
                continue
            else:
                passed = "FAILED"
                fail = result["Mismatches tested"]
                fail_m.append(fail)
                ass_no = result["Number of assemblies, in silico PCR was performed on"]
    if not found_dict:
        logger.error(
            f"No dicitionary was found in the iteration over test_res. This is not a dictionary but a {type(found_dict)}"
        )
        raise TypeError("No dictionary was found in the iteration over test_res.")

    fail_m = ", ".join(fail_m)

    return passed, ass_no, fail_m


def print_results(
    header: str,
    primer_frwd: str,
    primer_rev: str,
    primer_internal: str,
    sensi_pass: str,
    sensi_ass_no: int,
    sensi_m: list,
    speci_pass: str,
    speci_ass_no: int,
    speci_m: list,
    sensi: dict,
    speci: dict,
    res_target_str: str,
    primer_fold: Path,
    target_fold: Path,
    seqkit_fold: Path,
    source_folder: Path,
    qPCR: str,
    ampli_len: int,
    logger: Logger,
):
    """
    Does what it says on the tin:
    Print results to a text file and a PDF (currently not implemented). Also generate an html.

    Args:
        header (str): the header with type of primer and file name of primer
        primer_frwd (str): the forward primer
        primer_rev (str): the reverse primer
        primer_internal (str): the internal probe (qPCR only, otherwise set to NA)
        sensi_pass (str): did it pass the sensitivity tests? (either "passed" or "failed")
        speci_pass (str): did it pass the specificity tests? (either "passed" or "failed")
        speci_ass_no (int): number of assemblies tested for specificity
        speci_m (str): the number of mismatches allowed in the test for seqkit
        sensi (dict): the dictionary with all sensitivity testing results
        speci (dict): the dictionary with all the specificity testing results
        res_target_str (str): results of the blastx testing
        primer_fold (Path): the path object of the folder the primers are in
        target_fold (Path): the Path object of the folder the targets are in
        seqkit_fold (Path): the Path object of the folder with the in silico testing results
        source_fold (Path): the Path object of the folder the entire thing is in (parent folder of all DiPPER2 results)
        qPCR (str): toggle of y or n to change formatting according to whether it is qPCR or not
    Returns:
        None
    """
    # txt outfile
    result_out = source_folder / "Results.txt"
    logger.info(f"Writing results to the results file Results.txt in {source_folder}")
    with open(result_out, "a", encoding="utf-8") as file:
        file.write(header)
        file.write(primer_frwd)
        file.write(primer_rev)
        if qPCR == "y":
            file.write(primer_internal)
        else:
            pass
        file.write(f"\nAmplicon length: {ampli_len}\n")
        file.write("Sensitivity and Specificity testing:\n")
        file.write("Sensitivity:\n")
        file.write(f"Pass or fail:\t{sensi_pass}\n")
        file.write(
            f"Number of assemblies the in silico PCR was run for:\t {sensi_ass_no}\n"
        )
        file.write(
            f"Number of mismatches the sensitivity test failed for (m3 is tolerated):\t{sensi_m}\n"
        )
        file.write("Specificity:\n")
        file.write(f"Pass or fail:\t{speci_pass}\n")
        file.write(
            f"Number of assemblies the in silico PCR was run for:\t {speci_ass_no}\n"
        )
        file.write(
            f"Number of mismatches the in silico PCR generated amplicons for:\t{speci_m}\n\
                   (Note that if the test has passed, this means the amplicons had an incorrect size)\n"
        )
        file.write("BLASTX Target Testing:\n")
        file.write(res_target_str)
        file.write("Files&Folders:\n")
        file.write(f"Primers are found in {primer_fold}.\n")
        file.write(f"Target fastas and blastx results are found in {target_fold}.\n")
        file.write(f"If generated, bed files are present in {source_folder}.\n")
        file.write(f"In silico PCR results are found in {seqkit_fold}.\n")
        file.write("Detailed Sensitivity results:\n")
        # prettyprint allows to print a nested dictionary
        pprint(sensi, stream=file, depth=4)
        file.write("Detailed Specificity results:\n")
        pprint(speci, stream=file, depth=4)
        file.write(
            "\n################################################################################\n"
        )
    logger.info(f"Done!\n Generating html...")
    generate_html_jinja(
        header,
        primer_frwd,
        primer_rev,
        primer_internal,
        sensi_pass,
        sensi_ass_no,
        sensi_m,
        speci_pass,
        speci_ass_no,
        speci_m,
        sensi,
        speci,
        res_target_str,
        primer_fold,
        target_fold,
        seqkit_fold,
        source_folder,
        qPCR,
        ampli_len,
    )


def generate_results(
    destination_folder_primer: Path,
    destination_folder_tar: Path,
    destination_seqkit: Path,
    file_path: Path,
    count_target: int,
    count_neighbour: int,
    source_folder: Path,
    qPCR: str,
    logger: Logger,
):
    """
    Gets the primers & amplicon length, processes them for output, gets the results from the sensitivity and specificity tests
    and finally also gets the blast results. It then calls the print results function to generate an output file.

    Args:
        destination_folder_primer (Path): the folder that stores primers and primer data, as well as in silico tests
        destination_folder_target (Path): the folder that stores the target fastas and the blastx results
        destination_folder_seqkit (Path): subfolder of the primer folder where the in silico results are stored
        file_path (Path): a primer fasta file
        count_target (int): the number of targets
        count_neighbour (int): the number of neighbours
        source_folder (Path): one folder to rule them all (the folder with all DiPPER results of this run)
        qPCR (str): a toggle whether qPCR primers are wanted and have been generated or not

    Returns:
        None
    """
    logger.info("Retrieving Primers and obtaining amplicon length...")
    try:
        number, pr_frwd, pr_rev, pr_intern, ampli_len = extract_number_and_primers(
            file_path, destination_seqkit, logger
        )
    except Exception as e:
        logger.error(f"Error processing {file_path}:{e}", exc_info=1)
        raise Exception(f"Error processing {file_path}: {e}") from e
    if qPCR == "y":
        header = f"RESULTS FOR qPCR PRIMER {number} in {file_path.name}\n"
    else:
        header = f"RESULTS FOR CONV. PCR PRIMER {number} in {file_path.name}\n"
    primer_frwd = f">Forward Primer {number}\n{pr_frwd}"
    primer_rev = f">Reverse Primer {number}\n{pr_rev}"
    primer_internal = f">Internal Probe {number}\n{pr_intern}\n"

    logger.info(
        "Assessing the in silico PCR results and determining Specificity and Sensitivity..."
    )
    try:
        sensi = run_tests(
            destination_seqkit,
            file_path.name,
            ampli_len,
            count_target,
            "target",
            logger,
        )
    except Exception as e:
        logger.error(f"Sensitivity tests failed: {e}", exc_info=1)
        raise Exception(f"Sensitivity tests failed: {e}") from e

    try:
        speci = run_tests(
            destination_seqkit,
            file_path.name,
            ampli_len,
            count_neighbour,
            "neighbour",
            logger,
        )
    except Exception as e:
        logger.error(f"Specificity tests failed: {e}", exc_info=1)
        raise Exception(f"Specificity tests failed: {e}") from e

    logger.info(
        "Retrieving the blastx results of the targets, if entrez direct (eutils) is installed, also check NCBI annotation..."
    )
    try:
        res_target_str = handle_blasts_and_efetch(
            destination_folder_tar, number, logger
        )
    except Exception as e:
        logger.error(f"BLASTX Target Testing failed: {e}", exc_info=1)
        raise Exception(f"BLASTX Target Testing failed: {e}") from e

    results_folder = os.getcwd()
    # interpret the tests
    # Sensi
    sensi_pass, sensi_ass_no, sensi_m = interpret_and_reformat_sensi_speci_tests(
        sensi, "target", logger
    )
    # Speci
    speci_pass, speci_ass_no, speci_m = interpret_and_reformat_sensi_speci_tests(
        speci, "neighbour", logger
    )
    # Print results to txt file and html
    logger.info(
        f"Printing results to outfile. The results summary files 'Results.txt' & 'Results.html' is found in {results_folder}."
    )
    print_results(
        header,
        primer_frwd,
        primer_rev,
        primer_internal,
        sensi_pass,
        sensi_ass_no,
        sensi_m,
        speci_pass,
        speci_ass_no,
        speci_m,
        sensi,
        speci,
        res_target_str,
        destination_folder_primer,
        destination_folder_tar,
        destination_seqkit,
        source_folder,
        qPCR,
        ampli_len,
        logger,
    )


#####################
######  MAIN  #######
#####################
def main():
    """
    The main function. This is where the magic happens. Master of all puppets (functions). Also parses arguments, initializes logger, and calls all functions. Tadaa!
    """
    parser = argparse.ArgumentParser(
        prog="DiPPER2 Summarize Results Module",
        description="Summarize the results and generate html results",
        epilog="Theresa Wacker T.Wacker2@exeter.ac.uk",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.1.0-beta")
    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        required=True,
        help="Results folder name, which includes results folders from previous steps",
    )
    parser.add_argument(
        "-q",
        "--qPCR",
        default="n",
        type=str,
        help="Toggles between qPCR (y) or no qPCR (n) results. Default is 'n'",
    )
    parser.add_argument(
        "-V", "--verbose", action="store_true", help="increase logging verbosity"
    )
    # get args
    args = parser.parse_args()

    # one folder to rule them all
    source_folder = Path(args.folder)

    # configures the logger
    module_name = Path(__file__).name
    logger_instance = Logger(module_name, source_folder, args.verbose)
    logger = logger_instance.get_logger()
    logger.info("Starting to run the Summarize_results module. Logger is initialized.")
    # toggle qPCR
    qPCR = args.qPCR

    # Define folders
    destination_folder_pr = source_folder / "FUR.P3.PRIMERS"
    destination_folder_tar = source_folder / "FUR.P3.TARGETS"
    destination_folder_seqkit = source_folder / "FUR.P3.PRIMERS" / "in_silico_tests"

    # do they exist and are they not empty?
    check_folders(
        destination_folder_pr, destination_folder_tar, destination_folder_seqkit, logger=logger
    )

    # target and neighbour folders
    fur_target = source_folder / "FUR.target"
    fur_neighbour = source_folder / "FUR.neighbour"

    # do they exist and are they not empty?
    check_folders(fur_target, fur_neighbour, logger=logger)

    # define how many targets (equivalent to number of target assemblies, each assembly in one file) and neighbours (ditto targets) there are
    count_target = count_files(fur_target)
    count_neighbour = count_files(fur_neighbour)

    # info
    logger.info(
        "%s has %d entries and %s has %d entries.",
        fur_target,
        count_target,
        fur_neighbour,
        count_neighbour,
    )

    # to be able to loop through each primer of the 4 candidates, find all the files and generate a list of paths (it is a generator object and yields Path objects with name and path attributes)
    all_files = list(destination_folder_pr.glob("*"))

    # generate results
    for file_path in all_files:
        if file_path.is_file():
            generate_results(
                destination_folder_pr,
                destination_folder_tar,
                destination_folder_seqkit,
                file_path,
                count_target,
                count_neighbour,
                source_folder,
                qPCR,
                logger,
            )


if __name__ == "__main__":
    main()
