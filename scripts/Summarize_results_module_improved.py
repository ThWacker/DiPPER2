#!/usr/bin/python

import sys
import argparse
from pathlib import Path
import subprocess
from shutil import which
import re
from pprint import pprint
import os
from pprint import pformat
import jinja2
import pandas as pd # type: ignore



def generate_html_jinja(header: str, primer_frwd: str, primer_rev: str, primer_internal: str,sensi_pass:str,sensi_ass_no:int, sensi_m:str, speci_pass:str,\
                   speci_ass_no:int, speci_m:str, sensi: dict, speci: dict, res_target_str: str, primer_fold: Path, target_fold: Path, seqkit_fold: Path, source_folder:Path,qPCR:str):
    """Generate HTML content for the results."""
    HTML_OUTFILE = source_folder/ "Results.html"
    script_dir = os.path.dirname(os.path.realpath(__file__))
    environment = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.join(script_dir, "templates/")))
    if qPCR == "y":
        template = environment.get_template("results_qPCR.html")
    else:
        template = environment.get_template("results.html")
    content = template.render(
        header=header,
        primer_frwd=primer_frwd,
        primer_rev=primer_rev,
        primer_internal=primer_internal,
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
        source_folder=source_folder)

    with open(HTML_OUTFILE, "a", encoding="utf-8") as stream:
        stream.write(content)


def quit_program(message: str = None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    usage()
    sys.exit(1)

def is_tool(name: str) -> bool:
    """Check whether `name` is on PATH and marked as executable."""
    #which is a tool in shutil. 
    return which(name) is not None

def usage():
    """Display usage information."""
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Theresa Wacker T.Wacker2@exeter.ac.uk')
    print('')
    print('Module/ Script that summarizes results')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-f/--folder - results folder name, which includes results folders from previous steps')
    print('optional: -v/--version for the version & -q/--qPCR to toggle qPCR (y) or no qPCR (n) results (default: n)')
    print('CAUTION: needs to be in the folder with the assemblies to be sorted in target and neighbour')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')

def count_files(folda: Path) -> int:
    """Counts the number of files in a folder."""
    return sum(1 for x in folda.iterdir() if x.is_file())

def check_folders(*folders: Path):
    """Check if folders exist and are non-empty."""
    for folder in folders:
        if not folder.is_dir():
            sys.exit(f"The folder {folder} does not exist")
        if not any(folder.iterdir()):
            sys.exit(f"The folder {folder} is empty")

def extract_number_and_primers(file_path: Path, folder:Path) -> tuple[int, str, str, str, int]:
    """Extract the primer number, primer sequences, and amplicon length from a file."""
    number = extract_number_from_filename(file_path.name)
    pr_frwd, pr_rev, pr_intern = extract_primer_sequences(file_path)
    #NEED TO ADD HOW TO GET TO THE SEQKIT FILE AND EXTRACT THE AMPLICON LENGTH
    target_seq=f"*Primer_{number}.txt_seqkit_amplicon_against_target_m0.txt"
    file = next(folder.glob(target_seq), None)
    if file is None:
        raise FileNotFoundError(f"No file matches the pattern {target_seq}")
    ampli_len = get_amplicon_length_from_seq(file)

    return number, pr_frwd, pr_rev, pr_intern, ampli_len

def extract_number_from_filename(filename: str) -> int:
    """Extract the primer number from the file name."""
    match = re.search(r'_(\d+)\.txt', filename)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"No number found in the filename: {filename}")

def extract_primer_sequences(file: Path) -> tuple[str, str, str]:
    """Extract primer sequences from a file."""
    sequences = {"PRIMER_LEFT": None, "PRIMER_RIGHT": None, "PRIMER_INTERNAL": None}
    with file.open('r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            for key in sequences.keys():
                if key in line:
                    sequences[key] = lines[i + 1].strip()
    if not all(sequences.values()):
        quit_program("The primer headers are not correctly formatted and cannot be processed. Please change the headers accordingly.")
    return sequences["PRIMER_LEFT"], sequences["PRIMER_RIGHT"], sequences["PRIMER_INTERNAL"]

def get_amplicon_length_from_seq(file: Path) -> int:
    """Extract the amplicon length from the sequence file."""
    with file.open('r') as f:
        first_line = f.readline()
        amplicon = first_line.strip().split('\t')[6]
        return len(amplicon)

def run_tests(folder: Path, file: str, length: int, count: int, target_type: str) -> dict:
    """Interpret the Specificity and Sensitivity test results based on the in silico PCRs"""
    file_pattern = f"{file}_seqkit_amplicon_against_{target_type}_m*"
    files = list(folder.glob(file_pattern))
      #initialize counters
    amp_f_overall,passed, failed = 0,0,0
    results_dict={}
    if not files:
        results_dict["NA"] = {
                "Mismatches tested":"NA",
                "Did the test pass?": "NA",
                "Number of assemblies, in silico PCR was performed on": "NA",
                "Number of assemblies with correct size amplicon": "NA",
                "Number of assemblies with wrong size amplicon": "NA"
        }
        results_dict["Number of files that passed:"]="NA"
        results_dict["Number of files that failed:"]="NA"
        return results_dict
    no_files=len(files)
    for doc in files:
        #for each file initialize
        try:
            mismatch = re.search(r'_(m\d)\.txt', doc.name)
            if mismatch:
                m_no = mismatch.group(1)
                #print(f"Found match: {m_no}")
            else:
                raise ValueError("No match found in the document name.")
                
        except AttributeError:
            print("Error: Tried to access a match group on a NoneType object.")
        except ValueError as ve:
            print(ve)
        passed_amp=0
        failed_amp=0
        try:
            with doc.open('r') as fp:
                lines = fp.readlines()
                if len(lines) == count:
                    #passed the test for correct number of assemblies, test whether all of them are the right amplicon length
                    for line in lines:
                        amp_len = len(line.strip().split('\t')[6])
                        if amp_len == length:
                            passed_amp += 1
                            
                        else:
                            failed_amp += 1
                        amp_f_overall+=1
                    #did we count as many correct length amplicons as assemblies
                    if passed_amp==count:
                        passed+=1
                    if not failed_amp==0:
                        failed+=1
                #next file after this iteration. passed or failed should only incriment per file
                else:
                    failed += 1
        except OSError as e:
            quit_program(f"could not open or read file: {e}")
    #Define whether it is Specificity or Sensitivity test and fail or pass test accordingly
        if target_type=="target":
            note = "passed" if count == passed_amp else "failed"
            results_dict[doc.name] = {
                "Mismatches tested": m_no,
                "Did the test pass?": note,
                "Number of assemblies, in silico PCR was performed on": count,
                "Number of assemblies with correct size amplicon": passed_amp,
                "Number of assemblies with wrong size amplicon": failed_amp
            }
        else:
            note = "passed" if passed==0 else "failed"
            results_dict[doc.name] = {
                "Mismatches tested": m_no,
                "Did the test pass?": note,
                "Number of assemblies, in silico PCR was performed on": count,
                "Number of assemblies with correct size amplicon": passed_amp
            }
    results_dict["Number of files (in silico PCR result files for different number of primer mismatches) tested:"]=no_files
    if target_type=="target":
        results_dict["Number of files that passed:"]=passed
        results_dict["Number of files that failed:"]=failed
    else:
        results_dict["Number of files that passed:"]=failed
        results_dict["Number of files that failed:"]=passed
    return results_dict

def handle_blasts_and_efetch(destination_folder_tar: Path, number: int) -> str:
    """Handle BLAST results and use efetch to get the highest scoring accession's title (description)."""
    target_f = f"*Target_{number}.txt_blastx_1e-5.txt"
    found_target = list(destination_folder_tar.glob(target_f))
    if not found_target:
        return f"No files for {target_f} found."
    if any(file.stat().st_size == 0 for file in found_target):
        return f"One or more files for {target_f} are empty."

    try:
        id, evalue, bitscore = get_highest_scoring_accession(found_target[0])
        evalue=format(evalue, ".2e")
    except Exception as e:
        quit_program(f"Could not find highest scoring accession for blastx results: {e}")

    if id == "Blast did not find hits":
        return f"The primer {number}'s BLASTX search did not find any hits.\n"
    elif is_tool("efetch"):
        try:
            efetch_r = subprocess.Popen(["efetch", "-db", "protein", "-id", id, "-format", "docsum"], stdout=subprocess.PIPE, text=True)
            xtract_r = subprocess.Popen(["xtract", "-pattern", "DocumentSummary", "-element", "Title"],
                                        stdin=efetch_r.stdout, stdout=subprocess.PIPE, text=True)
            title, _ = xtract_r.communicate()
            accession="https://www.ncbi.nlm.nih.gov/protein/"+id+"/"
            return f"The primer {number}'s target with the highest bitscore {bitscore} & evalue {evalue} has the accession \
                {id} and codes for {title.strip()}.\
                \nFurther information on the accession: {accession} \n"
        except Exception as e:
            quit_program(f"efetch/xtract command failed: {e}")
    else:
        accession="https://www.ncbi.nlm.nih.gov/protein/"+id+"/"
        return f"The primer {number}'s target with the highest bitscore {bitscore} & evalue {evalue} has the accession {id}.\
        \nFurther information on the accession: {accession} \n"

def get_highest_scoring_accession(blastx_file: Path) -> tuple[str, float, int]:
    """Parse a BLASTX output file to find the highest scoring accession and its evalue."""
    if not blastx_file.exists() or blastx_file.stat().st_size == 0:
        raise FileNotFoundError(f"Blastx did not find hits. No file with the suffix {blastx_file} exists.")
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    try:
        df = pd.read_csv(blastx_file, sep='\t', comment='#', names=columns)
        if df.empty:
            raise ValueError(f"The BLASTX file {blastx_file} is empty or does not contain expected data.")
        highest_scoring_row = df.loc[df['bitscore'].idxmax()]
        return highest_scoring_row['sseqid'], highest_scoring_row['evalue'], highest_scoring_row['bitscore']
    except:
        return "Blast did not find hits", float('inf'), 0

def interpret_and_reformat_sensi_speci_tests(test_res: dict, flavour: str) -> tuple[str, int, str]:
    passed = "PASSED"
    ass_no = 0
    fail_m = []

    # Filtering keys to only process those that map to dictionaries
    for key in test_res.keys():
        result = test_res[key]
        
        # Skip keys that do not map to dictionaries
        if not isinstance(result, dict):
            continue
        
        if flavour == "target":
            if result['Did the test pass?'] == "passed":
                ass_no = result['Number of assemblies, in silico PCR was performed on']
                continue
            elif result['Mismatches tested'] == "m3":
                fail_m = result['Mismatches tested']  # This seems incorrect, should it be fail_m.append(...)
                ass_no = result['Number of assemblies, in silico PCR was performed on']
                continue
            else:
                passed = "FAILED"
                fail = result['Mismatches tested']
                fail_m.append(fail)
                ass_no = result['Number of assemblies, in silico PCR was performed on']
            
        else:
            if result['Did the test pass?'] == "passed":
                ass_no = result['Number of assemblies, in silico PCR was performed on']
                fail = result['Mismatches tested']
                print(f"{fail}")
                fail_m.append(fail)
                print(f"{fail_m}")
                continue
            else:
                passed = "FAILED"
                fail = result['Mismatches tested']
                print(f"{fail}")
                fail_m.append(fail)
                print(f"{fail_m}")
                ass_no = result['Number of assemblies, in silico PCR was performed on']
    fail_m=', '.join(fail_m)

    return passed, ass_no, fail_m
    

def print_results(header: str, primer_frwd: str, primer_rev: str, primer_internal: str,sensi_pass:str,sensi_ass_no:int, sensi_m:list\
                  , speci_pass:str, speci_ass_no:int, speci_m:list, sensi: dict, speci: dict, res_target_str: str,primer_fold:Path, \
                    target_fold:Path, seqkit_fold:Path, source_folder:Path, qPCR:str):
    """Does what it says on the tin:
       Print results to a text file and a PDF (currently not implemented)."""
    #txt outfile
    result_out= source_folder / "Results.txt"
    with open(result_out, 'a', encoding="utf-8") as file:
        file.write(header)
        file.write(primer_frwd)
        file.write(primer_rev)
        if qPCR =="y":
            file.write(primer_internal)
        else:
            pass
        file.write("\n Sensitivity and Specificity testing:\n")
        file.write("\n Sensitivity:\n")
        file.write(f"Pass or fail:\t{sensi_pass}\n")
        file.write(f"Number of assemblies the in silico PCR was run for:\t {sensi_ass_no}\n")
        file.write(f"Number of mismatches the sensitivity test failed for (m3 is tolerated):\t{sensi_m}\n")
        file.write("\nSpecificity:\n")
        file.write(f"Pass or fail:\t{speci_pass}\n")
        file.write(f"Number of assemblies the in silico PCR was run for:\t {speci_ass_no}\n")
        file.write(f"Number of mismatches the in silico PCR generated amplicons for:\t{speci_m}\n\
                   (Note that if the test has passed, this means the amplicons had an incorrect size)\n")
        file.write("\n BLASTX Target Testing:\n")
        file.write(res_target_str)
        file.write("\n Files&Folders:\n")
        file.write(f"Primers are found in {primer_fold}.\n")
        file.write(f"Target fastas and blastx results are found in {target_fold}.\n")
        file.write(f"If generated, bed files are present in {source_folder}.\n")
        file.write(f"In silico PCR results are found in {seqkit_fold}.\n")
        file.write("Detailed Sensitivity results:\n")
        #prettyprint allows to print a nested dictionary
        pprint(sensi, stream=file, depth=4)
        file.write("Detailed Specificity results:\n")
        pprint(speci, stream=file, depth=4)
        file.write("\n################################################################################\n")

    #generate the html
    # html_content = generate_html(header, primer_frwd, primer_rev, primer_internal,sensi_pass,sensi_ass_no, sensi_m, speci_pass, speci_ass_no, \
    #               speci_m, sensi, speci, res_target_str, primer_fold, target_fold, seqkit_fold, source_folder)

    generate_html_jinja(header, primer_frwd, primer_rev, primer_internal,sensi_pass,sensi_ass_no, sensi_m, speci_pass, speci_ass_no, \
                  speci_m, sensi, speci, res_target_str, primer_fold, target_fold, seqkit_fold, source_folder, qPCR)
    
    
    # with open("Results.html", 'a', encoding="utf-8") as file:
    #     file.write(html_content)

def generate_results(destination_folder_primer, destination_folder_tar, destination_seqkit,file_path, count_target, count_neighbour, source_folder, qPCR):
    """Gets the primers & amplicon length, processes them for output, gets the results from the sensitivity and specificity tests 
       and finally also gets the blast results. It then calls the print results function to generate an output file"""
    print("Retrieving Primers and obtaining amplicon length...")
    try:
        number, pr_frwd, pr_rev, pr_intern, ampli_len = extract_number_and_primers(file_path, destination_seqkit)
    except Exception as e:
        quit_program(f"Error processing {file_path}: {e}")
    if qPCR=="y":
        header = f'RESULTS FOR qPCR PRIMER {number} in {file_path.name}\n\n'
    else:
        header = f'RESULTS FOR CONV. PCR PRIMER {number} in {file_path.name}\n\n'
    primer_frwd = f'>Forward Primer {number}\n{pr_frwd}\n'
    primer_rev = f'>Reverse Primer {number}\n{pr_rev}\n'
    primer_internal = f'>Internal Probe {number}\n{pr_intern}\n\n'

    print("Assessing the in silico PCR results and determining Specificity and Sensitivity...")
    try:
        sensi = run_tests(destination_seqkit, file_path.name, ampli_len, count_target, "target")
    except Exception as e:
        quit_program(f"Sensitivity tests failed: {e}")

    try:
        speci = run_tests(destination_seqkit, file_path.name, ampli_len, count_neighbour, "neighbour")
    except Exception as e:
        quit_program(f"Specificity tests failed: {e}")

    print("Retrieving the blastx results of the targets, if entrez direct (eutils) is installed, also check NCBI annotation...")
    try:
        res_target_str = handle_blasts_and_efetch(destination_folder_tar, number)
    except Exception as e:
        quit_program(f"BLASTX Target Testing failed: {e}")

    results_folder=os.getcwd()
    #interpret the tests
    # Sensi
    sensi_pass,sensi_ass_no, sensi_m=interpret_and_reformat_sensi_speci_tests(sensi, "target")
    # Speci
    speci_pass,speci_ass_no, speci_m=interpret_and_reformat_sensi_speci_tests(speci, "neighbour")
    #Print results to txt file and html
    print (f"Printing results to outfile. The results summary files 'Results.txt' & 'Results.html' is found in {results_folder}.")
    print_results(header, primer_frwd, primer_rev, primer_internal,sensi_pass,sensi_ass_no, sensi_m, speci_pass, speci_ass_no, \
                  speci_m, sensi, speci, res_target_str, destination_folder_primer, destination_folder_tar, destination_seqkit, source_folder, qPCR)


#####################
######  MAIN  #######
#####################
def main():

    parser = argparse.ArgumentParser(description='Primer testing module')
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which includes results folders from previous steps')
    parser.add_argument('-q', '--qPCR', default="n", type=str, help="Toggles between qPCR (y) or no qPCR (n) results. Default is 'n'" )

    #get args
    args = parser.parse_args()

    #one folder to rule them all
    source_folder = Path(args.folder)
    
    #toggle qPCR
    qPCR=args.qPCR

    # Define folders
    destination_folder_pr = source_folder / "FUR.P3.PRIMERS"
    destination_folder_tar = source_folder / "FUR.P3.TARGETS"
    destination_folder_seqkit= source_folder/ "FUR.P3.PRIMERS"/ "in_silico_tests"

    # do they exist and are they not empty?
    check_folders(destination_folder_pr, destination_folder_tar, destination_folder_seqkit)

    # target and neighbour folders
    fur_target = source_folder / "FUR.target"
    fur_neighbour = source_folder / "FUR.neighbour"

    # do they exist and are they not empty?
    check_folders(fur_target, fur_neighbour)

    # define how many targets (equivalent to number of target assemblies, each assembly in one file) and neighbours (ditto targets) there are
    count_target = count_files(fur_target)
    count_neighbour = count_files(fur_neighbour)

    #info
    print(f"{fur_target} has {count_target} entries and {fur_neighbour} has {count_neighbour} entries.")

    #to be able to loop through each primer of the 4 candidates, find all the files and generate a list of paths (it is a generator object and yields Path objects with name and path attributes)
    all_files = list(destination_folder_pr.glob('*'))

    #generate results
    for file_path in all_files:
        if file_path.is_file():
            generate_results(destination_folder_pr,destination_folder_tar, destination_folder_seqkit, file_path, count_target, count_neighbour, source_folder,qPCR)

  

if __name__ == '__main__':
    main()
