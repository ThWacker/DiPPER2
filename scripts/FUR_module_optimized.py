#!/usr/bin/python

import sys
import argparse
import subprocess
from datetime import datetime
from pathlib import Path
import shutil
from logging_handler import Logger

def check_folders(*folders: Path,logger: Logger):
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

def make_fur_db(target_folder: Path, neighbour_folder: Path, source: Path, logger: Logger, reference=None):
    """
    Create the FUR database.
    
    Args:
        target_folder(Path): the folder with the target accessions
        neighbour_folder(Path): the folder with the neighbour accessions
        source(Path): one folder to rule them all (the folder with all results from DiPPER2 run)
    
    Returns:
        None
    
    Raises: 
        RunTimeError
    """
    try:
        logger.info("Making FUR DB.")
        fur_db = source / 'FUR.db'
        if reference:
            subprocess.run(['makeFurDb', '-t', str(target_folder), '-n', str(neighbour_folder), '-d', str(fur_db), '-r', str(reference)], check=True)
        else:
            subprocess.run(['makeFurDb', '-t', str(target_folder), '-n', str(neighbour_folder), '-d', str(fur_db)], check=True)
    except subprocess.CalledProcessError as e:
        logger.exception(f"Error when executing makeFurDb: {e}")
        raise RuntimeError(f"Error when executing makeFurDb: {e}") from e

def run_fur(option: str, outfile_prefix: str, source: Path, logger: Logger):
    """
    Run the FUR command and handle the output.
    
    Args: 
        option(str): one of the flags from FUR (u, U or m)
        outfile_prefix(str): the outfile prefix
        source(Path): one folder to rule them all (the folder with all results from this DiPPER2 run)
    
    Returns:
        The StdErr output of FUR, which contains basic stats. 
   
    Raises:
        RunTimeError
    """
    option_flag = f"-{option.strip()}" if option.strip() else ""
    fur_out = source / f"{outfile_prefix}_FUR.db.out.txt"
    fur_db = source / 'FUR.db'
    try:
        logger.info("Running FUR.")
        result = subprocess.run(['fur', '-d', str(fur_db), option_flag], check=True, capture_output=True, text=True)
        fur_output = result.stdout.strip()
        error=result.stderr.strip()
    except subprocess.CalledProcessError as e:
        logger.exception(f"FUR did not run to completion and a Runtime Error occured: {e}")
        raise RuntimeError(f"FUR did not run to completion and a Runtime Error occured: {e}") from e

    fur_out.write_text(fur_output, encoding="utf-8" )

    if not fur_out.exists() or fur_out.stat().st_size == 0:
        logger.error(f"Fur could not find unique regions or did not run successfully. {fur_out} is empty or does not exist.")
        sys.exit()

    return error

def clean_up(db, logger: Logger):
    """
    Clean up temporary files.
    
    Args:
        db(Path): the FUR database
    
    Returns:
        None
    """
    if db.exists():
        shutil.rmtree(db)
        logger.info(f"The directory {db} has been deleted.", file=sys.stderr)

def main():
    """
    The main function of the FUR module that is called by the wrapper. 
    Requires a folder, has defaults for everything else. 
    """
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y_%Hh%Mmin%Ss_%z")
    #Parse all arguments
    parser = argparse.ArgumentParser(prog="DiPPER2",description='This module runs FUR to find unique genomic regions', epilog="Bugs, suggestions, criticism, chocolate, puppies and praise to t.wacker2@exeter.ac.uk")
    parser.add_argument('-p', '--parameter', type=str, default=" ", help='Tells FUR how to run. -u -> only the first subtraction step; -U -> runs first Subtraction step and Intersection step; -m -> runs megablast instead of blastn in the last step; without argument runs FUR in full.')
    parser.add_argument('-o', '--outfile_prefix', default=dt_string, type=str, help='Outfile prefix. Default is date and time in d-m-y-h-m-s-tz format')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')
    parser.add_argument('-f', '--folder', type=str, required=True, help='Results folder name, which includes results folders from previous steps')
    parser.add_argument('-r', '--reference', type=str, default=None,help="Optional reference for the makeFurDb step" )

    args = parser.parse_args()
    
    # find the source folder, define the FUR.target and FUR.neighbour folders
    source_folder = Path(args.folder)
    target_folder = source_folder / "FUR.target"
    neighbour_folder = source_folder / "FUR.neighbour"
    
    # configures the logger
    module_name = Path(__file__).name
    logger_instance = Logger(module_name, source_folder, args.verbose)
    logger = logger_instance.get_logger()

    # Basic info
    logger.info(f'The parameters are:\n for the FUR parameter {args.parameter}\n for the outfile prefix: {args.outfile_prefix}', file=sys.stderr)
    
    #check if they exist
    check_folders(target_folder, neighbour_folder, logger=logger)

    # run makefurdb
    if args.reference:
        try:
            reference= args.reference
            make_fur_db(target_folder, neighbour_folder, source_folder,logger=logger, reference=reference)
        except (FileNotFoundError, OSError):
            logger.info(f"Could not find file {reference}. Defaulting to running without a reference")
            make_fur_db(target_folder, neighbour_folder, source_folder, logger)
    else:
        make_fur_db(target_folder, neighbour_folder, source_folder, logger)

    #run fur
    fur_out = run_fur(args.parameter, args.outfile_prefix, source_folder, logger)

    #print the stats to stderr & save in file(might not want that on server)
    print(f'{fur_out}')
    fur_file= source_folder/ "FUR_Summary_output.txt"
    try:
        with open(fur_file, 'w', encoding="utf-8") as f:
            f.write(fur_out)
    except OSError as e:
        logger.exception(f"Could not open or write the FUR summary output to {fur_file}", {e})
        raise OSError(f"Could not open or write the FUR summary output to {fur_file}") from e



    # clean up after yourself
    clean_up(source_folder / "FUR.db", logger)

    #Heureka
    logger.info('FUR_module.py ran to completion: exit status 0')

if __name__ == '__main__':
    main()
