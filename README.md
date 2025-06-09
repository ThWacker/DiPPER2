[![CodeQL Advanced](https://github.com/ThWacker/DiPPER2/actions/workflows/codeql.yml/badge.svg)](https://github.com/ThWacker/DiPPER2/actions/workflows/codeql.yml) [![Python Lint and Code Quality](https://github.com/ThWacker/DiPPER2/actions/workflows/pylint.yml/badge.svg)](https://github.com/ThWacker/DiPPER2/actions/workflows/pylint.yml)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

 <!----
  (![Static Badge](https://img.shields.io/badge/unittest_code_coverage-78%25-yellow))
  --->

![GitHub last commit](https://img.shields.io/github/last-commit/ThWacker/DiPPER2)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14699018.svg)](https://doi.org/10.5281/zenodo.14699018)
![](./logo/Dipper2_white.png) 
[![CodeFactor](https://www.codefactor.io/repository/github/thwacker/dipper2/badge/main)](https://www.codefactor.io/repository/github/thwacker/dipper2/overview/main)
# DiPPER2
DiPPER2:  Diagnostic Primer Picking and Evaluation pipeline for Reliability and Reproducibility

>__+++++NOTE - IMPORTANT:+++++__
>
>__This is a work in progress.__
>
>__++++++++++++++++++++++++++__

## Synopsis
__*This pipeline and modules are meant to facilitate reliable and reproducible finding of diagnostic targets and to make picking primers for those targets as user-friendly as possible. The approach taken is a phylogeny-driven and clade-specific approach.*__ 

The pipeline automatically finds unique genomic regions for a target species, pathovar or race with moderate user intervention, automatically picks primers for these regions that are amenable to both conventional PCR and quantitative PCR, and finally tests them for specificity and sensitivity.

The pipeline can also be run "stand-alone" to *only test* but not pick primers. In that case it will only perform tests for specificity and sensitivity for existing primers.

## Overview
DiPPER2 is a tool to find unique genomic regions of target species, build primers (conv. or qPCR) for them, validate them *in silico*, define the identity of the unique genomic regions (either by homology or genome coordinates) and produce a user-friendly html report (and a more machine-readable .txt report).

Standalone-DiPPER2 is modified pipeline workflow that will validate existing primers *in silico*, define the identity of their amplicons (either by homology or genome coordinates) and produce a user-friendly html report (and a more machine-readable .txt report).

## Features
Key features of DiPPER2 (&Dagger; only in full pipeline):

* generates a html and .txt report with the primers, the identity of the region the primers target and information of specificity and sensitivity.
* generates .txt files (in fasta format) with primers for the unique genomic regions found&Dagger;
* generates fasta files with the sequence of the unique genomic regions the primers target&Dagger;
* generates .txt files that contain information on the primer's characteristics (Tm, GC, product size *etc*.)&Dagger;
* generates .bed files of unique genomic regions that cannot be identified by homology (blastx). These can be used in IGV to find where the unique genomic region falls in a reference genome (either provided, or longest)
* generates .bed files for the *in silico* tests of primers

For more details on the outputs and the folder structure, compare the 'DiPPER2 results walkthrough' section.

## Full DiPPER2 pipeline workflow

This workflow will find unique genomic regions, 
Each script initializes a logger from the custom Logger class, defined in ```logging_handler.py```.

1. ```DiPPER2_wrapper.sh``` or ```DiPPER2_wrapper.py``` checks if mandatory arguments have been provided, checks which python version is found, finds the default python interpreter, checks if dipper2 environment exists and if not, creates one (upgrades pip, installs requirements, activates the python environment). It then changes directory to the directory that contains all assemblies. 
2. Called by ```DiPPER2_wrapper.sh/py```, ```target_move_module_optimized.py``` is run. This module creates two folders, FUR.target and FUR.neighbour, and sorts the target and neighbour assemblies from the current directory into the respective folders.
3. If the previous module successfully finished, the third party software ```FUR``` is run by calling ```FUR_module_optimized```. It checks if all folder are present and non-empty, then makes a ```FUR``` database and runs ```FUR```. ```FUR``` provides some statistics about its run, those are saved in ```FUR_Summary_output.txt```. Finally, the module deletes the ```FUR``` database, as it takes up a lot of space.
4. Once that has run successfully, ```DiPPER2_wrapper.sh/py``` calls ```Primer3_module_optimized.py```, which finds the ```FUR``` output file. ```FUR```s output, found in <outfile_prefix>_FUR.db.out.txt, does not have the right format for ```Primer3```, which will be run next. Therefore, the script imports functions to convert the output from ```fur2primer3.py```. It checks if it has all the primer3 parameters necessary (```-t``` option) and reformats the ```FUR``` output as a ```Primer3``` input. The resulting file is <outfile_prefix>_FUR.db.out.txt.primers.txt. Then ```primer3_core``` is run. Results from that run are in <outfile_prefix>_FUR.db.out.txt.primers.txt.primer3_out.txt. From there, it extracts the four [lowest penalty](https://www.primer3plus.com/primer3plusHelp.html#calculatePenalties) primers, puts them into primer files and also extracts their targets, as well as key primer parameters (Tm etc.) from ```Primer 3```'s output. Primers are identified by number which corresponds to the line number of the primer found in the Primer3 output file, and the files with the primers are moved into FUR.P3.PRIMERS. Target fastas are moved into FUR.P3.TARGETS. FUR.P3.PRIMERS has a subfolder primer_data which contains the data for each primer.
5. If this has successfully finished, the wrapper calls ```Primer_Testing_module_optimized.py```. This checks if it has all it needs (primers, targets etc), then concatenates all assemblies in the FUR.target folder and all assemblies in the FUR.neighbour folder in two files, respectively. It then runs first sensitivity tests using ```seqkit amplicon``` for *in silico* PCR on the concatenated targets. It does that iteratively up till 3 mismatches allowed. Then it runs specificity tests on the neighbours, iteratively for all four primers with up to 4 mismatches. This test can time out, the time out is set to 6 minutes, as often ```seqkit amplicon``` find huge (several kb) off-target amplicons when primer pairs randomly match in the genome of neighbours. After finishing that, it runs ```blastx``` on the targets. If ```blastx``` does not return any hits, it either takes the reference provided by the user (```-r``` option) or finds the longest assembly within the targets and runs ```seqkit locate``` with the primer target. The result is a bedfile with coordinates where the primer target is found in the assembly. Finally, it puts all *in silico* testing files into a separate subfolder of FUR.P3.PRIMERS call in_silico_tests.
6. Upon successful finish of this script, the ```Summarize_results_module_improved.py``` module is run. It will produce an html report for the user that will tell them whether the sets of primers failed or passed. Sensitivity tests are passed when the primer pair produces an amplicon of the correct size for all target assemblies. Off-target amplification or missing target assemblies is accepted for *in silico* tests with 3 mismatches. For the specifictity tests, the script checks whether no neighbours produced an amplicon with the correct size. No correct sized amplicons of any neighbour are accepted for up to 4 mismatches. The report also contains information on the targets, and lists which files are found where.

## Requirements
The following programs need to be in path for __the full DiPPER2 pipeline__:

* [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
* [seqkit](https://bioinf.shenwei.me/seqkit/download/)
* [FUR](https://github.com/EvolBioInf/fur)
* [Primer3](https://primer3.org/releases.html)
* [E-utilities](https://www.ncbi.nlm.nih.gov/home/tools/)

For __primer testing only standalone-pipeline__, ```FUR``` and ```Primer3``` are *not* needed. 

DiPPER2 was developed using python 3.9 and 3.12. It runs stably on both version 3.9 and 3.12, but a python version >=3.12 is recommended.

DiPPER2 automatically generates a python environment and installs the necessary python modules defined in the requirements file.

Seqkit, BLAST+ and primer3 can be installed via conda. primer3 can also be installed using homebrew. 


## Installation

To install the dependencies, please run the ```dependency_install.sh``` script. For both Linux and MacOS you can decide whether you want to install the entire pipeline (including ```FUR``` and ```Primer3```) or the standalone, primer-testing one.
By default it will install the full pipeline (```-s n```). To install the standalone-pipeline, set ```-s y```.

```Bash
Usage: dependency_install.sh -m n -s n
    Installation Script. Either installs the entire DiPPER2 pipeline dependencies or just dependencies of the standalone pipeline, which tests primers only.
      -m <Is this a MacOS?> [y/n; default n]
      -s <standalone or not?> [y/n; default n]
```
To install ```FUR``` successfully, this installation script might have to be run as sudo.

__NOTE: the MacOS install of ```FUR``` will fail!__</br>
Follow the following steps:
1. check which of the following might already be installed:
    - golang 
    - libdivsufsort-dev 
    - make
    - phylonium
    - homebrew
2. Install homebrew if necessary first:
```Bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)
```
3. Install all missing other packages using
```Bash
brew install <package>
```

Assuming all the required programs listed in Requirements are then installed and in path, DiPPER2 can be installed by cloning the repo:
```Bash
git clone https://github.com/ThWacker/DiPPER2.git
```
> *DiPPER2 should now be ready to go*

## Usage

> __*PLEASE NOTE THAT CURRENTLY RELATIVE PATHS ARE NOT RESOLVED PROPERLY. USE ABSOLUTE PATHS. You might have to run ```conda activate``` first for all programs to be in PATH.*__

### Minimal usage (non-parallel):

```Bash
/repos/DiPPER2/scripts/DiPPER2_wrapper.sh -d <folder with the assemblies> -f <name of the results folder> -o <results files prefix> -q <toggle if qPCR primers are wanted or not, default n> -l <list of targets>
```

#### Optional parameters:

```
-c <delete concatenated files. Default:1. Set to 0 if you don't want that> 
-p <FUR parameters>
-t <primer3 parameters> [default: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200 Oligo=1]
-q <qpcr (y) or conventional pcr (n)> [default: n] 
-r <reference for bed files>
-a <assembly used as reference for FUR>
```

### Minimal usage (parallel):
```Bash
./DiPPER_wrapper_parallel.sh -f <results folder> -d <folder with the assemblies> -l <list with targets> 
```

#### Optional parameters:
```
-o <outfile prefix> [default: date]
-p <FUR parameters> [default: none]
-t <primer3 parameters> [default: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200 Oligo=1] 
-q <qpcr (y) or conventional pcr (n)> [default: n]
-r <reference for bed files> [default: longest target assembly]
-a <assembly used as reference for FUR> [default: longest target assembly]
-m <max memory that a job can use before getting killed> [default: 16000000000]
-n <min memory that needs to be available to run a job> [default: 5000000000]
-c <check interval for when not enough memory is available to start a job. After <-c> seconds, another memory check is conducted> [default: 2]
-w <max number of workers in the pool of workers for jobs> [default: 6]
```


## Additional information on the parallelized version of DiPPER2

The parallelized version of DiPPER2 does not, like the unparallelized version, concatenate all targets and neighbour assembly fastas and run the *in silico* PCR on those, but uses the ```multiprocessing``` package to run ```seqkit amplicon ```, the tool which is used for *in-silico* PCR, on individual assemblies in parallel. To make sure this does not overwhelm the computer or server, the ```Primer_testing_parallelize.py``` script will check if there is sufficient memory available (currently set to 5 Gb) before starting each job. If no sufficient memory is available, it will wait 2 seconds and then check again. It does that up to 6h long before timing out completely (the script will then stop). During running the multiple processes, a Daemon thread will monitor the memory usage of each process and if it exceeds 16Gb, kills it. 

__*Notes*__:

__Multiprocessing__ uses "spawn" as start method to avoid file descriptors being passed on and to make it save for MacOS. However, this can lead to other problems, especially due to the following documented restrictions:

>__More picklability:__ Ensure that all arguments to Process.__init__() are picklable. Also, if you subclass Process then make sure that instances will be picklable when the Process.start method is called. </br>
>__Global variable:__ Bear in mind that if code run in a child process tries to access a global variable, then the value it sees (if any) may not be the same as the value in the parent process at the time that Process.start was called.
>However, global variables which are just module level constants cause no problems.

If that turns out to be a problem, change the ```set_start_method("spawn")``` in the main function to ```fork``` or comment that line out in ```Primer_testing_parallelize.py```. That might solve the problem. 

## Full DiPPER2 pipeline run results walkthrough

### The results folder contains:


- __subfolders:__
    - FUR.P3.PRIMERS
    - FUR.P3.PRIMERS/primer_data
    - FUR.P3.PRIMERS/in_silico_tests
    - FUR.P3.TARGETS

- __files:__
    - Results.txt and Results.html: </br> the results summary files
    - `<Prefix>`.targets.txt and `<Prefix>`.neighbours.txt: files listing the accessions/ files used as targets and neighbours
    - `<Prefix>`.FUR.db.out.txt: </br>  the FUR results file (a multifasta file)
    - <Prefix>.FUR.db.out.primers.txt: </br>  the FUR results file reformatted in the format required by Primer3 as an input
    - `<Prefix>`.FUR.db.out.primers.primer3_out.txt: </br>  the primer3 output results file
    - FUR_Summary_output.txt: </br> FUR output to STDERR which summarizes how many sequences were retained after each processing step, including the length in bp and the number of Ns
    - *optional:*</br>  bed files generated by seqkit locate


#### Within the subfolders, the following files can be found:


- __within FUR.P3.PRIMERS:__ </br>
    Primer fasta files (.txt ending, but fasta formatted). Primers are numbered. The number is the unique identity of the primer.
- __within FUR.P3.PRIMERS/primer_data:__</br>
    Here, textfiles with information about Tm, amplicon length, GC etc are found. These follow the Primer3 conventions, please compare here for explanation: ([https://primer3.org/manual.html#outputTags](https://primer3.org/manual.html#outputTags))
- __within FUR.P3.PRIMERS/in_silico_tests__ </br>
    Here, the results of the in silico PCR tests, performed using seqkit locate, can be found. Files names with 'target' in the name are files used for sensitivity testing, running an in silico PCR against the concatenated targets.
    File names with 'neighbour' in the name are files used for specificity testing, running an in silico PCR against the concatenated neighbours. The 'm' in the name refers to the number of allowed mismatches in the primer
- __within FUR.P3.TARGETS:__ </br>
    Here, we find target files with 'Target' in the name, followed by a number that matches the unique identifier of the primer pair, which are fasta files containing the target sequence.
    Also, blastx results are found here. These are tab-separated and the headers of the blastx files are "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

### Other information:

All relevant information are found in the results files, which will give you an idea whether a primer passed or failed. Currently, the results files do not contain amplicon lengths or Tms etc, please refer to the txt files in the FUR.P3.PRIMERS/primer_data folder for that. 

Up to 4 primers are picked, but sometimes less than 4 primers are generated by Primer3. Having less than 4 primers in the results does not mean the run/ pipeline did not complete.

Currently, the default target parameters for conventional primers are:</br> 
```primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=200 prodMaxSize=1000 Oligo=0```

The default target parameters for qPCR primers are:</br> 
```primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200```

This means that conventional primers have a target optimum Tm of 60 and a product size of 200-1000 and qPCR primers a Tm of 60 and an amplicon size of 100-200. For qPCR primers an internal probe is also picked with an optimal Tm of 65.

## Future developments

* FUR.target and FUR.neighbour assemblies are listed in target.txt and neighbour.txt and the folders are deleted
* ~~A working wrapper for the parallelized primer testing~~ &check;
* ~~A standalone pipeline that tests user provided primers (currently pre-release version in branch parallelize)~~&check;
* optimized blastx strategy to account for Ns in FURs output
* packaging for PyPi
* Docker containerization and conda packaging
* extension of testing suite, refactoring to pytest

  
## Contact

Either open an issue in this repo or contact Theresa Wacker at theresa[underscore]wacker[at]mailbox[dot]org.

## Acknowledgements

DiPPER2 uses FUR, primer3, seqkit, BLAST+ and eutilies, all of which have been published. When using DiPPER2, please also cite:

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K. and Madden, T.L. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10(1), p.421.

Haubold, B., Klötzl, F., Hellberg, L., Thompson, D. and Cavalar, M. (2021). Fur: Find unique genomic regions for diagnostic PCR. Bioinformatics, 37(15), pp.2081–2087. doi:https://doi.org/10.1093/bioinformatics/btab059.

Kans, J. (2016). Entrez Direct: E-utilities on the UNIX Command Line. [online] Available at: https://www.ncbi.nlm.nih.gov/books/NBK179288/ [Accessed 4 Dec. 2024].

Shen, W., Le, S., Li, Y. and Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), p.e0163962. doi:https://doi.org/10.1371/journal.pone.0163962.

Untergasser, A., Cutcutache, I., Koressaar, T., Ye, J., Faircloth, B.C., Remm, M. and Rozen, S.G. (2012). Primer3—new capabilities and interfaces. Nucleic Acids Research, 40(15), pp.e115–e115. doi:https://doi.org/10.1093/nar/gks596.


