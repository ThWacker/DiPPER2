[![CodeQL Advanced](https://github.com/ThWacker/DiPPER2/actions/workflows/codeql.yml/badge.svg)](https://github.com/ThWacker/DiPPER2/actions/workflows/codeql.yml) [![Python Lint and Code Quality](https://github.com/ThWacker/DiPPER2/actions/workflows/pylint.yml/badge.svg)](https://github.com/ThWacker/DiPPER2/actions/workflows/pylint.yml)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

 <!----
  (![Static Badge](https://img.shields.io/badge/unittest_code_coverage-78%25-yellow))
  --->

![GitHub last commit](https://img.shields.io/github/last-commit/ThWacker/DiPPER2)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14699018.svg)](https://doi.org/10.5281/zenodo.14699018)
![](./logo/Dipper2_white.png) 
# DiPPER2
DiPPER2:  Diagnostic Primer Picking and Evaluation pipeline for Reliability and Reproducibility

>__+++++NOTE - IMPORTANT:+++++__
>
>__This is a work in progress.__
>
>__++++++++++++++++++++++++++__

## Synopsis
__*This pipeline and modules are meant to facilitate reliable and reproducible finding of diagnostic targets and to make picking primers for those targets as user-friendly as possible. The approach taken is a phylogeny-driven and clade-specific approach.*__ 

The pipeline, once functional, automatically finds unique genomic regions for a target species, pathovar or race with moderate user intervention, automatically picks primers for these regions that are amenable to both conventional PCR and quantitative PCR, and finally tests them for specificity and sensitivity.

## Overview
DiPPER2 is a tool to find unique genomic regions of target species, build primers (conv. or qPCR) for them, validate them *in silico*, define the identity of the unique genomic regions (either by homology or genome coordinates) and produce a user-friendly html report (and a more machine-readable .txt report).

## Features
Key features of DiPPER2:

* generates a html and .txt report with the primers, the identity of the region the primers target and information of specificity and sensitivity.
* generates .txt files (in fasta format) with primers for the unique genomic regions found
* generates fasta files with the sequence of the unique genomic regions the primers target
* generates .txt files that contain information on the primer's characteristics (Tm, GC, product size *etc*.)
* generates .bed files of unique genomic regions that cannot be identified by homology (blastx). These can be used in IGV to find where the unique genomic region falls in a reference genome (either provided, or longest)
* generates .bed files for the *in silico* tests of primers

For more details on the outputs and the folder structure, compare the 'DiPPER2 results walkthrough' section.

## Requirements
The following programs need to be in path:

* [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
* [seqkit](https://bioinf.shenwei.me/seqkit/download/)
* [FUR](https://github.com/EvolBioInf/fur)
* [primer3_core](https://primer3.org/releases.html)
* [E-utilities](https://www.ncbi.nlm.nih.gov/home/tools/)

DiPPER2 was developed using python 3.9 and 3.12. It runs stably on both version 3.9 and 3.12, but a python version >=3.12 is recommended.

DiPPER2 automatically generates a python environment and installs the necessary python modules defined in the requirements file.

Seqkit, BLAST+ and primer3 can be installed via conda. primer3 can also be installed using homebrew. 


## Installation

Assuming all the required programs listed in Requirements are installed and in path, DiPPER2 can be installed by cloning the repo:
```Bash
git clone https://github.com/ThWacker/DiPPER2.git
```

## Usage

> __*PLEASE NOTE THAT CURRENTLY RELATIVE PATHS ARE NOT RESOLVED PROPERLY. USE ABSOLUTE PATHS.*__

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
> The wrapper for this is not yet available but will be made available shortly.

## Additional information on the parallelized version of DiPPER2

The parallelized version of DiPPER2 does not, like the unparallelized version, concatenate all targets and neighbour assembly fastas and run the *in silico* PCR on those, but uses the ```multiprocessing``` package to run ```seqkit amplicon ```, the tool which is used for *in-silico* PCR, on individual assemblies in parallel. To make sure this does not overwhelm the computer or server, the ```Primer_testing_parallelize.py``` script will check if there is sufficient memory available (currently set to 5 Gb) before starting each job. If no sufficient memory is available, it will wait 2 seconds and then check again. It does that up to 6h long before timing out completely (the script will then stop). During running the multiple processes, a Daemon thread will monitor the memory usage of each process and if it exceeds 16Gb, kills it. 

__Note__:

* ```-c 0``` is incompatible with the parallelized version of DiPPER2,as no concatenated files are generated in the parallelized version. 
* Currently, the __amount of workers/ cpus__ is set to __6__. This can be changed by going into the script ```Primer_testing_parallelize.py``` and changing the constant ```MAX_WORKERS``` in line 24.
* __The ```MAX_MEMORY_MB_PER_JOB```__ is not, as its variable name tries to imply, in Mb, but in bytes and is currently set to __16GB__. Again, this can be changed in the script ```Primer_testing_parallelize.py``` and changing the constant ```MAX_MEMORY_MB_PER_JOB``` in line 25.
* Finally, the __minimum memory required to start a parallel job__ is set to __5Gb__. This can be changed in the script ```Primer_testing_parallelize.py``` and changing the constant ```MIN_AVAILABLE_MEMORY``` in line 26.

## DiPPER2 results walkthrough

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

## Contact

Either open an issue in this repo or contact Theresa Wacker at t[dot]wacker2[at]exeter[dot]ac[dot]uk.

## Acknowledgements

DiPPER2 uses FUR, primer3, seqkit, BLAST+ and eutilies, all of which have been published. When using DiPPER2, please also cite:

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K. and Madden, T.L. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10(1), p.421.

Haubold, B., Klötzl, F., Hellberg, L., Thompson, D. and Cavalar, M. (2021). Fur: Find unique genomic regions for diagnostic PCR. Bioinformatics, 37(15), pp.2081–2087. doi:https://doi.org/10.1093/bioinformatics/btab059.

Kans, J. (2016). Entrez Direct: E-utilities on the UNIX Command Line. [online] Available at: https://www.ncbi.nlm.nih.gov/books/NBK179288/ [Accessed 4 Dec. 2024].

Shen, W., Le, S., Li, Y. and Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), p.e0163962. doi:https://doi.org/10.1371/journal.pone.0163962.

Untergasser, A., Cutcutache, I., Koressaar, T., Ye, J., Faircloth, B.C., Remm, M. and Rozen, S.G. (2012). Primer3—new capabilities and interfaces. Nucleic Acids Research, 40(15), pp.e115–e115. doi:https://doi.org/10.1093/nar/gks596.


