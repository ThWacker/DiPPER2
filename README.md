![](./logo/Dipper2_white.png)
# DiPPER2
DiPPER2:  Diagnostic Primer Picking and Evaluation pipeline for Reliability and Reproducibility

>__+++++++++++++++++++++++++++++NOTE - IMPORTANT:+++++++++++++++++++++++++++++++__
>
>__This is a work in progress which has not been extensively unit tested and will likely undergo substantial refactoring in the near future.__
>
>__+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++__

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

For more details on the outputs and the folder structure, compare RESULTS_MANUAL.md

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

## Contact

Either open an issue in this repo or contact Theresa Wacker at t[dot]wacker2[at]exeter[dot]ac[dot]uk.

## Acknowledgements

DiPPER2 uses FUR, primer3, seqkit, BLAST+ and eutilies, all of which have been published. When using DiPPER2, please also cite:

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K. and Madden, T.L. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10(1), p.421.

Haubold, B., Klötzl, F., Hellberg, L., Thompson, D. and Cavalar, M. (2021). Fur: Find unique genomic regions for diagnostic PCR. Bioinformatics, 37(15), pp.2081–2087. doi:https://doi.org/10.1093/bioinformatics/btab059.

Kans, J. (2016). Entrez Direct: E-utilities on the UNIX Command Line. [online] Available at: https://www.ncbi.nlm.nih.gov/books/NBK179288/ [Accessed 4 Dec. 2024].

Shen, W., Le, S., Li, Y. and Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), p.e0163962. doi:https://doi.org/10.1371/journal.pone.0163962.

Untergasser, A., Cutcutache, I., Koressaar, T., Ye, J., Faircloth, B.C., Remm, M. and Rozen, S.G. (2012). Primer3—new capabilities and interfaces. Nucleic Acids Research, 40(15), pp.e115–e115. doi:https://doi.org/10.1093/nar/gks596.


