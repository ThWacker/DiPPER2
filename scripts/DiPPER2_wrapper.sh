#!/usr/bin/env bash 
set -e

# Complaints, bugs, praise, requests, cute cat pics, complaints etc can be sent to theresa_wacker@mailbox.org or t.wacker2@exeter.ac.uk

die() { echo "$@" ; exit 1; }
diemsg() {
    echo "Usage: $0 -f <results folder> -d <folder with the assemblies> -l <list with targets> -o <outfile prefix -c <delete concatenated files. Default:1. Set to 0 if you don't want that> -p <FUR parameters> -t <primer3 parameters> -q <qpcr default or conventional pcr default> -r <reference>"  
    echo ""
    echo "Arguments -f, -d, and -l are mandatory."
    echo "FUR (https://github.com/EvolBioInf/fur/tree/master), primer3_core (https://github.com/primer3-org/primer3), BLAST+, seqkit, python and all dependencies must be installed in path!!!"
    echo "For the FUR parameters (-p), the default is to run all three steps: subtraction using macle, intersection using phylonium and another subtraction using blastn. If you want only the first step, select -u, the first and second: -U and if you want to use megablast instead of blastn, use -m."
    echo ""
    echo "-h for help"
    die ""
}

echo "Welcome! This is the DiPPER2 wrapper. This pipeline will find unique genomic regions in the target assemblies, build primers for those, test them in silico for specificity and sensitivity and identify the target using blastx"
echo ""
echo "-h for help"
echo ""

if [ "$1" == "-h" ] || [ $# -lt 4 ]; then
    diemsg
fi

# Variable assignment
FOLD=
ASSEM_F=
TARGET=
OUT=
FUR=
P3=
QPCR="n"
DEL=1
REF=

# Initialize variables to track whether mandatory options are provided
a_provided=false
b_provided=false
c_provided=false

# Gather input files and settings
while [ $# -gt 0 ]; do
    case "$1" in
    -f) FOLD="$2"; shift; a_provided=true;;
    -d) ASSEM_F="$2"; shift; b_provided=true;;
    -l) TARGET="$2"; shift; c_provided=true;;
    -o) OUT="$2"; shift;;
    -p) FUR="$2"; shift;;
    -t) P3="$2"; shift;;
    -q) QPCR="$2"; shift;;
    -c) DEL="$2"; shift;;
    -r) REF="$2"; shift;;
    -*) echo >&2 "Unknown option: $1"; diemsg;;
    *) break;;    # terminate while loop
    esac
    shift
done

##############
### CHECKS ###
##############
echo "Checking if mandatory arguments have been provided..."

# CHECK1: Check if mandatory arguments are provided
if [ "$a_provided" = false ] || [ "$b_provided" = false ] || [ "$c_provided" = false ]; then
    diemsg "Error: The values -f, -d, and -l are mandatory!"
fi

echo "Mandatory arguments have been provided!"

# CHECK2: Function to check if a given Python command exists and print its version
echo "Checking python version..."

check_python_version() {
    local potential_paths=(
        "$1"
        "/Library/Frameworks/Python.framework/Versions/3.9/bin/$1"
        "/usr/local/bin/$1"
        "/usr/bin/$1"
    )

    for path in "${potential_paths[@]}"; do
        if command -v "$path" >/dev/null 2>&1; then
            "$path" --version 2>&1
            return 0
        fi
    done

    echo "Command $1 not found."
    return 1
}

# Capture the output of the check_python_version function
python_version=$(check_python_version python || true)
python3_version=$(check_python_version python3 || true)

# Determine the default Python interpreter
if [[ $python_version == *"Python "* ]]; then
    echo "python is available: $python_version"
    default_python="python"
elif [[ $python3_version == *"Python "* ]]; then
    echo "python3 is available: $python3_version"
    default_python="python3"
else
    diemsg "Neither python nor python3 is available."
fi
############
### MAIN ###
############
# Get the script directory:
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Shift into assembly folder
cd "$ASSEM_F" || exit 1

# Shift all the files into target and neighbour folders
echo "Shifting all the target assemblies into FUR.target and all the neighbours into FUR.neighbour folders"
"$default_python" ~/repos/repos/python_scripts/target_move_module_optimized.py -t "$TARGET" -f "$FOLD"

# Run FUR
echo "Running FUR"
if [[ -n $OUT && -n $FUR ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -p "$FUR" -o "$OUT"
elif [[ -n $FUR ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -p "$FUR"
elif [[ -n $OUT ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -o "$OUT"
else
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD"
fi

echo "Finished running FUR"

# Run P3 picking
echo "Pick primers using Primer3"
convPCR="prodMinSize=200 prodMaxSize=1000"
if [[ -n $OUT && -n $P3 && $QPCR == "y" ]]; then
    "$default_python" "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" -p "$P3"
elif [[ -n $OUT && -n $P3 ]]; then
    "$default_python" "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" -p "$P3"
elif [[ -n $OUT && $QPCR == "y" ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" -p "$convPCR"
elif [[ -n $OUT ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT"
elif [[ -n $P3 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -p "$P3"
elif [[ $QPCR == "y" ]]; then
    "$default_python" "$SCRIPT_DIR"/scripts/Primer3_module_optimized.py -f "$FOLD" -p "$convPCR"
else
    "$default_python" "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD"
fi

echo "Finished primer picking."

# Run in silico tests and blastx
echo "Testing the primers for specificity and sensitivity in silico and determining the target"

if [[ -n $OUT && $DEL -eq 0 && -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD" -o "$OUT" -c "$DEL" -r "$REF"
elif [[ -n $OUT && -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD" -o "$OUT" -r "$REF"
elif [[ -n $OUT && $DEL -eq 0 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD" -o "$OUT" -c "$DEL"
elif [[ $DEL -eq 0 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"  -c "$DEL"
elif [[ -n $OUT ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"  -o "$OUT"
elif [[ -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"  -r "$REF"
else
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"
fi

echo "Finished in silico PCR and target definition."

# Summary
echo "Writing summary outfiles/"

if  [[ $QPCR == "y" ]]; then
    "$default_python"  "$SCRIPT_DIR"/Summarize_results_module_improved.py -f "$FOLD" -q "$QPCR"
else
    "$default_python"  "$SCRIPT_DIR"/Summarize_results_module_improved.py  -f "$FOLD"
fi

echo "Finished!"