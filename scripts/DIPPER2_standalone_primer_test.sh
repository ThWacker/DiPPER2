#!/usr/bin/env bash 
set -e

# Complaints, bugs, praise, requests, cute cat pics, complaints etc can be sent to theresa_wacker@mailbox.org or t.wacker2@exeter.ac.uk

die() { echo "$@" ; exit 1; }
diemsg() {
    echo "Usage: $0 -f <results folder> -d <folder with the assemblies> -l <list with targets> -F <forward primer> -R <reverse primer> 
    -o <outfile prefix 
    -c <delete concatenated files. Default:1. Set to 0 if you don't want that> 
    -r <reference for bed files> 
    -v <version>"  
    echo ""
    echo "Arguments -f, -F, -R, -d, and -l are mandatory."
    echo "BLAST+, seqkit, python and all dependencies must be installed in path!!!"
    echo "Primers must be strings without whitespaces. They can be quoted or unquoted."
    echo ""
    echo "-h for help"
    die ""
}

echo "Welcome! This is the DiPPER2 primer testing wrapper. This pipeline will test your primers in silico for specificity and sensitivity and identify the target using blastx"
echo ""
echo "-h for help"
echo ""

if [ "$1" == "-h" ] || [ $# -lt 5 ]; then
    diemsg
fi

# Variable assignment
FOLD=
ASSEM_F=
TARGET=
OUT=
FRW=
REV=
DEL=1
REF=
VERSION="1.1.0"

#Display version
if [ "$1" == "-v" ]; then
    echo "DiPPER2 version ${VERSION}"
fi


# Initialize variables to track whether mandatory options are provided
a_provided=false
b_provided=false
c_provided=false
d_provided=false
e_provided=false

# Gather input files and settings
while [ $# -gt 0 ]; do
    case "$1" in
    -f) FOLD="$2"; shift; a_provided=true;;
    -d) ASSEM_F="$2"; shift; b_provided=true;;
    -l) TARGET="$2"; shift; c_provided=true;;
    -o) OUT="$2"; shift;;
    -F) FRW="$2"; shift; d_provided=true;;
    -R) REV="$2"; shift; e_provided=true;;
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
if [ "$a_provided" = false ] || [ "$b_provided" = false ] || [ "$c_provided" = false ] || [ "$d_provided" = false ] || [ "$e_provided" = false ]; then
    echo "Error: The values -f, -d, -F, -R and -l are mandatory!"
    diemsg 
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

# get the parent of the script directory

PARENT_DIR="$(dirname "$SCRIPT_DIR")"
#check if environment exists or create one

if [ ! -d dipper2 ]; then
    "$default_python" -m venv dipper2
fi

dipper2/bin/pip install --upgrade pip
dipper2/bin/pip install -r "$PARENT_DIR"/requirements.txt
source dipper2/bin/activate

# Shift into assembly folder
cd "$ASSEM_F" || exit 1

# Shift all the files into target and neighbour folders
echo "Shifting all the target assemblies into FUR.target and all the neighbours into FUR.neighbour folders"
"$default_python" "$SCRIPT_DIR"/target_move_module_optimized.py -t "$TARGET" -f "$FOLD"


# Run in silico tests and blastx
echo "Testing the primers for specificity and sensitivity in silico and determining the target"

if [[ -n $OUT && $DEL -eq 0 && -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -F "$FRW" -R "$REV" -o "$OUT" -c "$DEL" -r "$REF" && true
    EXIT_STATUS="$?"
elif [[ -n $OUT && -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -F "$FRW" -R "$REV" -o "$OUT" -r "$REF"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT && $DEL -eq 0 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -o "$OUT" -c "$DEL" && true
    EXIT_STATUS="$?"
elif [[ -n $REF && $DEL -eq 0 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -F "$FRW" -R "$REV" -r "$REF" -c "$DEL" && true
    EXIT_STATUS="$?"
elif [[ $DEL -eq 0 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -F "$FRW" -R "$REV"  -c "$DEL"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -F "$FRW" -R "$REV"  -o "$OUT"&& true
    EXIT_STATUS="$?"
elif [[ -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -F "$FRW" -R "$REV"  -r "$REF"&& true
    EXIT_STATUS="$?"
else
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_standalone.py -f "$FOLD" -F "$FRW" -R "$REV" && true
    EXIT_STATUS="$?"
fi

#check if the command failed with an exit status other than 0
if [[ $EXIT_STATUS -ne 0 ]]; then
    echo "Primer testing failed with error ${EXIT_STATUS}. Check log for details!"
    exit 1
fi

echo "Finished in silico PCR and target definition."

# Summary
echo "Writing summary outfiles/"
"$default_python"  "$SCRIPT_DIR"/Summarize_results_module_standalone.py  -f "$FOLD"&& true
EXIT_STATUS="$?"


# Clean up FUR.target and FUR.neighbour after saving all the targets and neighbours into txt files
#check if the command failed with an exit status other than 0
if [[ $EXIT_STATUS -ne 0 ]]; then
    echo "Could not generate Results files. Summarize_results_module failed with error ${EXIT_STATUS}. Check log for details!"
    exit 1
fi
deactivate
echo "Finished!"