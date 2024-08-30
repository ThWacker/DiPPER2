#!/usr/bin/env bash 
set -e

# Complaints, bugs, praise, requests, cute cat pics, complaints etc can be sent to theresa_wacker@mailbox.org or t.wacker2@exeter.ac.uk

die() { echo "$@" ; exit 1; }
diemsg() {
    echo "Usage: $0 -f <results folder> -d <folder with the assemblies> -l <list with targets> -o <outfile prefix -c <delete concatenated files. Default:1. Set to 0 if you don't want that> 
    -p <FUR parameters>  -t <primer3 parameters> [default: primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200] -q <qpcr (y) or conventional pcr (n)> [default: n] 
    -r <reference for bed files> -a <assembly used as reference for FUR>"  
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
REF_FUR=

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
    -a) REF_FUR="$2"; shift;;
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

#check if environment exists or create one

if [ ! -d dipper2 ]; then
    "$default_python" -m venv dipper2
fi

dipper2/bin/pip install --upgrade pip
dipper2/bin/pip install -r "$SCRIPT_DIR"/requirements.txt
source dipper2/bin/activate
# Shift into assembly folder
cd "$ASSEM_F" || exit 1

# Shift all the files into target and neighbour folders
echo "Shifting all the target assemblies into FUR.target and all the neighbours into FUR.neighbour folders"
"$default_python" ~/repos/repos/python_scripts/target_move_module_optimized.py -t "$TARGET" -f "$FOLD"

# Run FUR
echo "Running FUR"
 
if [[ -n $OUT && -n $FUR && -n $REF_FUR ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -p "$FUR" -o "$OUT" -r "$REF_FUR"&& true
    EXIT_STATUS="$?"
elif [[ -n $FUR && -n $REF_FUR ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -p "$FUR" -r "$REF_FUR"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT && -n $REF_FUR ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -o "$OUT" -r "$REF_FUR"&& true
    EXIT_STATUS="$?"
elif [[ -n $FUR ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -p "$FUR"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT ]]; then
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD" -o "$OUT"&& true
    EXIT_STATUS="$?"
else
    "$default_python"  "$SCRIPT_DIR"/FUR_module_optimized.py -f "$FOLD"&& true
    EXIT_STATUS="$?"
fi

#check if the command failed with an exit status other than 0
if [[ $EXIT_STATUS -ne 0 ]]; then
    echo "FUR failed with error ${EXIT_STATUS}. Check log for details!"
    exit 1
fi
echo "Finished running FUR"

# Run P3 picking
echo "Pick primers using Primer3"
convPCR="primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=200 prodMaxSize=1000 Oligo=0"
if [[ -n $OUT && -n $P3 && $QPCR == "y" ]]; then
    "$default_python" "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" -p "$P3"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT && -n $P3 ]]; then
    "$default_python" "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" -p "$P3"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT && $QPCR == "n" ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" -p "$convPCR"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT && $QPCR == "y" ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" && true
    EXIT_STATUS="$?"
elif [[ -n $OUT ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -o "$OUT" -p "$convPCR"&& true
    EXIT_STATUS="$?"
elif [[ -n $P3 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -p "$P3"&& true
    EXIT_STATUS="$?"
elif [[ $QPCR == "y" ]]; then
    "$default_python" "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" && true
    EXIT_STATUS="$?"
else
    "$default_python" "$SCRIPT_DIR"/Primer3_module_optimized.py -f "$FOLD" -p "$convPCR"&& true
    EXIT_STATUS="$?"
fi

#check if the command failed with an exit status other than 0
if [[ $EXIT_STATUS -ne 0 ]]; then
    echo "Primer picking failed with error ${EXIT_STATUS}. Check log for details!"
    exit 1
fi
echo "Finished primer picking."

# Run in silico tests and blastx
echo "Testing the primers for specificity and sensitivity in silico and determining the target"

if [[ -n $OUT && $DEL -eq 0 && -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD" -o "$OUT" -c "$DEL" -r "$REF"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT && -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD" -o "$OUT" -r "$REF"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT && $DEL -eq 0 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD" -o "$OUT" -c "$DEL"&& true
    EXIT_STATUS="$?"
elif [[ $DEL -eq 0 ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"  -c "$DEL"&& true
    EXIT_STATUS="$?"
elif [[ -n $OUT ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"  -o "$OUT"&& true
    EXIT_STATUS="$?"
elif [[ -n $REF ]]; then
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"  -r "$REF"&& true
    EXIT_STATUS="$?"
else
    "$default_python"  "$SCRIPT_DIR"/Primer_Testing_module_optimized.py -f "$FOLD"&& true
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

if  [[ $QPCR == "y" ]]; then
    "$default_python"  "$SCRIPT_DIR"/Summarize_results_module_improved.py -f "$FOLD" -q "$QPCR"&& true
    EXIT_STATUS="$?"
else
    "$default_python"  "$SCRIPT_DIR"/Summarize_results_module_improved.py  -f "$FOLD"&& true
    EXIT_STATUS="$?"
fi

#check if the command failed with an exit status other than 0
if [[ $EXIT_STATUS -ne 0 ]]; then
    echo "Could not generate Results files. Summarize_results_module failed with error ${EXIT_STATUS}. Check log for details!"
    exit 1
fi
deactivate
echo "Finished!"