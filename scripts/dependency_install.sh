#!/usr/bin/env bash 
set -u

die() { echo "$@" ; exit 1; }
diemsg() {
    echo "Usage: $0 -m n -s n
    Installation Script. Either installs the entire DiPPER2 pipeline dependencies or just dependencies of the standalone pipeline, which tests primers only.
      -m <Is this a MacOS?> [y/n; default n]
      -s <standalone or not?> [y/n; default n]
    "  
    echo ""
    echo "-h for help"
    die ""
}

echo "Welcome! This is the DiPPER2 dependency installation wrapper."
echo ""
echo "-h for help"
echo ""

if [ "$1" == "-h" ] ; then
    diemsg
fi

# Variable assignment
ALONE="n"
MAC="n"
declare InstallArray

# Gather input files and settings
while [ $# -gt 0 ]; do
    case "$1" in
    -m) MAC="$2"; shift;;
    -s) ALONE="$2"; shift;;
    -*) echo >&2 "Unknown option: $1"; diemsg;;
    *) break;;    # terminate while loop
    esac
    shift
done

# Initialize logging
LOGFILE="$(pwd)/install.log"
echo "Detailed install log information" | tee -a "$LOGFILE"

# Utility: add to PATH
update_path() {
    local new_path="$1"
    echo "export PATH=$new_path:\${PATH}" >> "$HOME/.bashrc"
    echo "Added $new_path to PATH" | tee -a "$LOGFILE"
}

# install functions that return if they fail. 
install_primer3() {
    if [[ $MAC == "n" ]]; then
        echo "Installing Primer3..." | tee -a "$LOGFILE"
        git clone https://github.com/primer3-org/primer3.git || {
            echo "Primer3 install failed at git clone" | tee -a "$LOGFILE"; return; 
        }

        cd primer3/src || {
            echo "Primer3 directory not found" | tee -a "$LOGFILE"; return; 
        }

        make || { echo "Primer3 make failed" | tee -a "$LOGFILE"; return; }
        make test || { echo "Primer3 make test failed" | tee -a "$LOGFILE"; return; }

        update_path "$(pwd)"
        cd ../..
    else
        yes "" | mamba install bioconda::primer3
    fi
    echo "Primer3 installed successfully" | tee -a "$LOGFILE"
    InstallArray+=("Primer3")
}

install_fur() {
    echo "Installing FUR..." | tee -a "$LOGFILE"
    git clone https://github.com/evolbioinf/fur || { echo "FUR clone failed" | tee -a "$LOGFILE"; return; }
    cd fur || { echo "FUR directory missing" | tee -a "$LOGFILE"; return; }

    bash ./scripts/setup.sh || { echo "FUR setup.sh failed" | tee -a "$LOGFILE"; return; }
    make || { echo "FUR make failed" | tee -a "$LOGFILE"; return; }

    update_path "$(pwd)/bin"

    cd ..
    echo "FUR installed successfully" | tee -a "$LOGFILE"
    InstallArray+=("FUR")
}

install_edirect() {
    echo "Installing Entrez Direct E-utilities..." | tee -a "$LOGFILE"
    yes "" | sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" || {
        echo "Entrez Direct install failed" | tee -a "$LOGFILE"; return;
    }

    update_path "$HOME/edirect"

    echo "Entrez Direct installed successfully" | tee -a "$LOGFILE"
    InstallArray+=("Edirect")
}

install_blast() {
    echo "Installing ncbi-blast-2.2.31+..." | tee -a "$LOGFILE"
    if [[ $MAC == "n" ]]; then
        curl -s -f -Lo ncbi-blast-2.2.31+-x64-linux.tar.gz \
            https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz || {
            echo "Download of BLAST failed" | tee -a "$LOGFILE"; return;
        }

        tar -xzf ncbi-blast-2.2.31+-x64-linux.tar.gz || {
            echo "Extraction of BLAST failed" | tee -a "$LOGFILE"; return;
        }

        rm ncbi-blast-2.2.31+-x64-linux.tar.gz

        update_path "$(pwd)/ncbi-blast-2.2.31+/bin"
    else
        curl -s -f -Lo ncbi-blast-2.16.0+-x64-macosx.tar.gz  \
        https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-macosx.tar.gz || {
            echo "Download of BLAST failed" | tee -a "$LOGFILE"; return;
        }
        tar -xzf ncbi-blast-2.16.0+-x64-macosx.tar.gz || {
            echo "Extraction of BLAST failed" | tee -a "$LOGFILE"; return;
        }

        rm ncbi-blast-2.16.0+-x64-macosx.tar.gz

        update_path "$(pwd)/ncbi-blast-2.16.0+/bin"
    fi

    echo "BLAST+ installed successfully" | tee -a "$LOGFILE"
    InstallArray+=("BLAST+")
}

install_mamba() {
    echo "Installing mamba via Miniforge..." | tee -a "$LOGFILE"
    if [[ $MAC == "y" ]]; then
      curl -fsSLo Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh" || {
          echo "Failed to download Miniforge installer" | tee -a "$LOGFILE"; return;
      }
    else
      wget -O Miniforge3.sh \
          "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" || {
          echo "Failed to download Miniforge installer" | tee -a "$LOGFILE"; return;
      }
    fi 
    bash Miniforge3.sh -b -p "${HOME}/conda" || {
        echo "Miniforge installer failed" | tee -a "$LOGFILE"; return;
    }

    mamba shell init || {
        echo "Failed at mamba shell init" | tee -a "$LOGFILE"; return;
    }

    echo "Mamba installed successfully" | tee -a "$LOGFILE"
    InstallArray+=("Mamba")
}

install_seqkit() {
    echo "Installing seqkit..." | tee -a "$LOGFILE"
    if ! command -v mamba >/dev/null 2>&1; then
        install_mamba
    fi

    mamba install -y -c bioconda seqkit || {
        echo "Seqkit install via mamba failed" | tee -a "$LOGFILE"; return;
    }

    update_path "${HOME}/conda/bin"

    echo "Seqkit installed successfully" | tee -a "$LOGFILE"
    InstallArray+=("SeqKit")
}


# Test if programs are there and if not,install. Otherwise just notify that they are installed. 
if ! command -v efetch >/dev/null 2>&1; then install_edirect; else echo "Entrez Direct already installed" | tee -a "$LOGFILE"; fi
if ! command -v blastx >/dev/null 2>&1; then install_blast; else echo "BLAST+ already installed" | tee -a "$LOGFILE"; fi
if ! command -v mamba >/dev/null 2>&1; then install_mamba; else echo "Mamba already installed" | tee -a "$LOGFILE"; fi
if ! command -v seqkit >/dev/null 2>&1; then install_seqkit; else echo "Seqkit already installed" | tee -a "$LOGFILE"; fi
if [[ $ALONE == "n" ]]; then
  if ! command -v primer3_core >/dev/null 2>&1; then install_primer3; else echo "Primer3 already installed" | tee -a "$LOGFILE"; fi
  if ! command -v makeFurDb >/dev/null 2>&1; then install_fur; else echo "FUR already installed" | tee -a "$LOGFILE"; fi
fi

echo "All requested software installations attempted." | tee -a "$LOGFILE"

echo "The following software was installed:" | tee -a "$LOGFILE"
for str in "${InstallArray[@]}"; do
  echo "$str" | tee -a "$LOGFILE"
done


exit 0


