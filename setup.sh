#!/bin/bash

# ==============================================================================
# setup.sh 
# downloads and prepares the static databases required for the pipeline. 
# ==============================================================================

# strict mode: fail on errors, unassigned variables, and pipe failures
set -euo pipefail

# --- CONFIGURATION ------------------------------------------------------------
# (!!!) the links below may become obsolete, which would cause errors. 
# if this happens, feel free to create and issue in the repo. thank you!
KRAKEN_DIR="data/dbs/kraken2"
KRAKEN_DB_VERSION="k2_pluspf_08gb_20240112.tar.gz"
KRAKEN_URL="https://genome-idx.s3.amazonaws.com/kraken/${KRAKEN_DB_VERSION}"

# for when annotation is implemented:
# EGGNOG_DIR="data/eggnog_db"
# DBCAN_DIR="data/dbcan_db"
# ------------------------------------------------------------------------------

# helper functions for formatted output
print_info() { echo -e "\n[INFO] $1"; }
print_error() { echo -e "\n[ERROR] $1" >&2; }
print_success() { echo -e "\n[SUCCESS] $1"; }

cat << 'EOF'

=======================================================================


    o O       o O       o O       o O       o O       o O       o O
  o | | O   o | | O   o | | O   o | | O   o | | O   o | | O   o | | O 
O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O
  O | | o   O | | o   O | | o   O | | o   O | | o   O | | o   O | | o 
    O o       O o       O o       O o       O o       O o       O o   
 __  __                __  ___                             __    __     
 \ \/ /__  ____ ______/ /_/   |  _____________  ____ ___  / /_  / /__   
  \  / _ \/ __ `/ ___/ __/ /| | / ___/ ___/ _ \/ __ `__ \/ __ \/ / _ \  
  / /  __/ /_/ (__  ) /_/ ___ |(__  |__  )  __/ / / / / / /_/ / /  __/  
 /_/\___/\__,_/____/\__/_/  |_/____/____/\___/_/ /_/ /_/_.___/_/\___/   

    o O       o O       o O       o O       o O       o O       o O
  o | | O   o | | O   o | | O   o | | O   o | | O   o | | O   o | | O 
O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O
  O | | o   O | | o   O | | o   O | | o   O | | o   O | | o   O | | o 
    O o       O o       O o       O o       O o       O o       O o   


=======================================================================
                              SETTING UP...
=======================================================================
EOF


# ==============================================================================
# KRAKEN2 DATABASE SETUP
# ==============================================================================
if [ ! -f "${KRAKEN_DIR}/taxo.k2d" ]; then
    print_info "Kraken2 database not found. Starting download (~6GB)..."
    
    # create directory and navigate into it
    if ! mkdir -p "${KRAKEN_DIR}"; then
        print_error "Could not create directory: ${KRAKEN_DIR}"
        print_error "Suggestion: Check your write permissions in this folder."
        exit 1
    fi
    
    cd "${KRAKEN_DIR}"
    
    # download with error catching
    print_info "Downloading ${KRAKEN_DB_VERSION}..."
    if ! wget --show-progress "${KRAKEN_URL}"; then
        print_error "Failed to download the Kraken2 database."
        print_error "Suggestion: The link may be broken, or you might have network issues."
        print_error "Please verify if this URL is still valid: ${KRAKEN_URL}"
        exit 1
    fi
    
    # extraction with error catching
    print_info "Extracting the database..."
    if ! tar -xzf "${KRAKEN_DB_VERSION}"; then
        print_error "Failed to extract the Kraken2 database."
        print_error "Suggestion: The downloaded archive might be corrupted or incomplete."
        print_error "Try deleting the file '${KRAKEN_DIR}/${KRAKEN_DB_VERSION}' and running setup.sh again."
        exit 1
    fi
    
    # cleanup
    print_info "Cleaning up temporary files..."
    rm "${KRAKEN_DB_VERSION}"
    cd ../../..
    
    print_success "Kraken2 database configured successfully!"
else
    print_info "Kraken2 database already exists. Skipping download."
fi

# ==============================================================================
# MINIMAP2 HUMAN REFERENCE GENOME SETUP (for decontamination)
# ==============================================================================
HUMAN_DIR="data/dbs/human"
HUMAN_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
HUMAN_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/${HUMAN_FILE}"

if [ ! -f "${HUMAN_DIR}/${HUMAN_FILE}" ]; then
    print_info "Human reference genome not found. Starting download (~800MB)..."
    
    if ! mkdir -p "${HUMAN_DIR}"; then
        print_error "Could not create directory: ${HUMAN_DIR}"
        exit 1
    fi
    
    cd "${HUMAN_DIR}"
    
    if ! wget --show-progress "${HUMAN_URL}"; then
        print_error "Failed to download the Human reference genome."
        exit 1
    fi
    
    cd ../../..
    print_success "Human reference genome configured successfully!"
else
    print_info "Human reference genome already exists. Skipping download."
fi


# ==============================================================================
# BOWTIE2 INDEX SETUP (alternative for decontamination)
# ==============================================================================
BOWTIE2_DIR="data/dbs/human/bowtie2"
BOWTIE2_URL="https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip"
BOWTIE2_FILE="GRCh38_noalt_as.zip"

if [ ! -f "${BOWTIE2_DIR}/GRCh38_noalt_as.1.bt2" ]; then
    print_info "Bowtie2 index not found. Starting download (~4GB)..."
    
    if ! mkdir -p "${BOWTIE2_DIR}"; then
        print_error "Could not create directory: ${BOWTIE2_DIR}"
        exit 1
    fi
    
    # download with error catching
    print_info "Downloading Bowtie2 index archive..."
    if ! wget --show-progress "${BOWTIE2_URL}" -O "${BOWTIE2_DIR}/${BOWTIE2_FILE}"; then
        print_error "Failed to download the Bowtie2 index."
        exit 1
    fi
    
    # extraction with error catching
    print_info "Extracting the index..."
    if ! unzip -j "${BOWTIE2_DIR}/${BOWTIE2_FILE}" -d "${BOWTIE2_DIR}/"; then
        print_error "Failed to extract the Bowtie2 index."
        print_error "Suggestion: The downloaded zip archive might be corrupted."
        exit 1
    fi
    
    # cleanup
    print_info "Cleaning up temporary files..."
    rm "${BOWTIE2_DIR}/${BOWTIE2_FILE}"
    
    print_success "Bowtie2 index configured successfully!"
else
    print_info "Bowtie2 index already exists. Skipping download."
fi


# for when annotation is implemented:
# ==============================================================================
# EGGNOG-MAPPER DATABASE SETUP
# ==============================================================================

# ==============================================================================
# dbCAN3 DATABASE SETUP 
# ==============================================================================


echo "====================================================================="
echo " Setup Completed! The environment is ready."
echo " Now you might want to run: snakemake assemble --cores all"
echo "====================================================================="