#!/bin/bash

# ==============================================================================
# Setup Script 
# Downloads and prepares the static databases required for the pipeline.
# ==============================================================================

# Strict mode: fail on errors, unassigned variables, and pipe failures
set -euo pipefail

# --- CONFIGURATION (Update links and versions here if they become obsolete) ---
KRAKEN_DIR="data/dbs/kraken2"
KRAKEN_DB_VERSION="k2_pluspf_08gb_20240112.tar.gz"
KRAKEN_URL="https://genome-idx.s3.amazonaws.com/kraken/${KRAKEN_DB_VERSION}"

# EGGNOG_DIR="data/eggnog_db"
# DBCAN_DIR="data/dbcan_db"
# ------------------------------------------------------------------------------

# Helper functions for formatted output
print_info() { echo -e "\n[INFO] $1"; }
print_error() { echo -e "\n[ERROR] $1" >&2; }
print_success() { echo -e "\n[SUCCESS] $1"; }

cat << 'EOF'

O       o O       o O       o O       o O       o O       o O       o 
| O   o | | O   o | | O   o | | O   o | | O   o | | O   o | | O   o |
| | O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O | |
| o   O | | o   O | | o   O | | o   O | | o   O | | o   O | | o   O |
o       O o       O o       O o       O o       O o       O o       O
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


=====================================================================
 SETTING UP...
=====================================================================
EOF



# ==============================================================================
# 1. KRAKEN 2 DATABASE SETUP
# ==============================================================================
if [ ! -f "${KRAKEN_DIR}/taxo.k2d" ]; then
    print_info "Kraken2 database not found. Starting download (approx. 6GB)..."
    
    # Create directory and navigate into it
    if ! mkdir -p "${KRAKEN_DIR}"; then
        print_error "Could not create directory: ${KRAKEN_DIR}"
        print_error "Suggestion: Check your write permissions in this folder."
        exit 1
    fi
    
    cd "${KRAKEN_DIR}"
    
    # Download with error catching
    print_info "Downloading ${KRAKEN_DB_VERSION}..."
    if ! wget --show-progress "${KRAKEN_URL}"; then
        print_error "Failed to download the Kraken2 database."
        print_error "Suggestion: The link may be broken, or you might have network issues."
        print_error "Please verify if this URL is still valid: ${KRAKEN_URL}"
        exit 1
    fi
    
    # Extraction with error catching
    print_info "Extracting the database..."
    if ! tar -xzf "${KRAKEN_DB_VERSION}"; then
        print_error "Failed to extract the Kraken2 database."
        print_error "Suggestion: The downloaded archive might be corrupted or incomplete."
        print_error "Try deleting the file '${KRAKEN_DIR}/${KRAKEN_DB_VERSION}' and running setup.sh again."
        exit 1
    fi
    
    # Cleanup
    print_info "Cleaning up temporary files..."
    rm "${KRAKEN_DB_VERSION}"
    cd ../../..
    
    print_success "Kraken2 database configured successfully!"
else
    print_info "Kraken2 database already exists. Skipping download."
fi

# ==============================================================================
# 2. HUMAN REFERENCE GENOME SETUP (For Decontamination)
# ==============================================================================
HUMAN_DIR="data/dbs/human"
HUMAN_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
HUMAN_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/${HUMAN_FILE}"

if [ ! -f "${HUMAN_DIR}/${HUMAN_FILE}" ]; then
    print_info "Human reference genome not found. Starting download (approx. 800MB)..."
    mkdir -p "${HUMAN_DIR}"
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
# 3. EGGNOG-MAPPER DATABASE SETUP (Placeholder for next phase)
# ==============================================================================
# print_info "Checking eggNOG-mapper databases..."
# ... (code will be added here) ...

# ==============================================================================
# 4. dbCAN3 DATABASE SETUP (Placeholder for next phase)
# ==============================================================================
# print_info "Checking dbCAN3 databases..."
# ... (code will be added here) ...





echo "====================================================================="
echo " Setup Completed! The environment is ready."
echo " Now you might want to run: snakemake assemble --cores all"
echo "====================================================================="