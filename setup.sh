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
# BIOCONTAINERS PRE-FETCH (READ FROM config.yaml)
# pre-builds SIF container images dynamically from the pipeline configuration
# ==============================================================================
print_info "Pre-fetching all BioContainers SIF images from config.yaml..."

mkdir -p data/containers/cache
export APPTAINER_CACHEDIR="$PWD/data/containers"

# parse container URLs from config.yaml using python
python -c "
import yaml
with open('config.yaml', 'r') as f:
    cfg = yaml.safe_load(f)
for tool, url in cfg.get('containers', {}).items():
    print(f'{tool}\t{url}')
" | while IFS=$'\t' read -r TOOL URL; do
    echo "[INFO] Pulling and verifying container for ${TOOL}: ${URL}"
    apptainer exec --cleanenv "${URL}" true
done

print_success "All pipeline containers pre-fetched and ready!"

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
# HOSTILE DATABASE SETUP (HUMAN DECONTAMINATION)
# Pre-fetches the masked human reference genome (T2T-CHM13v2.0 + HLA)
# ==============================================================================
echo "[INFO] Fetching Hostile human T2T-HLA reference index (Bowtie2)..."

# ensure target directories exist
mkdir -p data/dbs/hostile
mkdir -p data/containers/cache

# enforce local Apptainer cache directory
export APPTAINER_CACHEDIR="$PWD/data/containers"

# download and cache the Bowtie2 index using the containerized binary (in a clean environment)
apptainer exec \
    --cleanenv \
    --bind "$PWD:$PWD" \
    --env HOSTILE_CACHE_DIR="$PWD/data/dbs/hostile" \
    docker://quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0 \
    bash -c "export HOSTILE_CACHE_DIR=$PWD/data/dbs/hostile && hostile index fetch --name human-t2t-hla --bowtie2"

echo "[INFO] Hostile database setup complete."


# for when annotation is implemented:
# ==============================================================================
# EGGNOG-MAPPER DATABASE SETUP
# ==============================================================================

# ==============================================================================
# dbCAN3 DATABASE SETUP 
# ==============================================================================


echo "====================================================================="
echo " Setup Completed! The environment is ready."
echo " Now you might want to run: snakemake assemble --cores all --use-apptainer"
echo "====================================================================="