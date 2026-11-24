#!/bin/bash
# =============================================================================
# setup_band.sh
#
# Prepare instruction files and submission script for one latitude band.
# Creates a band-specific directory, copies the .ins template files,
# patches all band-dependent paths, and copies the submission script
# with the correct configuration.
#
# Usage:
#   bash setup_band.sh <BAND_NAME>
#   bash setup_band.sh 60N_90N
#
# This is the first step in the per-band production workflow:
#   1. bash setup_band.sh <BAND>
#   2. cd <RUN_DIR> && bash submit.sh
#
# Part of: lpj-guess-ntd-vesta HPC workflow (MareNostrum 5)
# =============================================================================

set -euo pipefail

BAND=$1
if [ -z "$BAND" ]; then
    echo "Usage: bash setup_band.sh <BAND_NAME>"
    echo "Example: bash setup_band.sh 60N_90N"
    echo ""
    echo "Available bands:"
    echo "  60N_90N  50N_60N  40N_50N  30N_40N  20N_30N"
    echo "  10N_20N  00N_10N  10S_00N  20S_10S  30S_20S  60S_30S"
    exit 1
fi

# =============================================================================
# Configuration — UPDATE THESE PATHS for your HPC environment
# =============================================================================

# Base directory containing per-band climate/soil/gridlist files
# Expected structure: BANDS_DIR/{BAND}/temp_{BAND}.nc, prec_{BAND}.nc, etc.
BANDS_DIR="/gpfs/scratch/ehpc336/LPJ_GUESS/input/chelsa_bands"

# Directory where band-specific .ins files will be created
INS_DIR="/gpfs/scratch/ehpc336/LPJ_GUESS/ins/global_chelsa_v2_${BAND}"

# Run directory (where the job will be submitted from)
RUN_DIR="${INS_DIR}/results_production"

# Template directory containing the master .ins files
TEMPLATE_DIR="/gpfs/scratch/ehpc336/LPJ_GUESS/ins/global_chelsa_v2_template"

# Template scripts
TEMPLATE_SUBMIT="/gpfs/scratch/ehpc336/LPJ_GUESS/scripts/submit_chelsa_template.sh"
TEMPLATE_GUESSCMD="/gpfs/scratch/ehpc336/LPJ_GUESS/scripts/guess.cmd.template"

# =============================================================================
# Setup
# =============================================================================

# ---- Create directories ---
mkdir -p "$RUN_DIR"

# ---- Copy all .ins files ---
cp "${TEMPLATE_DIR}"/*.ins "${INS_DIR}/"

# ---- Patch main .ins file with band-specific paths ---
# The template .ins uses 60S_30S as the default band; substitute all
# band-specific references to point to the correct files.
sed -i \
    -e "s|gridlist_60S_30S.txt|gridlist_${BAND}.txt|g" \
    -e "s|temp_v2_60S_30S_reorder.nc|${BAND}/temp_${BAND}.nc|g" \
    -e "s|prec_v2_60S_30S_reorder.nc|${BAND}/prec_${BAND}.nc|g" \
    -e "s|rad_v2_60S_30S_reorder.nc|${BAND}/rad_${BAND}.nc|g" \
    -e "s|soilmap_60S_30S.dat|${BAND}/soilmap_${BAND}.dat|g" \
    "${INS_DIR}/global_p-traits_daily_cf.ins"

# ---- Copy and patch submit script ---
cp "$TEMPLATE_SUBMIT" "${RUN_DIR}/submit.sh"
sed -i \
    -e "s|^INSFILE=.*|INSFILE=${INS_DIR}/global_p-traits_daily_cf.ins|" \
    -e "s|^GRIDLIST=.*|GRIDLIST=${BANDS_DIR}/gridlist_${BAND}.txt|" \
    -e "s|^SOILFILE=.*|SOILFILE=${BANDS_DIR}/${BAND}/soilmap_${BAND}.dat|" \
    -e "s|^RUNNAME=.*|RUNNAME=global_chelsa_v2_${BAND}|" \
    -e "s|^CHELSA_TEMP=.*|CHELSA_TEMP=${BANDS_DIR}/${BAND}/temp_${BAND}.nc|" \
    -e "s|^CHELSA_PREC=.*|CHELSA_PREC=${BANDS_DIR}/${BAND}/prec_${BAND}.nc|" \
    -e "s|^CHELSA_RAD=.*|CHELSA_RAD=${BANDS_DIR}/${BAND}/rad_${BAND}.nc|" \
    "${RUN_DIR}/submit.sh"

cd "$RUN_DIR"

# ---- Summary ---
echo "=========================================="
echo "Band: $BAND"
echo "=========================================="
echo ""
echo "INS files copied:"
ls -1 "${INS_DIR}"/*.ins
echo ""
echo "Climate paths in .ins:"
grep -E "file_temp|file_prec|file_insol|file_gridlist|file_soildata" "${INS_DIR}/global_p-traits_daily_cf.ins"
echo ""
echo "Submit script:"
grep -E "^NPROCESS|^WALLTIME|^INSFILE|^GRIDLIST|^SOILFILE|^RUNNAME|^CHELSA_" "${RUN_DIR}/submit.sh"
echo ""
echo "To submit:  cd ${RUN_DIR} && bash submit.sh"
