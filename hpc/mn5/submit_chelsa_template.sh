#!/bin/bash
# =============================================================================
# submit_chelsa_template.sh
#
# Production submission script for the global CHELSA v2-based LPJ-GUESS-NTD
# campaign on MareNostrum 5 (BSC).
#
# This script:
#   1. Builds a combined .ins file (flattening all imports)
#   2. Copies guess.cmd.template and fills in PLACEHOLDER_ tokens
#   3. Submits the job to SLURM
#
# It is called indirectly via setup_band.sh, which patches the band-specific
# variables below before this script runs.
#
# Usage:
#   bash submit.sh                        # (after setup_band.sh has patched it)
#   bash submit.sh -i /path/to/alt.ins    # override instruction file
#
# Part of: lpj-guess-ntd-vesta HPC workflow (MareNostrum 5)
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration — these are patched by setup_band.sh for each band
# =============================================================================

# SLURM resources
NPROCESS=11200          # 100 nodes x 112 cores/node
WALLTIME=48:00:00
PARTITION=gpp
ACCOUNT=ehpc336
QOS=gp_ehpc

# LPJ-GUESS binary
GUESS_BIN=/gpfs/home/senc/senc582781/LPJ_GUESS/NTD/LPJ-GUESS/build/guess
INPUT_MODULE=cf

# Input files — patched per band by setup_band.sh
INSFILE=/gpfs/scratch/ehpc336/LPJ_GUESS/ins/global_chelsa_v2_60S_30S/global_p-traits_daily_cf.ins
GRIDLIST=/gpfs/scratch/ehpc336/LPJ_GUESS/input/chelsa_bands/gridlist_60S_30S.txt
SOILFILE=/gpfs/scratch/ehpc336/LPJ_GUESS/input/chelsa_bands/soilmap_60S_30S.dat

# Climate files (for verification printout only; paths in .ins are authoritative)
CHELSA_TEMP=/gpfs/scratch/ehpc336/LPJ_GUESS/input/chelsa_bands/temp_60S_30S.nc
CHELSA_PREC=/gpfs/scratch/ehpc336/LPJ_GUESS/input/chelsa_bands/prec_60S_30S.nc
CHELSA_RAD=/gpfs/scratch/ehpc336/LPJ_GUESS/input/chelsa_bands/rad_60S_30S.nc

RUNNAME=global_chelsa_v2_60S_30S

# Template for the SLURM batch script
GUESS_CMD_TEMPLATE=/gpfs/scratch/ehpc336/LPJ_GUESS/scripts/guess.cmd.template

# =============================================================================
# Argument handling
# =============================================================================

while getopts ":n:s:i:" opt; do
    case $opt in
        n ) name=$OPTARG ;;
        s ) submit_vars_file=$OPTARG ;;
        i ) ins=$OPTARG ;;
    esac
done

if [ -n "${submit_vars_file:-}" ]; then
    source "$submit_vars_file"
fi

if [ -n "${ins:-}" ]; then
    INSFILE=$ins
fi

INSFILE=$(readlink -f "$INSFILE")
COMBINED_INSFILE=$(realpath ${RUNNAME}_combined_guess.ins)
GRIDLIST_FILENAME=$(basename "$GRIDLIST")
INSDIR=$(dirname "$INSFILE")

# =============================================================================
# Build guess.cmd from template
# =============================================================================

cp "$GUESS_CMD_TEMPLATE" guess.cmd

sed -i \
    -e "s|PLACEHOLDER_NPROCESS|${NPROCESS}|g" \
    -e "s|PLACEHOLDER_WALLTIME|${WALLTIME}|g" \
    -e "s|PLACEHOLDER_RUNNAME|${RUNNAME}|g" \
    -e "s|PLACEHOLDER_PARTITION|${PARTITION}|g" \
    -e "s|PLACEHOLDER_ACCOUNT|${ACCOUNT}|g" \
    -e "s|PLACEHOLDER_QOS|${QOS}|g" \
    -e "s|PLACEHOLDER_COMBINED_INSFILE|${COMBINED_INSFILE}|g" \
    -e "s|PLACEHOLDER_INSDIR|${INSDIR}|g" \
    -e "s|PLACEHOLDER_INSFILE|${INSFILE}|g" \
    -e "s|PLACEHOLDER_GRIDLIST_FILENAME|${GRIDLIST_FILENAME}|g" \
    -e "s|PLACEHOLDER_GRIDLIST|${GRIDLIST}|g" \
    -e "s|PLACEHOLDER_SOILFILE|${SOILFILE}|g" \
    -e "s|PLACEHOLDER_GUESS_BIN|${GUESS_BIN}|g" \
    -e "s|PLACEHOLDER_INPUT_MODULE|${INPUT_MODULE}|g" \
    guess.cmd

echo ""
echo "============================================="
echo "Production run: $RUNNAME"
echo "Nodes:    100 ($NPROCESS cores)"
echo "Queue:    $QOS ($WALLTIME)"
echo "Gridlist: $GRIDLIST"
echo "Soil:     $SOILFILE"
echo "INS:      $INSFILE"
echo "============================================="
echo ""

sbatch -J "$RUNNAME" guess.cmd | awk '{print $NF}'
