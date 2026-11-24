#!/bin/bash
# =============================================================================
# gather_compressed.sh
#
# Download and concatenate gzipped LPJ-GUESS output from MN5 to a local
# machine. Run this from mrpharmacist (or any machine with sshfs) after a
# band finishes on the HPC.
#
# When LPJ-GUESS is compiled with -D COMPRESS_OUTPUT, each MPI process
# writes compressed output directly: a .out.gz data file and a .out.hdr
# header file. This script mounts the remote run directory via sshfs and
# concatenates the per-process .gz streams into single multi-stream .gz
# files (gzip natively supports this).
#
# This is the counterpart to join_compress.sh, which works on-cluster
# with uncompressed .out files. Use this script when COMPRESS_OUTPUT is
# enabled (the production configuration); use join_compress.sh when it
# is not.
#
# Prerequisites:
#   - sshfs (FUSE-based SSH filesystem)
#   - SSH key-based access to the HPC transfer node
#
# Usage:
#   bash gather_compressed.sh -d <remote_directory> [-u <user>] [-r <remote>] [-v <variables>] [-c]
#
# Examples:
#   # Download all variables from 60S_30S band:
#   bash gather_compressed.sh -d /gpfs/scratch/ehpc336/LPJ_GUESS/ins/global_chelsa_v2_60S_30S/results_production
#
#   # Download specific variables only:
#   bash gather_compressed.sh -d /path/to/results -v "cmass agpp lai"
#
#   # Download and clean up remote run directories after:
#   bash gather_compressed.sh -d /path/to/results -c
#
# Original authors: Etienne Tourigny, Matthew Forrest
# Part of: lpj-guess-ntd-vesta HPC workflow
# =============================================================================

if [[ "$1" == "-h" || "$1" == "--help" || $# == 0 ]]; then
    echo "Usage: bash gather_compressed.sh -d <remote_directory> [-u <user>] [-r <remote>] [-v <variables>] [-c]"
    echo "  -d : Remote run directory path"
    echo "  -u : Remote username (default: senc582781)"
    echo "  -r : Remote HPC address (default: transfer2.bsc.es)"
    echo "  -v : Variables to gather, space-separated (default: ALL)"
    echo "  -c : Clean up remote run directories after download"
    exit 0
fi

# =============================================================================
# Configuration — UPDATE THESE for your HPC environment
# =============================================================================

USERNAME="senc582781"
REMOTEHPC="transfer2.bsc.es"
VARS="ALL"
CLEANUP=false

# =============================================================================
# Argument parsing
# =============================================================================

if [[ $# == 1 ]]; then
    DIR="$1"
else
    while getopts d:u:r:v:c flag; do
        case "${flag}" in
            u) USERNAME=${OPTARG};;
            r) REMOTEHPC=${OPTARG};;
            v) VARS=${OPTARG};;
            d) DIR=${OPTARG};;
            c) CLEANUP=true;;
        esac
    done
fi

echo "=========================================="
echo "Gathering LPJ-GUESS output"
echo "  Remote: ${USERNAME}@${REMOTEHPC}:${DIR}"
echo "  Variables: ${VARS}"
echo "  Cleanup: ${CLEANUP}"
echo "=========================================="

# =============================================================================
# Mount remote directory via sshfs
# =============================================================================

MOUNTPOINT="/tmp/remote_$$"
mkdir -p "$MOUNTPOINT"
sshfs "${USERNAME}@${REMOTEHPC}:${DIR}" "${MOUNTPOINT}"

echo "Mounted remote directory."

# =============================================================================
# Build file list
# =============================================================================

if [ "$VARS" == "ALL" ]; then
    pushd "${MOUNTPOINT}/run1" &> /dev/null
    outfiles_expanded=$(echo *.out.gz)
    popd &> /dev/null
else
    outfiles_expanded=""
    for var in $VARS; do
        outfiles_expanded="$outfiles_expanded ${var}.out.gz"
    done
fi

echo "Files to gather: ${outfiles_expanded}"

# Copy .ins files for reference
echo "Copying instruction files..."
cp "${MOUNTPOINT}"/*.ins . 2>/dev/null || true
echo "Done."

# Get number of run directories
NRUNDIRS=$(ls -d "${MOUNTPOINT}"/run*/ 2>/dev/null | wc -l)
echo "Found ${NRUNDIRS} run directories."
echo ""

set -e

# =============================================================================
# Concatenate per-process compressed output
# =============================================================================
# With COMPRESS_OUTPUT, each process writes:
#   - <var>.out.hdr  (plain-text header line)
#   - <var>.out.gz   (gzipped data rows, no header)
#
# Gzip supports multi-stream concatenation: cat file1.gz file2.gz > combined.gz
# decompresses to the concatenation of their contents. So we gzip the header
# and then append all per-process .gz files sequentially.

function append_files {
    local number_of_jobs=$1
    local file=$2

    # Compress the header file and use it as the base
    cp "${MOUNTPOINT}/run1/${file/.gz/.hdr}" "${file/.gz/}"
    gzip "${file/.gz/}"

    # Append all per-process .gz streams
    local i=""
    for ((i=1; i <= number_of_jobs; i++)); do
        if (( i % 500 == 0 )); then
            echo -n "${i} "
        fi
        if [ -f "${MOUNTPOINT}/run${i}/${file}" ]; then
            cat "${MOUNTPOINT}/run${i}/${file}" >> "$file"
        fi
    done
}

for file in $outfiles_expanded; do
    # Skip if already downloaded
    if [ -f "$file" ]; then
        echo "  **** Skipping: ${file} (already exists)"
        continue
    fi
    echo "  **** Gathering: ${file}"
    TSTART=$(date +%s)
    append_files "$NRUNDIRS" "$file"
    TEND=$(date +%s)
    SIZE=$(du -sh "$file" | cut -f1)
    echo ""
    echo "  **** Done: ${SIZE} in $((TEND-TSTART))s"
    echo ""
done

# =============================================================================
# Cleanup and unmount
# =============================================================================

if [ "$CLEANUP" = true ]; then
    echo "Cleaning up remote run directories..."
    rm -rf "${MOUNTPOINT}"/run*/
    echo "Remote cleanup complete."
fi

sleep 5
fusermount -u "$MOUNTPOINT"
rmdir "$MOUNTPOINT"

echo ""
echo "=========================================="
echo "Gather complete."
echo "=========================================="
