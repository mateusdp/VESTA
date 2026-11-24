#!/bin/bash
# =============================================================================
# join_compress.sh
#
# Join per-process LPJ-GUESS output files and compress with gzip.
# Run from the job's submit directory (where run1/, run2/, ... exist).
#
# For each .out file found in run1/:
#   1. Copies run1's file (with header) as the base
#   2. Appends data rows (skipping headers) from run2..runN
#   3. Compresses the joined file with gzip
#   4. Deletes the per-process copies to free disk space
#
# Usage:
#   bash join_compress.sh [NPROCESS]
#
#   NPROCESS defaults to the number of run*/ directories found.
#
# This can be run as a standalone script after the GUESS job finishes,
# or submitted as a SLURM job for large runs.
#
# Part of: lpj-guess-ntd-vesta HPC workflow
# =============================================================================

set -euo pipefail

# Determine NPROCESS
if [ -n "${1:-}" ]; then
    NPROCESS=$1
else
    NPROCESS=$(ls -d run*/ 2>/dev/null | wc -l)
fi

if [ "$NPROCESS" -eq 0 ]; then
    echo "ERROR: No run directories found." >&2
    exit 1
fi

echo "=========================================="
echo "Join + compress output files"
echo "  Run directories: $NPROCESS"
echo "  Start: $(date)"
echo "=========================================="

# Find output files from run1
pushd run1 &> /dev/null
OUTFILES=$(ls *.out 2>/dev/null || true)
popd &> /dev/null

if [ -z "$OUTFILES" ]; then
    echo "No .out files found in run1/. Nothing to do."
    exit 0
fi

for file in $OUTFILES; do
    echo -n "  $file: appending..."
    cp "run1/$file" "$file"

    for ((i=2; i<=NPROCESS; i++)); do
        [ -f "run$i/$file" ] && tail -n +2 "run$i/$file"
    done >> "$file"

    RAW=$(du -sh "$file" | cut -f1)
    echo -n " ($RAW) compressing..."
    gzip "$file"
    GZ=$(du -sh "$file.gz" | cut -f1)

    echo -n " ($GZ) cleaning..."
    for ((i=1; i<=NPROCESS; i++)); do
        rm -f "run$i/$file"
    done
    echo " done."
done

echo ""
echo "Cleaning up empty run directories..."
rm -rf run*/

echo ""
echo "=========================================="
echo "Done at $(date)"
echo "=========================================="
