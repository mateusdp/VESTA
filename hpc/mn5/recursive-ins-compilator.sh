#!/bin/bash
# =============================================================================
# recursive-ins-compilator.sh
#
# Recursively compile a tree of nested LPJ-GUESS instruction files
# (linked by "import" directives) into a single combined .ins file.
#
# LPJ-GUESS .ins files use:
#   import "other_file.ins"
# to include parameters from other files. When running with MPI across
# many processes that may start from different directories, it is safer
# to flatten everything into one self-contained file.
#
# This script calls itself recursively for each import directive found.
# It appends to the output file, so the caller must clear it first:
#   echo -n > combined.ins
#   bash recursive-ins-compilator.sh main.ins /absolute/path/to/combined.ins
#
# Original author: Johan Nord, 2019-09-15
# Part of: lpj-guess-ntd-vesta HPC workflow
#
# Usage:
#   bash recursive-ins-compilator.sh <input.ins> <output.ins>
#
#   input.ins  — the top-level LPJ-GUESS instruction file
#   output.ins — absolute path to the combined output file (appended to)
#
# Note: output path ($2) MUST be absolute, since the script changes
# directories to resolve relative import paths.
# =============================================================================

currentinstructionfile=$1
newinstrucionfile=$2

echo "$currentinstructionfile"

script=$(realpath "$0")

cd "$(dirname "$currentinstructionfile")"
echo "! ********** HERE STARTS IMPORTED INS FILE $(pwd -P)/$(basename "$currentinstructionfile")" >> "$newinstrucionfile"
echo >> "$newinstrucionfile"

while IFS= read -r line; do

  if [[ $line == "import"* ]]; then
    echo "! ********** $line **********" >> "$newinstrucionfile"
    echo >> "$newinstrucionfile"

    nextins=$(echo "$line" | cut -d\" -f2-3 | sed 's/\"//')
    bash "$script" "$nextins" "$newinstrucionfile"

    echo "! //******** GOING BACK TO PARENT INS FILE $(pwd -P)/$(basename "$currentinstructionfile")" >> "$newinstrucionfile"
    echo >> "$newinstrucionfile"

  else
    echo "$line" >> "$newinstrucionfile"

  fi

done < "$(basename "$currentinstructionfile")"

echo >> "$newinstrucionfile"
echo "! //******** HERE ENDS IMPORTED INS FILE $(pwd -P)/$(basename "$currentinstructionfile")" >> "$newinstrucionfile"
