#!/bin/bash
# =============================================================================
# progress.sh
#
# Check progress of a running LPJ-GUESS job by sampling guess.log files
# from a subset of run directories. Shows the last line from sampled
# processes, giving an overview of simulation progress.
#
# Usage:
#   bash progress.sh [N_SAMPLE]
#
#   N_SAMPLE: number of processes to sample (default: 10)
#
# Part of: lpj-guess-ntd-vesta HPC workflow
# =============================================================================

N_SAMPLE=${1:-10}

NPROCESS=$(ls -d run*/ 2>/dev/null | wc -l)

if [ "$NPROCESS" -eq 0 ]; then
    echo "No run directories found."
    exit 1
fi

echo "=========================================="
echo "LPJ-GUESS progress check"
echo "  Total processes: $NPROCESS"
echo "  Sampling: $N_SAMPLE processes"
echo "  Time: $(date)"
echo "=========================================="

# Sample evenly spaced processes
STEP=$((NPROCESS / N_SAMPLE))
if [ "$STEP" -lt 1 ]; then STEP=1; fi

for ((i=1; i<=NPROCESS; i+=STEP)); do
    if [ -f "run${i}/guess.log" ]; then
        echo ""
        echo "--- run${i}/guess.log ---"
        tail -3 "run${i}/guess.log"
    fi
done

echo ""
echo "=========================================="
# Count completed processes (those with "Finished" in their log)
FINISHED=$(grep -l "Finished" run*/guess.log 2>/dev/null | wc -l)
echo "Completed: $FINISHED / $NPROCESS"
echo "=========================================="
