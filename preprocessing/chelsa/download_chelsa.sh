#!/bin/bash
# ============================================================================
# download_chelsa.sh
# Download CHELSA v2.1 monthly climate data (GeoTIFF) from ENVICLOUD
#
# Variables: tas (temperature), pr (precipitation), rsds (radiation)
# Resolution: 30 arcsec (~1 km)
# Period:     1979-02 to 2021-12 (CHELSA v2.1 availability)
#
# NOTE: CHELSA v2.1 does NOT include January 1979. By default this script
# also downloads a CHELSA v1.2 substitute for January 1979, which is needed
# to have a complete 1979-01 to 2021-12 timeseries. The preprocessing R
# script detects v1 files automatically (different scale factor).
# Use -s to skip the v1 January 1979 download.
#
# Usage:
#   ./download_chelsa.sh -o /path/to/output [-v tas,pr,rsds] [-y 1979-2021]
#
# The script skips files that already exist (resume-safe) and retries
# failed downloads. Designed for large batch downloads on HPC login nodes
# or data transfer machines.
#
# Part of: lpj-guess-ntd-vesta preprocessing pipeline
# ============================================================================

set -euo pipefail

# --- defaults ---------------------------------------------------------------
BASE_URL="https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly"
VARIABLES="tas,pr,rsds"
YEAR_START=1979
YEAR_END=2021
OUTDIR=""
MAX_RETRIES=5
RETRY_WAIT=30          # seconds between retries
PARALLEL=1             # number of parallel downloads (wget default)
LOG_FILE=""
DRY_RUN=false
SKIP_JAN1979=false

# --- usage ------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Download CHELSA v2.1 monthly climate GeoTIFFs from ENVICLOUD.

Required:
  -o DIR        Output base directory (created if needed).
                Files are saved as DIR/{variable}/{year}/CHELSA_{var}_{MM}_{YYYY}_V.2.1.tif

Optional:
  -v VARS       Comma-separated variables to download (default: tas,pr,rsds)
  -y RANGE      Year range as START-END (default: 1979-2021)
  -r N          Max retries per file (default: 5)
  -w SEC        Wait between retries in seconds (default: 30)
  -l FILE       Log file path (default: OUTDIR/download_chelsa.log)
  -s            Skip January 1979 v1 substitute download
  -n            Dry run: print URLs without downloading
  -h            Show this help

Examples:
  # Download all three variables, full period
  $(basename "$0") -o /scratch/chelsa_raw

  # Download only temperature, 2010-2021
  $(basename "$0") -o /scratch/chelsa_raw -v tas -y 2010-2021

  # Dry run to check URLs
  $(basename "$0") -o /scratch/chelsa_raw -n
EOF
    exit 0
}

# --- parse arguments --------------------------------------------------------
while getopts "o:v:y:r:w:l:snh" opt; do
    case "$opt" in
        o) OUTDIR="$OPTARG" ;;
        v) VARIABLES="$OPTARG" ;;
        y) YEAR_START="${OPTARG%-*}"; YEAR_END="${OPTARG#*-}" ;;
        r) MAX_RETRIES="$OPTARG" ;;
        w) RETRY_WAIT="$OPTARG" ;;
        l) LOG_FILE="$OPTARG" ;;
        s) SKIP_JAN1979=true ;;
        n) DRY_RUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    echo "ERROR: Output directory (-o) is required." >&2
    usage
fi

# --- setup ------------------------------------------------------------------
mkdir -p "$OUTDIR"
if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="${OUTDIR}/download_chelsa.log"
fi

log() {
    local msg
    msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg" | tee -a "$LOG_FILE"
}

log "=========================================="
log "CHELSA v2.1 download started"
log "Variables:  $VARIABLES"
log "Period:     ${YEAR_START}-${YEAR_END}"
log "Output:     $OUTDIR"
log "Dry run:    $DRY_RUN"
log "=========================================="

# --- counters ---------------------------------------------------------------
TOTAL=0
SKIPPED=0
DOWNLOADED=0
FAILED=0

# --- download loop ----------------------------------------------------------
IFS=',' read -ra VARS <<< "$VARIABLES"

for var in "${VARS[@]}"; do
    for year in $(seq "$YEAR_START" "$YEAR_END"); do

        # CHELSA v2.1 starts at February 1979
        if [[ "$year" -eq 1979 ]]; then
            MONTH_START=2
        else
            MONTH_START=1
        fi

        for month in $(seq "$MONTH_START" 12); do
            MM=$(printf "%02d" "$month")
            FILENAME="CHELSA_${var}_${MM}_${year}_V.2.1.tif"
            URL="${BASE_URL}/${var}/${year}/${FILENAME}"
            DESTDIR="${OUTDIR}/${var}/${year}"
            DESTFILE="${DESTDIR}/${FILENAME}"

            TOTAL=$((TOTAL + 1))

            # --- dry run ---
            if $DRY_RUN; then
                echo "$URL"
                continue
            fi

            # --- skip existing files ---
            if [[ -f "$DESTFILE" ]]; then
                # Check file is not empty / truncated (at least 1 MB for a CHELSA tile)
                FSIZE=$(stat -c%s "$DESTFILE" 2>/dev/null || echo 0)
                if [[ "$FSIZE" -gt 1048576 ]]; then
                    SKIPPED=$((SKIPPED + 1))
                    continue
                else
                    log "WARN: ${FILENAME} exists but is too small (${FSIZE} bytes), re-downloading"
                    rm -f "$DESTFILE"
                fi
            fi

            mkdir -p "$DESTDIR"

            # --- download with retries ---
            SUCCESS=false
            for attempt in $(seq 1 "$MAX_RETRIES"); do
                if wget -q --timeout=120 --tries=1 -O "$DESTFILE" "$URL" 2>>"$LOG_FILE"; then
                    # Verify downloaded file is not empty
                    FSIZE=$(stat -c%s "$DESTFILE" 2>/dev/null || echo 0)
                    if [[ "$FSIZE" -gt 1048576 ]]; then
                        SUCCESS=true
                        break
                    else
                        log "WARN: ${FILENAME} downloaded but too small (${FSIZE} bytes), attempt ${attempt}/${MAX_RETRIES}"
                        rm -f "$DESTFILE"
                    fi
                else
                    log "WARN: Download failed for ${FILENAME}, attempt ${attempt}/${MAX_RETRIES}"
                    rm -f "$DESTFILE"
                fi

                if [[ "$attempt" -lt "$MAX_RETRIES" ]]; then
                    sleep "$RETRY_WAIT"
                fi
            done

            if $SUCCESS; then
                DOWNLOADED=$((DOWNLOADED + 1))
                # Progress every 50 files
                if [[ $((DOWNLOADED % 50)) -eq 0 ]]; then
                    log "Progress: ${DOWNLOADED} downloaded, ${SKIPPED} skipped, ${FAILED} failed (of ${TOTAL} processed)"
                fi
            else
                FAILED=$((FAILED + 1))
                log "ERROR: Failed to download ${FILENAME} after ${MAX_RETRIES} attempts"
                echo "$URL" >> "${OUTDIR}/failed_downloads.txt"
            fi
        done
    done
done

# --- January 1979 from CHELSA v1.2 -----------------------------------------
# CHELSA v2.1 starts at February 1979. For a complete 1979-01 to 2021-12
# timeseries, we download January 1979 from v1.2.
# The R preprocessing script detects these files by their unscaled values.
V1_BASE_URL="https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly"

if ! $SKIP_JAN1979 && [[ "$YEAR_START" -le 1979 ]]; then
    log "Downloading January 1979 substitute files from CHELSA v1.2..."
    for var in "${VARS[@]}"; do
        # v1.2 naming convention differs slightly from v2.1
        V1_FILENAME="CHELSA_${var}_01_1979_V1.2_land.tif"
        V1_URL="${V1_BASE_URL}/${var}/1979/${V1_FILENAME}"
        DESTDIR="${OUTDIR}/${var}/1979"
        DESTFILE="${DESTDIR}/${V1_FILENAME}"

        TOTAL=$((TOTAL + 1))

        if $DRY_RUN; then
            echo "$V1_URL"
            continue
        fi

        if [[ -f "$DESTFILE" ]]; then
            FSIZE=$(stat -c%s "$DESTFILE" 2>/dev/null || echo 0)
            if [[ "$FSIZE" -gt 1048576 ]]; then
                SKIPPED=$((SKIPPED + 1))
                log "  ${var} Jan 1979 (v1.2): already exists, skipping"
                continue
            fi
        fi

        mkdir -p "$DESTDIR"
        log "  Downloading ${V1_FILENAME}..."

        SUCCESS=false
        for attempt in $(seq 1 "$MAX_RETRIES"); do
            if wget -q --timeout=120 --tries=1 -O "$DESTFILE" "$V1_URL" 2>>"$LOG_FILE"; then
                FSIZE=$(stat -c%s "$DESTFILE" 2>/dev/null || echo 0)
                if [[ "$FSIZE" -gt 1048576 ]]; then
                    SUCCESS=true
                    break
                fi
            fi
            rm -f "$DESTFILE"
            if [[ "$attempt" -lt "$MAX_RETRIES" ]]; then
                sleep "$RETRY_WAIT"
            fi
        done

        if $SUCCESS; then
            DOWNLOADED=$((DOWNLOADED + 1))
            log "  ${var} Jan 1979 (v1.2): OK"
        else
            FAILED=$((FAILED + 1))
            log "ERROR: Failed to download ${V1_FILENAME}"
            log "  You may need to find this file manually from CHELSA v1.2"
            echo "$V1_URL" >> "${OUTDIR}/failed_downloads.txt"
        fi
    done
fi

# --- summary ----------------------------------------------------------------
log "=========================================="
log "Download complete"
log "Total files:  $TOTAL"
log "Downloaded:   $DOWNLOADED"
log "Skipped:      $SKIPPED (already exist)"
log "Failed:       $FAILED"
if [[ "$FAILED" -gt 0 ]]; then
    log "Failed URLs saved to: ${OUTDIR}/failed_downloads.txt"
    log "Re-run the script to retry failed downloads"
fi
log "=========================================="
