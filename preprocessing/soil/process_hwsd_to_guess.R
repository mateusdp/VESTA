#!/usr/bin/env Rscript
# =============================================================================
# process_hwsd_to_guess.R
#
# Create per-band soil input files for LPJ-GUESS from the Harmonized World
# Soil Database (HWSD v1.21).
#
# For each latitude band, reads lon/lat from the band's climate NetCDF
# (produced by process_chelsa_to_guess.R), extracts HWSD soil properties at
# those EXACT coordinates, and writes a soilmap .dat file. This guarantees
# coordinate matching — no searchradius needed in the GUESS instruction file.
#
# Output per band:
#   soilmap_{band}.dat — tab-separated: lon lat clay silt sand orgC CN pH cellfraction
#
# Dependencies:
#   R packages: terra, data.table, ncdf4
#
# Data requirements:
#   1. HWSD raster (hwsd.bil + hwsd.hdr)
#      Download from: https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v12/en/
#      The .bil/.hdr files contain mapping unit IDs at 30 arcsec resolution.
#
#   2. HWSD attribute table (hwsd_data.csv)
#      Extracted from the HWSD Microsoft Access database (.mdb) or available
#      as a pre-extracted CSV. Contains soil properties per mapping unit.
#      Required columns: MU_GLOBAL, S_USDA_TEX_CLASS, SHARE,
#                        S_CLAY, S_SILT, S_SAND, S_OC, S_PH_H2O
#
#   3. Per-band climate NetCDFs (temp_{band}.nc)
#      From the CHELSA preprocessing step — used only to read the grid
#      coordinates, not climate data.
#
# Part of: lpj-guess-ntd-vesta preprocessing pipeline
# =============================================================================

rm(list = ls())
stime <- Sys.time()

library(terra)
library(data.table)
library(ncdf4)

# ---- Configuration ----------------------------------------------------------
hwsd_raster_path <- "/home/mdantas/map_data/hwsd/hwsd.bil"
hwsd_csv_path    <- "/home/mdantas/map_data/hwsd/hwsd_data.csv"

# Directory containing band subdirectories with climate NetCDFs
bands_climate_dir <- "/home/mdantas/chelsa_processing/global_v2"

# Output directory for soil files (default: same band subdirectories)
soil_output_dir <- bands_climate_dir

# Band definitions — must match process_chelsa_to_guess.R
band_defs <- list(
  c("60N_90N",   60,  90),
  c("50N_60N",   50,  60),
  c("40N_50N",   40,  50),
  c("30N_40N",   30,  40),
  c("20N_30N",   20,  30),
  c("10N_20N",   10,  20),
  c("00N_10N",    0,  10),
  c("10S_00N",  -10,   0),
  c("20S_10S",  -20, -10),
  c("30S_20S",  -30, -20),
  c("60S_30S",  -60, -30)
)

# =============================================================================
# STEP 1: Load HWSD raster
# =============================================================================
cat("\n=== Step 1: Loading HWSD raster ===\n")

hwsd <- rast(hwsd_raster_path)
crs(hwsd) <- "EPSG:4326"
cat("  Resolution:", res(hwsd), "\n")

# =============================================================================
# STEP 2: Process HWSD attribute table (once, shared across all bands)
# =============================================================================
cat("\n=== Step 2: Processing HWSD attribute table ===\n")

records <- fread(hwsd_csv_path,
                 select = c("MU_GLOBAL", "S_USDA_TEX_CLASS", "SHARE",
                            "S_CLAY", "S_SILT", "S_SAND", "S_OC", "S_PH_H2O"))

# For each mapping unit, select texture class with highest total share
share_summary <- records[, .(total_share = sum(SHARE, na.rm = TRUE)),
                         by = .(MU_GLOBAL, S_USDA_TEX_CLASS)]
best_tex <- share_summary[, .SD[which.max(total_share)], by = MU_GLOBAL]
records_keyed <- records[best_tex, on = .(MU_GLOBAL, S_USDA_TEX_CLASS), mult = "first"]

soil_lookup <- unique(records_keyed[, .(MU_GLOBAL,
                                        clay = S_CLAY / 100,
                                        silt = S_SILT / 100,
                                        sand = S_SAND / 100,
                                        orgC = S_OC / 100,
                                        pH   = S_PH_H2O)],
                      by = "MU_GLOBAL")

cat("  Lookup table:", nrow(soil_lookup), "unique mapping units\n")

rm(records, share_summary, best_tex, records_keyed)
gc()

# =============================================================================
# STEP 3: Process each band
# =============================================================================
cat("\n=== Step 3: Processing bands ===\n")

for (bd in band_defs) {
  bname     <- bd[1]
  lat_south <- as.numeric(bd[2])
  lat_north <- as.numeric(bd[3])

  # Find the band's climate NetCDF (use temp file for coordinates)
  temp_nc <- file.path(bands_climate_dir, bname, paste0("temp_", bname, ".nc"))
  if (!file.exists(temp_nc)) {
    cat("  Band", bname, "- climate file not found, skipping\n")
    next
  }

  # Check if output already exists
  outfile <- file.path(soil_output_dir, bname, paste0("soilmap_", bname, ".dat"))
  if (file.exists(outfile) && file.size(outfile) > 1e6) {
    cat("  Band", bname, "- soilmap already exists, skipping\n")
    next
  }

  cat("\n--- Band:", bname, "---\n")

  # Read EXACT coordinates from the band's climate NetCDF
  nc <- nc_open(temp_nc)
  nc_lon <- ncvar_get(nc, "lon")
  nc_lat <- ncvar_get(nc, "lat")
  nc_close(nc)

  n_cells <- length(nc_lon)
  cat("  Cells from NetCDF:", format(n_cells, big.mark = ","), "\n")

  # Extract HWSD mapping units at these exact locations (in chunks to manage memory)
  chunk_size <- 5000000
  n_chunks <- ceiling(n_cells / chunk_size)
  mu_values <- integer(n_cells)

  for (ch in seq_len(n_chunks)) {
    idx_start <- (ch - 1) * chunk_size + 1
    idx_end   <- min(ch * chunk_size, n_cells)
    cat("\r  Extracting chunk", ch, "/", n_chunks, "...")
    flush.console()

    pts <- vect(
      data.frame(lon = nc_lon[idx_start:idx_end],
                 lat = nc_lat[idx_start:idx_end]),
      geom = c("lon", "lat"), crs = "EPSG:4326"
    )

    extracted <- terra::extract(hwsd, pts)
    mu_values[idx_start:idx_end] <- extracted[, 2]
  }
  cat("\n")

  # Build data.table with EXACT NetCDF coordinates
  soil_grid <- data.table(
    lon       = nc_lon,
    lat       = nc_lat,
    MU_GLOBAL = as.integer(mu_values)
  )

  # Join with soil properties
  soil_grid <- soil_lookup[soil_grid, on = "MU_GLOBAL"]

  # Fill NAs with band-level means (cells over water, ice, or missing HWSD data)
  n_na <- sum(is.na(soil_grid$clay))
  mean_clay <- soil_grid[, mean(clay, na.rm = TRUE)]
  mean_sand <- soil_grid[, mean(sand, na.rm = TRUE)]
  mean_silt <- 1 - mean_clay - mean_sand
  mean_orgC <- soil_grid[, mean(orgC, na.rm = TRUE)]
  mode_pH   <- soil_grid[!is.na(pH), .N, by = pH][which.max(N), pH]

  soil_grid[is.na(clay), clay := mean_clay]
  soil_grid[is.na(silt), silt := mean_silt]
  soil_grid[is.na(sand), sand := mean_sand]
  soil_grid[is.na(orgC), orgC := mean_orgC]
  soil_grid[is.na(pH),   pH   := mode_pH]

  # Assemble output with fixed values:
  #   CN = 12          — default soil C:N ratio
  #   cellfraction = 0.45 — fraction of cell that is vegetated
  output <- soil_grid[, .(lon, lat, clay, silt, sand, orgC,
                          CN = 12, pH, cellfraction = 0.45)]
  output <- output[is.finite(clay) & is.finite(silt) &
                     is.finite(sand) & is.finite(orgC) & is.finite(pH)]

  # Write tab-separated file
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  write.table(format(output, digits = 22), file = outfile, sep = "\t",
              quote = FALSE, row.names = FALSE)

  cat("  Valid soil:", format(n_cells - n_na, big.mark = ","),
      " | Filled:", format(n_na, big.mark = ","), "\n")
  cat("  Written:", outfile,
      "(", round(file.size(outfile) / 1e6), "MB )\n")

  rm(soil_grid, output, mu_values, nc_lon, nc_lat)
  gc()
}

cat("\n=== Done ===\n")
cat("Total time:", format(Sys.time() - stime), "\n")
