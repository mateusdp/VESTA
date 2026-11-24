#!/usr/bin/env Rscript
# ============================================================================
# process_chelsa_to_guess.R
#
# Convert CHELSA v2.1 monthly GeoTIFFs (tas, pr, rsds) into per-band
# NetCDF climate inputs for LPJ-GUESS.
#
# Output per band:
#   temp_{band}.nc   — air temperature [K]
#   prec_{band}.nc   — precipitation [kg m-2 month-1]
#   rad_{band}.nc    — shortwave radiation [W m-2]
#   gridlist_{band}.txt — grid cell list for GUESS instruction file
#
# NetCDF layout:
#   Dimensions:  grid_ID (UNLIMITED), time (fixed = nyears*12)
#   Variables:   var(grid_ID, time)  — each cell's full timeseries contiguous
#   Calendar:    365_day (no leap years), "days since 1979-01-01"
#
#   grid_ID as the record (UNLIMITED) dimension is what LPJ-GUESS expects
#   for efficient per-MPI-rank reading. Do NOT swap to time-unlimited.
#
# Dependencies: R packages terra, stringr, data.table, RNetCDF
#
# CHELSA v2.1 data availability:
#   - tas and pr: February 1979 – December 2021 (515 monthly files each)
#   - rsds: February 1979 – December 2021 (515 monthly files)
#   - January 1979 is NOT available in v2.1. A substitute from CHELSA v1.2
#     must be placed in the same directories. The script detects v1 files
#     by their unscaled values (mean > 500 for temperature) and applies
#     the correct scale factor.
#
# Land mask:
#   Based on GMTED2010 mean elevation at 30 arcsec. Cells with elevation
#   <= 0 are excluded. This removes ocean and some below-sea-level land
#   (Netherlands, Dead Sea, etc.) — acceptable for the global 1 km campaign
#   but should be noted.
#
# Part of: lpj-guess-ntd-vesta preprocessing pipeline
# ============================================================================

library(terra)
library(stringr)
library(data.table)
library(RNetCDF)

rm(list = ls())
graphics.off()

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG — edit these paths for your system
# ─────────────────────────────────────────────────────────────────────────────
mainDir   <- "/home/mdantas/chelsa_processing/global_v2"
chelsaV2  <- "/home/mdantas/chelsa_processing/chelsa_v2.1/monthly"
gmtedFile <- "/home/mdantas/map_data/GMTED2010_mean_30arcsec.nc4"

# Latitude extent: exclude Antarctica
no_antarctica_ext <- ext(-180, 180, -60, 90)

# Simulation period
years          <- 1979:2021
nyears         <- length(years)
ntime          <- nyears * 12   # 516 (includes Jan 1979 from v1 substitute)

# ─────────────────────────────────────────────────────────────────────────────
# LAND MASK FROM GMTED2010
# ─────────────────────────────────────────────────────────────────────────────
# GMTED2010 mean elevation at 30 arcsec (~1 km).
# Download from: https://topotools.cr.usgs.gov/gmted_viewer/viewer.htm
# or:  https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Grid_ZipFiles/mn30_grd.zip
# Convert to NetCDF if needed: gdal_translate -of netCDF mn30_grd GMTED2010_mean_30arcsec.nc4
cat("Loading GMTED2010 land mask...\n")
gmted <- rast(gmtedFile)
gmted_crop  <- crop(gmted, no_antarctica_ext)
land_mask   <- gmted_crop
land_mask[land_mask <= 0] <- NA      # ocean + below-sea-level (see header note)
land_mask[!is.na(land_mask)] <- 1
cat("  Land mask ready.\n")

# ─────────────────────────────────────────────────────────────────────────────
# TIME AXIS (365-day calendar, no leap years)
# ─────────────────────────────────────────────────────────────────────────────
days_per_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
cum_start      <- cumsum(c(0, days_per_month[-12]))

daysInMonth <- numeric(ntime)
for (j in seq_along(years))
  for (m in 1:12)
    daysInMonth[(j - 1L) * 12L + m] <- (j - 1L) * 365L + cum_start[m] + 1L

# ─────────────────────────────────────────────────────────────────────────────
# LATITUDE BANDS
# ─────────────────────────────────────────────────────────────────────────────
# Each band is processed independently (download → preprocess → run → gather
# → cleanup → next band), driven by storage constraints on HPC systems.
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

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

#' Create a NetCDF file for one climate variable and one latitude band.
#'
#' Dimension layout: grid_ID (UNLIMITED), time (fixed).
#' Variable layout:  var(grid_ID, time) — in ncdump notation.
#'
#' In RNetCDF, dimensions are specified in REVERSED order relative to ncdump:
#'   c("time", "grid_ID") in R  →  var(grid_ID, time) in the file.
#' This makes each cell's full timeseries contiguous on disk, which is what
#' LPJ-GUESS reads per MPI rank — much faster than time-contiguous layout.
open_lpj_nc <- function(filename, varname, standard_name, units_str,
                        ngrids, daysInMonth, years) {
  nc <- create.nc(filename, format = "netcdf4")

  # grid_ID = UNLIMITED (record dimension) — correct for GUESS I/O
  # time    = FIXED (length known upfront)
  dim.def.nc(nc, "grid_ID", unlim = TRUE)
  dim.def.nc(nc, "time",    dimlength = length(daysInMonth))

  var.def.nc(nc, "grid_ID", "NC_INT",    "grid_ID")
  var.def.nc(nc, "lon",     "NC_FLOAT",  "grid_ID")
  var.def.nc(nc, "lat",     "NC_FLOAT",  "grid_ID")
  var.def.nc(nc, "time",    "NC_DOUBLE", "time")

  att.put.nc(nc, "lon",  "standard_name", "NC_CHAR", "longitude")
  att.put.nc(nc, "lon",  "units",         "NC_CHAR", "degrees_east")
  att.put.nc(nc, "lat",  "standard_name", "NC_CHAR", "latitude")
  att.put.nc(nc, "lat",  "units",         "NC_CHAR", "degrees_north")
  att.put.nc(nc, "time", "units",         "NC_CHAR",
             paste0("days since ", years[1], "-01-01"))
  att.put.nc(nc, "time", "calendar",      "NC_CHAR", "365_day")

  # RNetCDF reverses dimension order vs ncdump:
  # c("time", "grid_ID") in R → var(grid_ID, time) in file
  var.def.nc(nc, varname, "NC_FLOAT", c("time", "grid_ID"))
  att.put.nc(nc, varname, "standard_name", "NC_CHAR", standard_name)
  att.put.nc(nc, varname, "units",         "NC_CHAR", units_str)
  att.put.nc(nc, varname, "_FillValue",    "NC_FLOAT", 1e+20)
  att.put.nc(nc, varname, "coordinates",   "NC_CHAR", "lat lon")

  var.put.nc(nc, "time",    daysInMonth)
  var.put.nc(nc, "grid_ID", 0L:(ngrids - 1L))

  nc
}

write_nc_coords <- function(nc, lon.grids, lat.grids) {
  var.put.nc(nc, "lon", lon.grids)
  var.put.nc(nc, "lat", lat.grids)
}

write_nc_timestep <- function(nc, varname, values, tidx) {
  var.put.nc(nc, varname, values,
             start = c(tidx, 1L), count = c(1L, length(values)))
}


# ─────────────────────────────────────────────────────────────────────────────
# FILE LISTS
# ─────────────────────────────────────────────────────────────────────────────
# CHELSA v2.1 naming: CHELSA_tas_MM_YYYY_V.2.1.tif
#                      CHELSA_pr_MM_YYYY_V.2.1.tif
#                      CHELSA_rsds_MM_YYYY_V.2.1.tif
# January 1979 (v1):  CHELSA_tas_01_1979_V1.2_land.tif  (or similar)

sort_chelsa_files <- function(dir, varname) {
  files <- list.files(dir, full.names = TRUE, pattern = "\\.tif$")
  years  <- as.integer(str_extract(basename(files), "\\d{4}"))
  months <- as.integer(str_extract(basename(files),
                                   paste0("(?<=", varname, "_)\\d{2}")))
  files[order(years, months)]
}

temp_files <- sort_chelsa_files(file.path(chelsaV2, "tas"),  "tas")
prec_files <- sort_chelsa_files(file.path(chelsaV2, "pr"),   "pr")
rsds_files <- sort_chelsa_files(file.path(chelsaV2, "rsds"), "rsds")

cat("Temperature files:", length(temp_files), "(expected", ntime, ")\n")
cat("Precipitation files:", length(prec_files), "(expected", ntime, ")\n")
cat("Radiation files:", length(rsds_files), "(expected", ntime, ")\n")

stopifnot(length(temp_files) == ntime,
          length(prec_files) == ntime,
          length(rsds_files) == ntime)

cat("First 15 temp files:\n")
cat(basename(head(temp_files, 15)), sep = "\n")

setwd(mainDir)


# ─────────────────────────────────────────────────────────────────────────────
# OUTER LOOP OVER LATITUDE BANDS
# ─────────────────────────────────────────────────────────────────────────────
for (band in band_defs) {

  band_name <- band[1]
  lat_min   <- as.numeric(band[2])
  lat_max   <- as.numeric(band[3])

  cat("\n", strrep("=", 60), "\n")
  cat(" Band:", band_name, " [", lat_min, "to", lat_max, "]\n")
  cat(strrep("=", 60), "\n")

  bandDir <- file.path(mainDir, band_name)

  # ── Resume: skip bands whose output files already exist ──────────────────
  expected_files <- c(
    file.path(bandDir, paste0("temp_", band_name, ".nc")),
    file.path(bandDir, paste0("prec_", band_name, ".nc")),
    file.path(bandDir, paste0("rad_",  band_name, ".nc")),
    file.path(bandDir, paste0("gridlist_", band_name, ".txt"))
  )
  if (all(file.exists(expected_files))) {
    sizes <- file.size(expected_files[1:3])
    if (all(sizes > 1e6)) {   # at least 1 MB each
      cat("  SKIPPING — output files already exist.\n")
      cat("  Delete them to reprocess this band.\n")
      next
    }
  }

  dir.create(bandDir, showWarnings = FALSE)
  band_ext <- ext(-180, 180, lat_min, lat_max)

  # ── Determine grid structure from first temp file ──────────────────────────
  r0             <- crop(rast(temp_files[1]), band_ext)
  land_mask_band <- crop(land_mask, band_ext)
  r_crop         <- crop(r0, no_antarctica_ext)
  r_land         <- mask(r_crop, land_mask_band)

  df0       <- as.data.frame(r_land, xy = TRUE, na.rm = TRUE)
  ngrids    <- nrow(df0)
  lon.grids <- df0$x
  lat.grids <- df0$y
  rm(df0)

  cat("  Grid cells:", ngrids, "\n")
  cat("  Estimated NC size per variable:",
      round(ngrids * ntime * 4 / 1e9, 2), "GB\n")

  # ── GRIDLIST ────────────────────────────────────────────────────────────────
  grid_ID  <- 0L:(ngrids - 1L)
  gridlist <- data.table(grid_ID,
                         paste("Cell", seq_len(ngrids), ngrids, sep = "-"))
  fwrite(gridlist,
         file.path(bandDir, paste0("gridlist_", band_name, ".txt")),
         sep = "\t", quote = FALSE, col.names = FALSE)

  # ── TEMPERATURE ─────────────────────────────────────────────────────────────
  # CHELSA v2.1: Scale=0.1, units=K → terra applies Scale → values in K
  # CHELSA v1.2 substitute (Jan 1979): terra may NOT apply Scale → detect
  #   by mean > 500 and divide by 10
  # LPJ-GUESS expects: temperature in K (converted internally)
  cat("\n  Pass 1: Temperature\n")
  nc_temp <- open_lpj_nc(
    file.path(bandDir, paste0("temp_", band_name, ".nc")),
    "temp", "air_temperature", "K", ngrids, daysInMonth, years)
  write_nc_coords(nc_temp, lon.grids, lat.grids)

  for (k in seq_along(temp_files)) {
    if (k %% 50 == 0 || k == 1 || k == ntime)
      cat("\r    File", k, "/", ntime, " ")
    flush.console()

    r      <- crop(rast(temp_files[k]), band_ext)
    r_crop <- crop(r, no_antarctica_ext)
    r_land <- mask(r_crop, land_mask_band)
    df     <- as.data.frame(r_land, xy = TRUE, na.rm = TRUE)
    vals   <- df[[3]]

    # Detect CHELSA v1 substitute (unscaled): values ~2700 instead of ~270
    if (mean(vals, na.rm = TRUE) > 500) {
      vals <- vals / 10
    }

    write_nc_timestep(nc_temp, "temp", vals, k)
  }
  close.nc(nc_temp)
  cat("\n  Done.\n")

  # ── PRECIPITATION ───────────────────────────────────────────────────────────
  # CHELSA v2.1: no scale, units = kg m-2 month-1 → use directly
  # LPJ-GUESS expects: precipitation in kg m-2 month-1
  cat("\n  Pass 2: Precipitation\n")
  nc_prec <- open_lpj_nc(
    file.path(bandDir, paste0("prec_", band_name, ".nc")),
    "prec", "precipitation_amount", "kg m-2", ngrids, daysInMonth, years)
  write_nc_coords(nc_prec, lon.grids, lat.grids)

  for (k in seq_along(prec_files)) {
    if (k %% 50 == 0 || k == 1 || k == ntime)
      cat("\r    File", k, "/", ntime, " ")
    flush.console()

    r      <- crop(rast(prec_files[k]), band_ext)
    r_crop <- crop(r, no_antarctica_ext)
    r_land <- mask(r_crop, land_mask_band)
    df     <- as.data.frame(r_land, xy = TRUE, na.rm = TRUE)
    vals   <- df[[3]]
    write_nc_timestep(nc_prec, "prec", vals, k)
  }
  close.nc(nc_prec)
  cat("\n  Done.\n")

  # ── SOLAR RADIATION ─────────────────────────────────────────────────────────
  # CHELSA v2.1: Scale=0.01157407407 (=1/86.4)
  # terra applies the Scale automatically on read.
  # Values are stored as W m-2 — matching LPJ-GUESS expectations.
  #
  # This is a TIMESERIES (one file per month per year), not a 12-month
  # climatology — same loop structure as temperature and precipitation.
  cat("\n  Pass 3: Solar radiation (timeseries)\n")
  nc_srad <- open_lpj_nc(
    file.path(bandDir, paste0("rad_", band_name, ".nc")),
    "insol", "surface_downwelling_shortwave_flux", "W m-2",
    ngrids, daysInMonth, years)
  write_nc_coords(nc_srad, lon.grids, lat.grids)

  for (k in seq_along(rsds_files)) {
    if (k %% 50 == 0 || k == 1 || k == ntime)
      cat("\r    File", k, "/", ntime, " ")
    flush.console()

    r      <- crop(rast(rsds_files[k]), band_ext)
    r_crop <- crop(r, no_antarctica_ext)
    r_land <- mask(r_crop, land_mask_band)
    df     <- as.data.frame(r_land, xy = TRUE, na.rm = TRUE)
    vals   <- df[[3]]
    write_nc_timestep(nc_srad, "insol", vals, k)
  }
  close.nc(nc_srad)
  cat("\n  Done.\n")

  cat("\n  Band", band_name, "complete.\n")
  gc()
}

cat("\nAll bands complete.\n")
