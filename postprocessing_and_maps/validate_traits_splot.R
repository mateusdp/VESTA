# =============================================================================
# validate_traits_splot.R — CWM traits: LPJ-GUESS vs sPlotOpen
# =============================================================================
#
# Validates modeled community-weighted mean (CWM) traits against sPlotOpen
# vegetation-plot data (~95K plots worldwide).
#
# Traits validated: SLA, WSG, Leaf C:N, Leaf C:P
#
# Pipeline:
#   1. Load sPlotOpen header + CWM/CWV files, derive C:N and C:P from
#      ln-transformed LeafC, LeafN, LeafP
#   2. Read cwm_traits.out.gz using awk to select columns + filter years
#      (keeps pipe under 2 GB limit for large output files)
#   3. Cache 0.5° model rasters per trait for fast reuse
#   4. Extract model values at sPlotOpen plot locations
#   5. Compute metrics, produce scatter/bias/zonal plots
#
# Unit conversions (model → observation units):
#   SLA: m² kgC⁻¹  × 0.47  → m² kg⁻¹ DM
#   WSG: kgC m⁻³   × 1/470 → g cm⁻³
#   C:N: model C:N  × model_conv → g C / g N   (see note below)
#   C:P: model C:P  × model_conv → g C / g P   (see note below)
#
# NOTE on model_conv for C:N and C:P:
#   The model's daily growth routine may accumulate C:N over nday_leafon days.
#   If the averaging missed a factor of 2 (nday_leafon/2 correction), set
#   model_conv = 2. If the output already accounts for this, set model_conv = 1.
#   Check by comparing the model/observation ratio — it should be ~1.0.
#
# sPlotOpen data:
#   - 3474_76_sPlotOpen_header.txt     (plot coordinates and metadata)
#   - 3474_76_sPlotOpen_CWM_CWV.txt    (community-weighted means/variances)
#   - CWM traits are ln-transformed; C:N = exp(LeafC) / exp(LeafN)
#
# Output per trait:
#   - Scatter plot (base R, colored by latitude band)
#   - Spatial bias map (ggplot2, diverging)
#   - Zonal mean comparison (base R)
#   - Console metrics (global + per latitude band)
#
# Requires: terra, tidyterra, ggplot2, data.table
# =============================================================================

library(terra)
library(tidyterra)
library(ggplot2)
library(data.table)

# ---- CONFIG -----------------------------------------------------------------

base_dir  <- "/path/to/global_v2"           # Parent directory containing band folders
splot_dir <- "/path/to/splot"               # Directory with sPlotOpen files
out_dir   <- "validation_traits_splot"
dir.create(out_dir, showWarnings = FALSE)

year_start <- 2017
year_end   <- 2021

bands <- c("60S_30S", "30S_20S", "20S_10S", "10S_00N", "00N_10N",
           "10N_20N", "20N_30N", "30N_40N", "40N_50N", "50N_60N", "60N_90N")

# Trait configuration
# model_conv: conversion factor applied to model output AFTER extraction
# splot_ln:   TRUE if sPlotOpen stores the trait as ln(value)
trait_config <- list(
  SLA = list(
    model_col  = "SLA",
    model_conv = 0.47,                 # m²/kgC → m²/kg DM
    invert     = FALSE,
    splot_col  = "SLA_CWM",
    splot_ln   = TRUE,
    unit       = "m² kg⁻¹"
  ),
  WSG = list(
    model_col  = "WSG",
    model_conv = 1/470,                # kgC/m³ → g/cm³
    invert     = FALSE,
    splot_col  = "StemDens_CWM",
    splot_ln   = TRUE,
    unit       = "g cm⁻³"
  ),
  LeafCN = list(
    model_col  = "CN_leaf",
    model_conv = 2,                    # ×2 bugfix: daily growth nday_leafon/2
    invert     = FALSE,                # Set to 1 if model output is already corrected
    splot_col  = "CN_derived",
    splot_ln   = FALSE,
    unit       = "g C / g N"
  ),
  LeafCP = list(
    model_col  = "CP_leaf",
    model_conv = 2,                    # ×2 bugfix: daily growth nday_leafon/2
    invert     = FALSE,                # Set to 1 if model output is already corrected
    splot_col  = "CP_derived",
    splot_ln   = FALSE,
    unit       = "g C / g P"
  )
)

# =============================================================================
# 1. LOAD sPlotOpen
# =============================================================================

cat("Loading sPlotOpen...\n")

plots <- fread(file.path(splot_dir, "3474_76_sPlotOpen_header.txt"), sep = "\t")
cwm   <- fread(file.path(splot_dir, "3474_76_sPlotOpen_CWM_CWV.txt"), sep = "\t")
splot <- merge(cwm, plots, by = "PlotObservationID")

splot <- splot[!is.na(Longitude) & !is.na(Latitude)]
splot <- splot[Latitude >= -60 & Latitude <= 90]
cat("  Plots:", nrow(splot), "\n")

# Derive C:N and C:P from ln-transformed leaf nutrient CWMs
# C:N = exp(ln(LeafC)) / exp(ln(LeafN)) = LeafC / LeafN (mass basis)
splot[, CN_derived := exp(LeafC_perdrymass_CWM) / exp(LeafN_CWM)]
splot[, CP_derived := exp(LeafC_perdrymass_CWM) / exp(LeafP_CWM)]
cat("  Derived C:N range:", round(range(splot$CN_derived, na.rm = TRUE), 2), "\n")
cat("  Derived C:P range:", round(range(splot$CP_derived, na.rm = TRUE), 2), "\n")

# =============================================================================
# 2. BUILD AWK COMMAND FOR COLUMN SELECTION
# =============================================================================

# Identify which model columns we need
trait_cols <- unique(sapply(trait_config, function(cfg) cfg$model_col))
keep_cols  <- c("Lon", "Lat", "Year", trait_cols)

# awk script: find column indices from header, then print only needed columns.
# This cuts output from ~13+ columns to ~7, staying under the 2 GB pipe limit.
awk_header_print <- paste0('"', keep_cols, '"', collapse = ',')
awk_data_print   <- paste0('$c["', keep_cols, '"]', collapse = ',')

awk_template <- paste0(
  'BEGIN{OFS="\\t"} ',
  'NR==1{for(i=1;i<=NF;i++) c[$i]=i; print ', awk_header_print, '; next} ',
  '$c["Year"]==YEARPLACEHOLDER{print ', awk_data_print, '}'
)

cat("  AWK selects columns:", paste(keep_cols, collapse = ", "), "\n")

# =============================================================================
# 3. READ MODEL OUTPUT AND CACHE 0.5° RASTERS
# =============================================================================

all_cached <- all(sapply(names(trait_config), function(t)
  file.exists(file.path(out_dir, paste0("model_05deg_", t, ".tif")))))

if (!all_cached) {
  cat("\nReading cwm_traits.out.gz — year by year, columns selected by awk...\n\n")

  model_all <- rbindlist(lapply(bands, function(b) {
    gz_file <- file.path(base_dir, b, "cwm_traits.out.gz")
    if (!file.exists(gz_file)) return(NULL)
    cat("  Reading:", b, "...")

    yearly <- rbindlist(lapply(year_start:year_end, function(y) {
      awk_cmd <- gsub("YEARPLACEHOLDER", y, awk_template)
      cmd <- paste("zcat", gz_file, "| awk '", awk_cmd, "'")
      dt  <- fread(cmd = cmd)
      # Force trait columns to numeric (handles "nan", "-inf" in output)
      for (col in intersect(trait_cols, names(dt))) {
        if (!is.numeric(dt[[col]])) dt[, (col) := as.numeric(get(col))]
      }
      dt
    }))

    cat("", nrow(yearly), "rows...")
    result <- yearly[, lapply(.SD, mean, na.rm = TRUE),
                     by = .(Lon, Lat), .SDcols = trait_cols]
    cat("", nrow(result), "cells\n")
    rm(yearly); gc()
    result
  }))

  cat("\n  Total cells:", nrow(model_all), "\n")

  # Aggregate to 0.5°
  model_all[, Lon_05 := round(Lon * 2) / 2 + 0.25]
  model_all[, Lat_05 := round(Lat * 2) / 2 + 0.25]

  template_05 <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                      res = 0.5, crs = "EPSG:4326")

  for (trait_name in names(trait_config)) {
    cfg <- trait_config[[trait_name]]
    cat("  Caching:", trait_name, "...")

    tmp <- model_all[, .(Lon = Lon_05, Lat = Lat_05, value = get(cfg$model_col))]
    tmp[, value := value * cfg$model_conv]

    if (cfg$invert) {
      tmp[value > 0, value := cfg$model_conv / value]
      tmp[value == 0, value := NA]
    }

    # Remove zeros (WSG = 0 for grass-only cells, etc.)
    tmp <- tmp[value > 0 & is.finite(value)]
    tmp_05 <- tmp[, .(mod = mean(value, na.rm = TRUE)), by = .(Lon, Lat)]

    cat(" range:", round(range(tmp_05$mod, na.rm = TRUE), 3), cfg$unit, "...")

    pts <- vect(tmp_05, geom = c("Lon", "Lat"), crs = "EPSG:4326")
    r   <- rasterize(pts, template_05, field = "mod", fun = mean)

    cache_file <- file.path(out_dir, paste0("model_05deg_", trait_name, ".tif"))
    writeRaster(r, cache_file, datatype = "FLT4S", overwrite = TRUE,
                gdal = c("COMPRESS=DEFLATE", "ZLEVEL=9"))
    cat(" saved\n")
  }

  rm(model_all); gc()

} else {
  cat("\nAll trait rasters cached. Skipping .gz read.\n")
}

# =============================================================================
# 4. METRICS FUNCTION
# =============================================================================

compute_metrics <- function(obs, mod, name = "", unit = "") {
  valid <- is.finite(obs) & is.finite(mod) & obs > 0 & mod > 0
  obs <- obs[valid]
  mod <- mod[valid]
  if (length(obs) < 10) {
    cat(sprintf("\n==== %s ==== Too few valid pairs (%d)\n", name, length(obs)))
    return(NULL)
  }

  residuals <- mod - obs
  metrics <- list(
    n        = length(obs),
    obs_mean = mean(obs),
    mod_mean = mean(mod),
    bias     = mean(residuals),
    bias_pct = 100 * mean(residuals) / mean(obs),
    mae      = mean(abs(residuals)),
    rmse     = sqrt(mean(residuals^2)),
    nrmse    = sqrt(mean(residuals^2)) / sd(obs),
    r        = cor(obs, mod),
    r2       = cor(obs, mod)^2,
    mef      = 1 - sum(residuals^2) / sum((obs - mean(obs))^2)
  )

  cat(sprintf("\n==== %s ====\n", name))
  cat(sprintf("  N plots:     %d\n",          metrics$n))
  cat(sprintf("  Obs mean:    %.3f %s\n",     metrics$obs_mean, unit))
  cat(sprintf("  Model mean:  %.3f %s\n",     metrics$mod_mean, unit))
  cat(sprintf("  Bias:        %.3f (%.1f%%)\n", metrics$bias, metrics$bias_pct))
  cat(sprintf("  MAE:         %.3f\n",         metrics$mae))
  cat(sprintf("  RMSE:        %.3f\n",         metrics$rmse))
  cat(sprintf("  NRMSE:       %.2f\n",         metrics$nrmse))
  cat(sprintf("  r:           %.3f\n",         metrics$r))
  cat(sprintf("  R²:          %.3f\n",         metrics$r2))
  cat(sprintf("  MEF:         %.3f\n",         metrics$mef))

  metrics
}

# =============================================================================
# 5. VALIDATE EACH TRAIT
# =============================================================================

all_metrics  <- list()
template_05  <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                     res = 0.5, crs = "EPSG:4326")

for (trait_name in names(trait_config)) {

  cfg <- trait_config[[trait_name]]

  cat("\n##############################################\n")
  cat("  TRAIT:", trait_name, "\n")
  cat("##############################################\n")

  # Load cached model raster
  cache_file <- file.path(out_dir, paste0("model_05deg_", trait_name, ".tif"))
  model_rast <- rast(cache_file)
  cat("  Loaded cached raster\n")

  # Check sPlotOpen column exists
  if (!(cfg$splot_col %in% names(splot))) {
    cat("  WARNING: Column", cfg$splot_col, "not found. Skipping.\n")
    next
  }

  # Extract sPlotOpen observations
  splot_trait <- splot[!is.na(get(cfg$splot_col)),
                       .(PlotObservationID, Lon = Longitude, Lat = Latitude,
                         obs_raw = get(cfg$splot_col))]

  if (cfg$splot_ln) {
    splot_trait[, obs := exp(obs_raw)]
  } else {
    splot_trait[, obs := obs_raw]
  }

  splot_trait <- splot_trait[is.finite(obs) & obs > 0]
  cat("  sPlotOpen plots:", nrow(splot_trait), "\n")
  cat("  sPlotOpen range:", round(range(splot_trait$obs, na.rm = TRUE), 3), cfg$unit, "\n")

  # Extract model values at plot locations
  splot_pts <- vect(splot_trait, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  extracted <- as.numeric(extract(model_rast, splot_pts)[, 2])
  splot_trait[, mod := extracted]

  paired <- splot_trait[is.finite(obs) & is.finite(mod) & obs > 0 & mod > 0]
  cat("  Paired plots:", nrow(paired), "\n")

  obs_vals <- paired$obs
  mod_vals <- paired$mod
  lats     <- paired$Lat

  ratio <- mean(mod_vals, na.rm = TRUE) / mean(obs_vals, na.rm = TRUE)
  cat(sprintf("  Model/Obs ratio: %.2f\n", ratio))

  # ---- Global metrics ----
  m_global <- compute_metrics(obs_vals, mod_vals,
                              paste(trait_name, "- Global (sPlotOpen)"), cfg$unit)
  all_metrics[[trait_name]] <- m_global
  if (is.null(m_global)) next

  # ---- Per latitude band ----
  lat_bands <- list(
    "Tropics (23S-23N)"  = abs(lats) <= 23,
    "Subtropics (23-35)" = abs(lats) > 23 & abs(lats) <= 35,
    "Temperate (35-55)"  = abs(lats) > 35 & abs(lats) <= 55,
    "Boreal (55-70)"     = abs(lats) > 55 & abs(lats) <= 70
  )

  for (band_name in names(lat_bands)) {
    idx <- lat_bands[[band_name]]
    if (sum(idx) > 10) {
      compute_metrics(obs_vals[idx], mod_vals[idx],
                      paste(trait_name, "-", band_name), cfg$unit)
    }
  }

  # ---- Scatter plot (base R, colored by latitude) ----
  vmax <- quantile(c(obs_vals, mod_vals), 0.99, na.rm = TRUE)

  cols <- rep("grey50", length(lats))
  cols[abs(lats) <= 23]                        <- "red"
  cols[abs(lats) > 23 & abs(lats) <= 35]       <- "blue"
  cols[abs(lats) > 35 & abs(lats) <= 55]       <- "green3"
  cols[abs(lats) > 55 & abs(lats) <= 70]       <- "purple"
  cols[abs(lats) > 70]                          <- "orange"

  png(file.path(out_dir, paste0("splot_cwm_scatter_", trait_name, ".png")),
      width = 2400, height = 2400, res = 300)

  plot(obs_vals, mod_vals,
       pch = 16, cex = 0.3, col = adjustcolor(cols, alpha.f = 0.4),
       xlim = c(0, vmax), ylim = c(0, vmax),
       xlab = paste("sPlotOpen CWM", trait_name, "(", cfg$unit, ")"),
       ylab = paste("LPJ-GUESS CWM", trait_name, "(", cfg$unit, ")"),
       main = paste0(trait_name, ": LPJ-GUESS vs sPlotOpen"),
       sub  = sprintf("R² = %.3f | RMSE = %.3f | Bias = %.1f%% | n = %d",
                      m_global$r2, m_global$rmse, m_global$bias_pct, m_global$n),
       asp = 1)
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  abline(lm(mod_vals ~ obs_vals), col = "black", lwd = 1)
  legend("topleft",
         legend = c("Tropics", "Subtropics", "Temperate", "Boreal"),
         col = c("red", "blue", "green3", "purple"),
         pch = 16, cex = 0.8, bty = "n")
  dev.off()

  # ---- Bias map ----
  paired[, bias := mod - obs]
  bias_05 <- paired[, .(bias = mean(bias, na.rm = TRUE)),
                    by = .(Lon = round(Lon * 2) / 2 + 0.25,
                           Lat = round(Lat * 2) / 2 + 0.25)]

  pts_bias <- vect(bias_05, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  r_bias   <- rasterize(pts_bias, template_05, field = "bias", fun = mean)

  bias_lim <- quantile(abs(values(r_bias, na.rm = TRUE)), 0.95)

  p_bias <- ggplot() +
    geom_spatraster(data = r_bias) +
    scale_fill_gradient2(
      name     = paste("Bias\n(", cfg$unit, ")"),
      low      = "brown", mid = "white", high = "darkgreen",
      midpoint = 0,
      limits   = c(-bias_lim, bias_lim),
      oob      = scales::squish,
      na.value = "transparent"
    ) +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 90)) +
    labs(title = paste0(trait_name, " bias: LPJ-GUESS minus sPlotOpen")) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

  ggsave(file.path(out_dir, paste0("splot_cwm_bias_", trait_name, ".png")),
         p_bias, width = 16, height = 8, dpi = 300)

  # ---- Zonal mean (base R) ----
  dt_zonal <- data.table(obs = obs_vals, mod = mod_vals, lat = lats)
  dt_zonal[, lat_bin := round(lat / 5) * 5]
  zonal <- dt_zonal[, .(sPlotOpen  = mean(obs, na.rm = TRUE),
                        LPJ_GUESS = mean(mod, na.rm = TRUE),
                        n = .N), by = lat_bin]
  zonal <- zonal[n >= 5]
  zonal <- zonal[order(lat_bin)]

  png(file.path(out_dir, paste0("splot_cwm_zonal_", trait_name, ".png")),
      width = 2100, height = 3000, res = 300)

  xmax <- max(c(zonal$sPlotOpen, zonal$LPJ_GUESS), na.rm = TRUE) * 1.1
  plot(zonal$sPlotOpen, zonal$lat_bin,
       type = "l", col = "blue", lwd = 2,
       xlim = c(0, xmax),
       xlab = cfg$unit, ylab = "Latitude",
       main = paste("Zonal mean", trait_name, "(sPlotOpen vs LPJ-GUESS)"))
  lines(zonal$LPJ_GUESS, zonal$lat_bin, col = "red", lwd = 2)
  points(zonal$sPlotOpen, zonal$lat_bin, col = "blue", pch = 16, cex = 0.8)
  points(zonal$LPJ_GUESS, zonal$lat_bin, col = "red",  pch = 16, cex = 0.8)
  legend("topright", legend = c("sPlotOpen", "LPJ-GUESS"),
         col = c("blue", "red"), lwd = 2, pch = 16)
  dev.off()

  cat("\n  Plots saved for", trait_name, "\n")
  rm(paired, splot_trait); gc()
}

# =============================================================================
# 6. SUMMARY TABLE
# =============================================================================

cat("\n\n========== TRAIT VALIDATION SUMMARY (sPlotOpen) ==========\n")
cat(sprintf("%-10s %8s %10s %10s %8s %8s %8s\n",
            "Trait", "N plots", "Obs mean", "Mod mean", "Bias%", "R²", "MEF"))
cat("-------------------------------------------------------------------\n")
for (trait_name in names(all_metrics)) {
  m <- all_metrics[[trait_name]]
  if (!is.null(m)) {
    cat(sprintf("%-10s %8d %10.3f %10.3f %7.1f%% %8.3f %8.3f\n",
                trait_name, m$n, m$obs_mean, m$mod_mean,
                m$bias_pct, m$r2, m$mef))
  }
}
cat("===================================================================\n")

cat("\nCached rasters in:", out_dir, "\n")
for (t in names(trait_config)) {
  f <- file.path(out_dir, paste0("model_05deg_", t, ".tif"))
  if (file.exists(f)) cat("  ", f, "\n")
}

cat("\nDone.\n")
