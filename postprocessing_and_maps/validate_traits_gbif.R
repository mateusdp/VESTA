# =============================================================================
# validate_traits_gbif.R — CWM traits: LPJ-GUESS vs GBIF-TRY gapfilled maps
# =============================================================================
#
# Validates modeled community-weighted mean (CWM) traits against gridded
# trait maps derived from GBIF occurrences and TRY gapfilled trait values
# (Schiller et al. / Kattge et al.).
#
# Traits validated: SLA, WSG, SRL
# Resolution: 0.5° (model aggregated from 30 arc-sec)
#
# Unit conversions (model → observation units):
#   SLA: m² kgC⁻¹ × 0.47 → m² kg⁻¹ DM   (carbon fraction ~0.47)
#   WSG: kgC m⁻³  × 1/470 → g cm⁻³        (carbon fraction ~0.47, unit scaling)
#   SRL: model units × 0.47/1000            (carbon fraction + g→kg)
#        GBIF reference: cm g⁻¹ × 1/100 → m g⁻¹
#
# Reference data format: .grd raster stacks with a "mean" layer.
#
# Output per trait:
#   - Scatter plot (base R, colored by latitude)
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

base_dir <- "/path/to/global_v2"             # Parent directory containing band folders
out_dir  <- "validation_traits"
dir.create(out_dir, showWarnings = FALSE)

# GBIF-TRY gapfilled reference rasters (0.5° .grd files with "mean" layer)
ref_files <- list(
  SLA = "/path/to/GBIF_TRYgapfilled_X11_05deg.grd",
  WSG = "/path/to/GBIF_TRYgapfilled_X4_05deg.grd",
  SRL = "/path/to/GBIF_TRYgapfilled_X1080_05deg.grd"
)

# Trait configuration: model output file, column, unit conversion
trait_config <- list(
  SLA = list(
    file       = "sla.out.gz",         # CWM across all PFTs
    col        = "Total",
    unit       = "m² kg⁻¹",
    conversion = 0.47                  # m²/kgC → m²/kg DM
  ),
  WSG = list(
    file       = "wsg.out.gz",         # CWM (tree PFTs dominate)
    col        = "Total",
    unit       = "g cm⁻³",
    conversion = 1/470                 # kgC/m³ → g/cm³
  ),
  SRL = list(
    file           = "srl.out.gz",     # CWM across all PFTs
    col            = "Total",
    unit           = "m g⁻¹",
    conversion     = 0.47/1000,        # model units → m/g
    ref_conversion = 1/100             # GBIF cm/g → m/g
  )
)

year_start <- 2017
year_end   <- 2021

bands <- c("60S_30S", "30S_20S", "20S_10S", "10S_00N", "00N_10N",
           "10N_20N", "20N_30N", "30N_40N", "40N_50N", "50N_60N", "60N_90N")

# =============================================================================
# METRICS FUNCTION
# =============================================================================

compute_metrics <- function(obs, mod, name = "", unit = "") {
  valid <- is.finite(obs) & is.finite(mod) & obs > 0 & mod > 0
  obs <- obs[valid]
  mod <- mod[valid]
  if (length(obs) < 10) return(NULL)

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
  cat(sprintf("  N cells:     %d\n",          metrics$n))
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
# PROCESS EACH TRAIT
# =============================================================================

all_metrics <- list()

for (trait_name in names(trait_config)) {

  cfg <- trait_config[[trait_name]]

  cat("\n##############################################\n")
  cat("  TRAIT:", trait_name, "\n")
  cat("##############################################\n")

  # ---- Load reference ----
  cat("Loading reference:", ref_files[[trait_name]], "\n")
  ref_rast <- rast(ref_files[[trait_name]])
  ref_mean <- ref_rast[["mean"]]
  cat("  Reference range:", range(values(ref_mean), na.rm = TRUE), cfg$unit, "\n")

  # ---- Load model ----
  cat("Loading model output...\n")
  model_mean <- rbindlist(lapply(bands, function(b) {
    gz_file <- file.path(base_dir, b, cfg$file)
    if (!file.exists(gz_file)) {
      cat("  SKIP:", b, "\n")
      return(NULL)
    }
    cat("  Reading:", b, "...")
    dt <- fread(cmd = paste("zcat", gz_file,
                            "| awk 'NR==1 || ($3 >=", year_start,
                            "&& $3 <=", year_end, ")'"))
    cat("", nrow(dt), "rows, averaging...")
    result <- dt[, .(value = mean(get(cfg$col), na.rm = TRUE)), by = .(Lon, Lat)]
    cat("", nrow(result), "cells\n")
    rm(dt); gc()
    result
  }))

  # ---- Unit conversion ----
  if (cfg$conversion != 1) {
    cat(sprintf("  Applying conversion factor: × %.6f\n", cfg$conversion))
    model_mean[, value := value * cfg$conversion]
  }

  cat("  Model range:", range(model_mean$value, na.rm = TRUE), cfg$unit, "\n")
  cat("  Model mean:",  mean(model_mean$value, na.rm = TRUE), cfg$unit, "\n")

  # ---- Aggregate model to 0.5° ----
  cat("Aggregating to 0.5°...\n")
  model_mean[, Lon_05 := round(Lon * 2) / 2 + 0.25]
  model_mean[, Lat_05 := round(Lat * 2) / 2 + 0.25]
  model_05 <- model_mean[, .(value = mean(value, na.rm = TRUE)),
                         by = .(Lon = Lon_05, Lat = Lat_05)]
  cat("  0.5° cells:", nrow(model_05), "\n")

  # ---- Match grids ----
  ref_dt <- as.data.table(as.data.frame(ref_mean, xy = TRUE))
  names(ref_dt) <- c("Lon", "Lat", "ref_value")
  ref_dt <- ref_dt[!is.na(ref_value) & ref_value > 0]

  # Convert reference units if needed
  if (!is.null(cfg$ref_conversion)) {
    ref_dt[, ref_value := ref_value * cfg$ref_conversion]
  }

  ref_dt[,   Lon := round(Lon, 2)]
  ref_dt[,   Lat := round(Lat, 2)]
  model_05[, Lon := round(Lon, 2)]
  model_05[, Lat := round(Lat, 2)]

  paired <- merge(model_05, ref_dt, by = c("Lon", "Lat"))

  model_vals <- as.numeric(paired$value)
  ref_vals   <- as.numeric(paired$ref_value)
  lats       <- paired$Lat

  cat("  Paired cells:", length(model_vals), "\n")

  # Quick unit sanity check
  ratio <- mean(model_vals, na.rm = TRUE) / mean(ref_vals, na.rm = TRUE)
  cat(sprintf("  Model/Ref ratio: %.2f (should be ~1.0; if ~2 or ~0.5, check unit conversion)\n", ratio))

  # ---- Global metrics ----
  m_global <- compute_metrics(ref_vals, model_vals, paste(trait_name, "- Global"), cfg$unit)
  all_metrics[[trait_name]] <- m_global

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
      compute_metrics(ref_vals[idx], model_vals[idx],
                      paste(trait_name, "-", band_name), cfg$unit)
    }
  }

  # ---- Scatter plot (base R) ----
  vmax <- quantile(c(model_vals, ref_vals), 0.99, na.rm = TRUE)

  png(file.path(out_dir, paste0("trait_scatter_", trait_name, ".png")),
      width = 2400, height = 2400, res = 300)

  plot(ref_vals, model_vals,
       pch = ".", col = rgb(0, 0.3, 0.6, 0.05),
       xlim = c(0, vmax), ylim = c(0, vmax),
       xlab = paste("GBIF-TRY", trait_name, "(", cfg$unit, ")"),
       ylab = paste("LPJ-GUESS", trait_name, "(", cfg$unit, ")"),
       main = paste0(trait_name, ": LPJ-GUESS vs GBIF-TRY (",
                     year_start, "-", year_end, ")"),
       sub  = sprintf("R² = %.3f | RMSE = %.3f | Bias = %.1f%%",
                      m_global$r2, m_global$rmse, m_global$bias_pct),
       asp = 1)
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  abline(lm(model_vals ~ ref_vals), col = "black", lwd = 1)
  dev.off()

  # ---- Bias map ----
  template_05 <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                      res = 0.5, crs = "EPSG:4326")
  pts_mod <- vect(paired[, .(Lon, Lat, value)],     geom = c("Lon", "Lat"), crs = "EPSG:4326")
  pts_ref <- vect(paired[, .(Lon, Lat, ref_value)], geom = c("Lon", "Lat"), crs = "EPSG:4326")

  r_mod  <- rasterize(pts_mod, template_05, field = "value",     fun = mean)
  r_ref  <- rasterize(pts_ref, template_05, field = "ref_value", fun = mean)
  r_bias <- r_mod - r_ref

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
    labs(title = paste0(trait_name, " bias: LPJ-GUESS minus GBIF-TRY")) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

  ggsave(file.path(out_dir, paste0("trait_bias_", trait_name, ".png")),
         p_bias, width = 16, height = 8, dpi = 300)

  # ---- Zonal mean (base R) ----
  dt_zonal <- data.table(obs = ref_vals, mod = model_vals, lat = lats)
  dt_zonal[, lat_bin := round(lat)]
  zonal <- dt_zonal[, .(GBIF_TRY  = mean(obs, na.rm = TRUE),
                        LPJ_GUESS = mean(mod, na.rm = TRUE),
                        n = .N), by = lat_bin]
  zonal <- zonal[n >= 20]
  zonal <- zonal[order(lat_bin)]

  png(file.path(out_dir, paste0("trait_zonal_", trait_name, ".png")),
      width = 2100, height = 3000, res = 300)

  plot(zonal$GBIF_TRY, zonal$lat_bin,
       type = "l", col = "blue", lwd = 2,
       xlim = c(0, max(c(zonal$GBIF_TRY, zonal$LPJ_GUESS), na.rm = TRUE)),
       xlab = cfg$unit, ylab = "Latitude",
       main = paste("Zonal mean", trait_name))
  lines(zonal$LPJ_GUESS, zonal$lat_bin, col = "red", lwd = 2)
  legend("topright", legend = c("GBIF-TRY", "LPJ-GUESS"),
         col = c("blue", "red"), lwd = 2)
  dev.off()

  cat("\n  Plots saved for", trait_name, "\n")

  rm(model_mean, model_05, paired, ref_dt, ref_rast, ref_mean)
  gc()
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n\n========== TRAIT VALIDATION SUMMARY (GBIF-TRY) ==========\n")
cat(sprintf("%-8s %8s %8s %8s %8s %8s %8s\n",
            "Trait", "R²", "RMSE", "Bias%", "MAE", "MEF", "Ratio"))
cat("--------------------------------------------------------------\n")
for (trait_name in names(all_metrics)) {
  m <- all_metrics[[trait_name]]
  if (!is.null(m)) {
    cat(sprintf("%-8s %8.3f %8.3f %7.1f%% %8.3f %8.3f %8.2f\n",
                trait_name, m$r2, m$rmse, m$bias_pct, m$mae, m$mef,
                m$mod_mean / m$obs_mean))
  }
}
cat("===============================================\n")
cat("\nNOTE: If Model/Ref ratio is ~2 or ~0.5, check unit conversions.\n")
cat("  SLA: model uses m²/kgC, reference uses m²/kg dry mass\n")
cat("  WSG: both should be g/cm³ after conversion\n")
cat("  SRL: check model units vs reference units\n")

cat("\nAll trait validation plots saved to:", out_dir, "\n")
