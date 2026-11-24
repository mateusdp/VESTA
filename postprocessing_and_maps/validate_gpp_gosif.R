# =============================================================================
# validate_gpp_gosif.R — Compare LPJ-GUESS AGPP against GOSIF GPP at 0.5°
# =============================================================================
#
# Validates modeled annual GPP against the GOSIF satellite product (v4.0).
# Both datasets are averaged over the temporal overlap period and compared
# at 0.5° resolution.
#
# Units:
#   GOSIF GPP v4.0: g C m-2 month-1 (monthly totals)
#     → Annual sum of 12 months / 1000 → kg C m-2 yr-1
#   LPJ-GUESS AGPP: kg C m-2 yr-1
#
# Output:
#   - Rasterized model TIF at 30 arc-sec
#   - Scatter plot (colored by latitude band)
#   - Density scatter (hex bins, log-scale)
#   - Spatial bias map (model − reference)
#   - Zonal mean comparison
#   - Summary metrics table (global + per latitude band)
#
# Requires: terra, tidyterra, ggplot2, data.table, ncdf4
# =============================================================================

library(terra)
library(tidyterra)
library(ggplot2)
library(data.table)
library(ncdf4)

# ---- CONFIG -----------------------------------------------------------------

base_dir   <- "/path/to/global_v2"          # Parent directory containing band folders
ref_file   <- "/path/to/GOSIF_GPP_monthly_0.5deg_v4.0.nc"
out_dir    <- "."                            # Output directory for plots

variable   <- "agpp"                         # LPJ-GUESS output file stem
col_name   <- "Total"                        # Column in .out.gz
var_name   <- "AGPP"                         # Display name
ref_name   <- "GOSIF GPP"                    # Reference dataset name
unit_label <- "kg C m⁻² yr⁻¹"

# Temporal overlap between LPJ-GUESS and GOSIF
year_start <- 2001
year_end   <- 2021

bands <- c("60S_30S", "30S_20S", "20S_10S", "10S_00N", "00N_10N",
           "10N_20N", "20N_30N", "30N_40N", "40N_50N", "50N_60N", "60N_90N")

# =============================================================================
# 1. CHECK GOSIF METADATA
# =============================================================================

cat("Checking GOSIF GPP metadata...\n")
nc <- nc_open(ref_file)
cat("  Variable:", names(nc$var), "\n")
for (att in c("units", "long_name", "scale_factor", "add_offset")) {
  val <- tryCatch(ncatt_get(nc, names(nc$var)[1], att), error = function(e) NULL)
  if (!is.null(val) && val$hasatt) cat(sprintf("  %s: %s\n", att, val$value))
}
time_vals <- ncvar_get(nc, "time")
origin <- as.Date("1970-01-01")
cat("  Time range:", format(origin + time_vals[1]), "to",
    format(origin + time_vals[length(time_vals)]), "\n")
cat("  Timesteps:", length(time_vals), "\n")
nc_close(nc)

# =============================================================================
# 2. LOAD GOSIF GPP AND COMPUTE ANNUAL MEAN
# =============================================================================

cat("\nLoading GOSIF GPP...\n")
ref_raw <- rast(ref_file)

# GOSIF v4.0 starts 2001-01. Determine month indices for overlap period.
idx_start <- (year_start - 2001) * 12 + 1
idx_end   <- (year_end   - 2001) * 12 + 12
cat("  Using layers", idx_start, "to", idx_end, "\n")
ref_subset <- subset(ref_raw, idx_start:idx_end)

# Convert monthly totals (g C m-2 month-1) to annual GPP (kg C m-2 yr-1)
n_years <- year_end - year_start + 1
annual_list <- list()
for (y in 1:n_years) {
  month_layers <- ((y - 1) * 12 + 1):(y * 12)
  annual_list[[y]] <- sum(subset(ref_subset, month_layers)) / 1000
}

ref_annual <- rast(annual_list)
ref_mean   <- mean(ref_annual, na.rm = TRUE)

# Clamp small negatives to zero (GOSIF artefact in low-productivity areas)
ref_mean <- clamp(ref_mean, lower = 0)

cat("  GOSIF mean GPP:", sprintf("%.3f", global(ref_mean, "mean", na.rm = TRUE)$mean),
    unit_label, "\n")

# =============================================================================
# 3. LOAD LPJ-GUESS AGPP
# =============================================================================

cat("\nLoading LPJ-GUESS AGPP...\n")

model_mean <- rbindlist(lapply(bands, function(b) {
  gz_file <- file.path(base_dir, b, paste0(variable, ".out.gz"))
  if (!file.exists(gz_file)) {
    cat("  SKIP:", b, "\n")
    return(NULL)
  }
  cat("  Reading:", b, "...")
  dt <- fread(cmd = paste("zcat", gz_file,
                          "| awk 'NR==1 || ($3 >=", year_start,
                          "&& $3 <=", year_end, ")'"))
  cat("", nrow(dt), "rows, averaging...")
  result <- dt[, .(value = mean(get(col_name), na.rm = TRUE)), by = .(Lon, Lat)]
  cat("", nrow(result), "cells\n")
  rm(dt); gc()
  result
}))

cat("  Total cells:", nrow(model_mean), "\n")
cat("  LPJ-GUESS mean AGPP:", sprintf("%.3f", mean(model_mean$value, na.rm = TRUE)),
    unit_label, "\n")

# Save model as 30 arc-sec GeoTIFF
template <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                 res = 1/120, crs = "EPSG:4326")
pts <- vect(model_mean, geom = c("Lon", "Lat"), crs = "EPSG:4326")
r   <- rasterize(pts, template, field = "value", fun = mean)

tif_file <- file.path(out_dir, paste0(variable, "_global_", year_start, "_",
                                       year_end, "_mean.tif"))
writeRaster(r, tif_file, overwrite = TRUE)
cat("  Saved:", tif_file, "\n")

# =============================================================================
# 4. AGGREGATE MODEL TO 0.5° AND MATCH GRIDS
# =============================================================================

cat("\nAggregating model to 0.5° grid...\n")

# Round to 0.5° grid centers
model_mean[, Lon_05 := round(Lon * 2) / 2 + 0.25]
model_mean[, Lat_05 := round(Lat * 2) / 2 + 0.25]
model_05 <- model_mean[, .(value = mean(value, na.rm = TRUE)),
                       by = .(Lon = Lon_05, Lat = Lat_05)]
cat("  0.5° cells:", nrow(model_05), "\n")

# Rasterize at 0.5°
template_05 <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                    res = 0.5, crs = "EPSG:4326")
pts_05     <- vect(model_05, geom = c("Lon", "Lat"), crs = "EPSG:4326")
model_agg  <- rasterize(pts_05, template_05, field = "value", fun = mean)

# Resample to exact GOSIF grid
model_agg <- resample(model_agg, ref_mean, method = "bilinear")

# Mask to common coverage
mask_common  <- !is.na(model_agg) & !is.na(ref_mean)
model_masked <- mask(model_agg, mask_common, maskvalue = FALSE)
ref_masked   <- mask(ref_mean,  mask_common, maskvalue = FALSE)

# =============================================================================
# 5. EXTRACT PAIRED VALUES AND COMPUTE METRICS
# =============================================================================

model_vals <- values(model_masked, na.rm = TRUE)
ref_vals   <- values(ref_masked,   na.rm = TRUE)

valid      <- is.finite(model_vals) & is.finite(ref_vals) & ref_vals >= 0
model_vals <- model_vals[valid]
ref_vals   <- ref_vals[valid]

coords <- xyFromCell(ref_masked, which(!is.na(values(ref_masked))))
lats   <- coords[valid, 2]

cat("\nPaired grid cells:", length(model_vals), "\n")

# ---- Metrics function ----
compute_metrics <- function(obs, mod, name = "") {
  residuals <- mod - obs
  metrics <- list(
    n         = length(obs),
    obs_mean  = mean(obs),
    mod_mean  = mean(mod),
    bias      = mean(residuals),
    bias_pct  = 100 * mean(residuals) / mean(obs),
    mae       = mean(abs(residuals)),
    rmse      = sqrt(mean(residuals^2)),
    nrmse     = sqrt(mean(residuals^2)) / sd(obs),
    r         = cor(obs, mod),
    r2        = cor(obs, mod)^2,
    mef       = 1 - sum(residuals^2) / sum((obs - mean(obs))^2),
    slope     = coef(lm(mod ~ obs))[2],
    intercept = coef(lm(mod ~ obs))[1]
  )

  cat(sprintf("\n==== %s ====\n", name))
  cat(sprintf("  N cells:     %d\n",          metrics$n))
  cat(sprintf("  Obs mean:    %.3f %s\n",     metrics$obs_mean, unit_label))
  cat(sprintf("  Model mean:  %.3f %s\n",     metrics$mod_mean, unit_label))
  cat(sprintf("  Bias:        %.3f (%.1f%%)\n", metrics$bias, metrics$bias_pct))
  cat(sprintf("  MAE:         %.3f\n",         metrics$mae))
  cat(sprintf("  RMSE:        %.3f\n",         metrics$rmse))
  cat(sprintf("  NRMSE:       %.2f\n",         metrics$nrmse))
  cat(sprintf("  r:           %.3f\n",         metrics$r))
  cat(sprintf("  R²:          %.3f\n",         metrics$r2))
  cat(sprintf("  MEF:         %.3f\n",         metrics$mef))
  cat(sprintf("  Slope:       %.3f\n",         metrics$slope))
  cat(sprintf("  Intercept:   %.3f\n",         metrics$intercept))

  metrics
}

# Global metrics
metrics_global <- compute_metrics(ref_vals, model_vals, "Global")

# Per latitude band
lat_bands <- list(
  "Tropics (23S-23N)"   = abs(lats) <= 23,
  "Subtropics (23-35)"  = abs(lats) > 23 & abs(lats) <= 35,
  "Temperate (35-55)"   = abs(lats) > 35 & abs(lats) <= 55,
  "Boreal (55-70)"      = abs(lats) > 55 & abs(lats) <= 70,
  "Arctic (>70)"        = abs(lats) > 70
)

metrics_bands <- list()
for (band_name in names(lat_bands)) {
  idx <- lat_bands[[band_name]]
  if (sum(idx) > 10) {
    metrics_bands[[band_name]] <- compute_metrics(ref_vals[idx], model_vals[idx], band_name)
  }
}

# =============================================================================
# 6. PLOTS
# =============================================================================

cat("\nGenerating plots...\n")

dt <- data.table(obs = ref_vals, mod = model_vals, lat = lats)
dt[, lat_band := cut(abs(lat), breaks = c(0, 23, 35, 55, 70, 90),
                     labels = c("Tropics", "Subtropics", "Temperate", "Boreal", "Arctic"))]

max_val <- quantile(c(dt$obs, dt$mod), 0.99, na.rm = TRUE)

# ---- Scatter plot (colored by latitude band) ----
if (nrow(dt) > 50000) {
  dt_plot <- dt[sample(.N, 50000)]
} else {
  dt_plot <- dt
}

p_scatter <- ggplot(dt_plot, aes(x = obs, y = mod, color = lat_band)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
  scale_color_brewer(palette = "Set1", name = "Latitude band") +
  labs(
    title = paste("LPJ-GUESS", var_name, "vs", ref_name,
                  paste0("(", year_start, "-", year_end, ")")),
    subtitle = sprintf("R² = %.3f | RMSE = %.3f | Bias = %.3f (%.1f%%) | MEF = %.3f",
                       metrics_global$r2, metrics_global$rmse,
                       metrics_global$bias, metrics_global$bias_pct, metrics_global$mef),
    x = paste(ref_name, "(", unit_label, ")"),
    y = paste("LPJ-GUESS", var_name, "(", unit_label, ")")
  ) +
  coord_fixed(xlim = c(0, max_val), ylim = c(0, max_val)) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(out_dir, "validation_gpp_scatter.png"),
       p_scatter, width = 9, height = 8, dpi = 300)

# ---- Density scatter (hex bins) ----
p_density <- ggplot(dt, aes(x = obs, y = mod)) +
  geom_hex(bins = 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  labs(
    title = paste("LPJ-GUESS", var_name, "vs", ref_name, "(density)"),
    subtitle = sprintf("R² = %.3f | RMSE = %.3f | Bias = %.3f (%.1f%%)",
                       metrics_global$r2, metrics_global$rmse,
                       metrics_global$bias, metrics_global$bias_pct),
    x = paste(ref_name, "(", unit_label, ")"),
    y = paste("LPJ-GUESS", var_name, "(", unit_label, ")")
  ) +
  coord_fixed(xlim = c(0, max_val), ylim = c(0, max_val)) +
  theme_minimal()

ggsave(file.path(out_dir, "validation_gpp_density.png"),
       p_density, width = 9, height = 8, dpi = 300)

# ---- Spatial bias map ----
bias_raster <- model_masked - ref_masked
bias_lim    <- quantile(abs(values(bias_raster, na.rm = TRUE)), 0.95)

p_bias <- ggplot() +
  geom_spatraster(data = bias_raster) +
  scale_fill_gradient2(
    name     = paste("Bias\n(", unit_label, ")"),
    low      = "brown", mid = "white", high = "darkgreen",
    midpoint = 0,
    limits   = c(-bias_lim, bias_lim),
    oob      = scales::squish,
    na.value = "transparent"
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 90)) +
  labs(title = paste(var_name, "bias: LPJ-GUESS minus", ref_name,
                     paste0("(", year_start, "-", year_end, ")"))) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))

ggsave(file.path(out_dir, "validation_gpp_bias_map.png"),
       p_bias, width = 16, height = 8, dpi = 300)

# ---- Zonal mean comparison ----
dt_zonal <- data.table(obs = ref_vals, mod = model_vals, lat = lats)
dt_zonal[, lat_bin := round(lat * 2) / 2]  # 0.5° bins
zonal <- dt_zonal[, .(GOSIF = mean(obs, na.rm = TRUE),
                      LPJ_GUESS = mean(mod, na.rm = TRUE),
                      n = .N), by = lat_bin]
zonal <- zonal[n >= 10]
zonal_long <- melt(zonal[, .(lat_bin, GOSIF, LPJ_GUESS)],
                   id.vars = "lat_bin", variable.name = "source", value.name = "GPP")

p_zonal <- ggplot(zonal_long, aes(x = GPP, y = lat_bin, color = source)) +
  geom_path(linewidth = 0.8) +
  scale_color_manual(values = c("GOSIF" = "blue", "LPJ_GUESS" = "red"), name = "") +
  labs(title = paste("Zonal mean", var_name, paste0("(", year_start, "-", year_end, ")")),
       x = unit_label, y = "Latitude") +
  theme_minimal()

ggsave(file.path(out_dir, "validation_gpp_zonal.png"),
       p_zonal, width = 7, height = 10, dpi = 300)

# =============================================================================
# 7. SUMMARY TABLE
# =============================================================================

cat("\n\n========== SUMMARY TABLE ==========\n")
cat(sprintf("%-25s %8s %8s %8s %8s %8s %8s\n",
            "Region", "R²", "RMSE", "Bias", "Bias%", "MAE", "MEF"))
cat(sprintf("%-25s %8.3f %8.3f %8.3f %7.1f%% %8.3f %8.3f\n",
            "Global", metrics_global$r2, metrics_global$rmse,
            metrics_global$bias, metrics_global$bias_pct,
            metrics_global$mae, metrics_global$mef))
for (bn in names(metrics_bands)) {
  m <- metrics_bands[[bn]]
  cat(sprintf("%-25s %8.3f %8.3f %8.3f %7.1f%% %8.3f %8.3f\n",
              bn, m$r2, m$rmse, m$bias, m$bias_pct, m$mae, m$mef))
}
cat("===================================\n")

cat("\nValidation complete.\n")
cat("Saved: validation_gpp_scatter.png, validation_gpp_density.png,",
    "validation_gpp_bias_map.png, validation_gpp_zonal.png\n")
