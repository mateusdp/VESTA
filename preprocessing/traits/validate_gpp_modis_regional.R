# =============================================================================
# validate_gpp_modis_regional.R — Regional comparison: LPJ-GUESS vs MODIS GPP
#                                  at 1 km resolution
# =============================================================================
#
# Demonstrates the added value of 1 km resolution by comparing LPJ-GUESS AGPP
# against MODIS MOD17A3HGF (500 m, aggregated to 1 km) in selected regions
# with strong topographic or land-cover gradients.
#
# Prerequisites:
#   - LPJ-GUESS AGPP already saved as GeoTIFF at 30 arc-sec
#     (from plot_global.R workflow)
#   - MODIS MOD17A3HGF annual GPP downloaded from AppEEARS as GeoTIFF
#     (https://appeears.earthdatacloud.nasa.gov/)
#     Select: MOD17A3HGF.061, Gpp_500m, for each region
#
# MODIS units: kg C m-2 (AppEEARS GeoTIFFs typically have scale applied)
#
# Output per region:
#   - Side-by-side maps: LPJ-GUESS | MODIS | Difference (patchwork)
#   - Density scatter (hex bins)
#   - E-W transect at mid-latitude
#
# Requires: terra, tidyterra, ggplot2, data.table, patchwork
# =============================================================================

library(terra)
library(tidyterra)
library(ggplot2)
library(data.table)
library(patchwork)

# ---- CONFIG -----------------------------------------------------------------

# LPJ-GUESS model output (30 arc-sec GeoTIFF from plot_global.R)
model_tif <- "agpp_global_2017_2021_mean.tif"

# Directory containing MODIS GPP GeoTIFFs from AppEEARS.
# Organize as: modis_dir/<region_id>/MOD17A3HGF_*.tif
modis_dir <- "/path/to/validation/modis"

# Output directory
out_dir <- "validation_1km"
dir.create(out_dir, showWarnings = FALSE)

# Year range (must match what's in model_tif)
year_start <- 2017
year_end   <- 2021

# Regions to compare: uncomment those with available MODIS data.
# Each region: id (must match a subdirectory in modis_dir),
#              name, extent, description, plot dimensions.
regions <- list(
  andes = list(
    name = "Andes", ext = ext(-80, -75, -2, 4),
    desc = "Amazon lowland to Andean páramo",
    width = 8, height = 10
  ),
  alps = list(
    name = "European Alps", ext = ext(6, 14, 45, 48),
    desc = "Valley-to-peak gradients",
    width = 12, height = 6
  ),
  fjords = list(
    name = "Norwegian Fjords", ext = ext(4, 9, 59, 63),
    desc = "Coastal vs inland productivity",
    width = 8, height = 8
  ),
  borneo = list(
    name = "Borneo", ext = ext(108, 119, -4, 7),
    desc = "Tropical island heterogeneity",
    width = 10, height = 10
  ),
  se_australia = list(
    name = "SE Australia", ext = ext(143, 154, -39, -28),
    desc = "Great Dividing Range gradient",
    width = 8, height = 10
  ),
  amazon_arc = list(
    name = "Amazon Arc", ext = ext(-65, -50, -15, -5),
    desc = "Deforestation frontier",
    width = 10, height = 8
  )
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load MODIS MOD17A3HGF GeoTIFFs for one region, compute multi-year mean.
#'
#' Expects AppEEARS naming convention: *doyYYYYDDD* in filenames.
#' Clamps fill values (GPP > 0.5 kg C m-2 per 8-day composite) to NA.
load_modis_region <- function(modis_path, yr_start, yr_end) {
  tifs  <- list.files(modis_path, pattern = "\\.tif$", full.names = TRUE)
  if (length(tifs) == 0) return(NULL)

  years <- as.integer(gsub(".*doy(\\d{4})\\d{3}.*", "\\1", basename(tifs)))
  keep  <- years >= yr_start & years <= yr_end
  tifs  <- tifs[keep]
  years <- years[keep]

  cat("    Files:", length(tifs), "\n")

  annual_list <- list()
  for (y in sort(unique(years))) {
    year_tifs <- tifs[years == y]
    cat("      Year", y, ":", length(year_tifs), "files\n")
    r_year <- rast(year_tifs)

    # Remove fill values: real 8-day GPP is 0–0.1; fill values are ~3.2766
    r_year <- clamp(r_year, lower = 0, upper = 0.5, values = FALSE)
    annual_list[[as.character(y)]] <- sum(r_year, na.rm = TRUE)
  }

  r <- mean(rast(annual_list), na.rm = TRUE)

  cat("    Annual mean GPP range:",
      round(min(values(r), na.rm = TRUE), 3), "-",
      round(max(values(r), na.rm = TRUE), 3), "kg C/m2/yr\n")
  r
}

#' Compute validation metrics for paired vectors.
compute_metrics <- function(obs, mod) {
  valid <- is.finite(obs) & is.finite(mod) & obs > 0
  obs <- obs[valid]
  mod <- mod[valid]
  if (length(obs) < 10) return(NULL)

  residuals <- mod - obs
  list(
    n        = length(obs),
    obs_mean = mean(obs),
    mod_mean = mean(mod),
    bias     = mean(residuals),
    bias_pct = 100 * mean(residuals) / mean(obs),
    mae      = mean(abs(residuals)),
    rmse     = sqrt(mean(residuals^2)),
    r        = cor(obs, mod),
    r2       = cor(obs, mod)^2,
    mef      = 1 - sum(residuals^2) / sum((obs - mean(obs))^2)
  )
}

# =============================================================================
# LOAD MODEL (global raster, then crop per region)
# =============================================================================

cat("Loading LPJ-GUESS model...\n")
model_global <- rast(model_tif)

# =============================================================================
# PROCESS EACH REGION
# =============================================================================

summary_rows <- list()

for (reg_id in names(regions)) {

  reg <- regions[[reg_id]]
  cat("\n========================================\n")
  cat("Region:", reg$name, "-", reg$desc, "\n")
  cat("========================================\n")

  # Crop model
  model_crop <- crop(model_global, reg$ext)

  # Load MODIS
  modis_path <- file.path(modis_dir, reg_id)
  modis_raw  <- load_modis_region(modis_path, year_start, year_end)
  if (is.null(modis_raw)) {
    cat("  MODIS data not found at:", modis_path, "\n")
    cat("  Skipping. Download from AppEEARS for this region.\n")
    next
  }

  # Crop MODIS and resample to model grid
  modis_crop <- crop(modis_raw, reg$ext)
  if (res(modis_crop)[1] < res(model_crop)[1]) {
    cat("  Aggregating MODIS from", round(res(modis_crop)[1] * 111, 2),
        "km to", round(res(model_crop)[1] * 111, 2), "km...\n")
    modis_resamp <- resample(modis_crop, model_crop, method = "average")
  } else {
    modis_resamp <- resample(modis_crop, model_crop, method = "bilinear")
  }

  # Mask to common coverage
  common_mask  <- !is.na(model_crop) & !is.na(modis_resamp)
  model_masked <- mask(model_crop,   common_mask, maskvalue = FALSE)
  modis_masked <- mask(modis_resamp, common_mask, maskvalue = FALSE)

  # Extract paired values
  mod_vals <- as.numeric(values(model_masked, na.rm = TRUE))
  obs_vals <- as.numeric(values(modis_masked, na.rm = TRUE))
  cat("  Paired cells:", length(mod_vals), "\n")

  # Metrics
  m <- compute_metrics(obs_vals, mod_vals)
  if (is.null(m)) {
    cat("  Too few valid cells. Skipping.\n")
    next
  }
  summary_rows[[reg_id]] <- c(name = reg$name, m)

  cat(sprintf("  R² = %.3f | RMSE = %.3f | Bias = %.1f%% | MEF = %.3f\n",
              m$r2, m$rmse, m$bias_pct, m$mef))

  # ---- Side-by-side maps (model | MODIS | difference) ----
  all_vals <- c(values(model_masked, na.rm = TRUE), values(modis_masked, na.rm = TRUE))
  vmax <- quantile(all_vals, 0.99, na.rm = TRUE)

  p_model <- ggplot() +
    geom_spatraster(data = model_masked) +
    scale_fill_viridis_c(name = "GPP\n(kg C m⁻² yr⁻¹)",
                         limits = c(0, vmax), oob = scales::squish,
                         na.value = "transparent") +
    coord_sf() + labs(title = "LPJ-GUESS (1 km)") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"))

  p_modis <- ggplot() +
    geom_spatraster(data = modis_masked) +
    scale_fill_viridis_c(name = "GPP\n(kg C m⁻² yr⁻¹)",
                         limits = c(0, vmax), oob = scales::squish,
                         na.value = "transparent") +
    coord_sf() + labs(title = "MODIS MOD17 (500 m → 1 km)") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"))

  diff_rast <- model_masked - modis_masked
  dlim <- quantile(abs(values(diff_rast, na.rm = TRUE)), 0.95)

  p_diff <- ggplot() +
    geom_spatraster(data = diff_rast) +
    scale_fill_gradient2(name = "Bias\n(kg C m⁻² yr⁻¹)",
                         low = "brown", mid = "white", high = "darkgreen",
                         midpoint = 0, limits = c(-dlim, dlim),
                         oob = scales::squish, na.value = "transparent") +
    coord_sf() + labs(title = "Difference (Model − MODIS)") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"))

  p_maps <- p_model + p_modis + p_diff +
    plot_layout(ncol = 3) +
    plot_annotation(
      title    = paste0(reg$name, ": ", reg$desc),
      subtitle = sprintf("R² = %.3f | RMSE = %.3f | Bias = %.1f%% | n = %d cells",
                         m$r2, m$rmse, m$bias_pct, m$n)
    )

  ggsave(file.path(out_dir, paste0("regional_maps_", reg_id, ".png")),
         p_maps, width = reg$width * 3, height = reg$height, dpi = 300,
         limitsize = FALSE)

  # ---- Density scatter ----
  dt_scatter <- data.table(obs = obs_vals, mod = mod_vals)
  dt_scatter <- dt_scatter[is.finite(obs) & is.finite(mod) & obs > 0]
  if (nrow(dt_scatter) > 30000) dt_scatter <- dt_scatter[sample(.N, 30000)]

  p_scatter <- ggplot(dt_scatter, aes(x = obs, y = mod)) +
    geom_hex(bins = 80) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_fill_viridis_c(name = "Count", trans = "log10") +
    labs(
      title    = paste0(reg$name, ": LPJ-GUESS vs MODIS GPP (1 km)"),
      subtitle = sprintf("R² = %.3f | RMSE = %.3f | Bias = %.3f (%.1f%%)",
                         m$r2, m$rmse, m$bias, m$bias_pct),
      x = "MODIS GPP (kg C m⁻² yr⁻¹)",
      y = "LPJ-GUESS AGPP (kg C m⁻² yr⁻¹)"
    ) +
    coord_fixed(xlim = c(0, vmax), ylim = c(0, vmax)) +
    theme_minimal()

  ggsave(file.path(out_dir, paste0("regional_scatter_", reg_id, ".png")),
         p_scatter, width = 8, height = 8, dpi = 300)

  # ---- E-W transect at mid-latitude ----
  mid_lat <- mean(c(ext(reg$ext)[3], ext(reg$ext)[4]))
  coords_transect <- cbind(
    seq(ext(reg$ext)[1], ext(reg$ext)[2], length.out = 500),
    rep(mid_lat, 500)
  )

  transect_model <- extract(model_masked, coords_transect)
  transect_modis <- extract(modis_masked, coords_transect)

  dt_trans <- data.frame(
    lon       = coords_transect[, 1],
    LPJ_GUESS = as.numeric(transect_model[, ncol(transect_model)]),
    MODIS     = as.numeric(transect_modis[, ncol(transect_modis)])
  )
  dt_trans <- na.omit(dt_trans)
  dt_trans_long <- tidyr::pivot_longer(dt_trans, cols = c("LPJ_GUESS", "MODIS"),
                                       names_to = "source", values_to = "GPP")
  dt_trans_long <- dt_trans_long[!is.na(dt_trans_long$GPP), ]

  if (nrow(dt_trans_long) > 10) {
    p_transect <- ggplot(dt_trans_long, aes(x = lon, y = GPP, color = source)) +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      scale_color_manual(values = c("LPJ_GUESS" = "red", "MODIS" = "blue"), name = "") +
      labs(title = paste0(reg$name, ": E-W transect at ", round(mid_lat, 1), "°"),
           x = "Longitude", y = "GPP (kg C m⁻² yr⁻¹)") +
      theme_minimal()

    ggsave(file.path(out_dir, paste0("regional_transect_", reg_id, ".png")),
           p_transect, width = 10, height = 5, dpi = 300)
  }

  cat("  Plots saved.\n")
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n\n========== REGIONAL SUMMARY ==========\n")
cat(sprintf("%-20s %8s %8s %8s %8s %8s\n",
            "Region", "R²", "RMSE", "Bias%", "MAE", "MEF"))
cat("--------------------------------------------------------------\n")
for (row in summary_rows) {
  cat(sprintf("%-20s %8.3f %8.3f %7.1f%% %8.3f %8.3f\n",
              row$name, row$r2, row$rmse, row$bias_pct, row$mae, row$mef))
}
cat("==============================================\n")

cat("\nAll regional validation plots saved to:", out_dir, "\n")
