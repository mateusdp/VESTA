# =============================================================================
# plot_global.R — Read per-band simulation output, produce GeoTIFF + PNG maps
# =============================================================================
#
# Pipeline:
#   1. Read per-band .out.gz files, filter to target year range with awk
#   2. Compute temporal mean per grid cell
#   3. Rasterize to 30 arc-sec (1/120°) GeoTIFF (plain + DEFLATE-compressed)
#   4. Produce global map (PNG) with ggplot2/tidyterra
#   5. Optionally produce regional crop maps for publication figures
#
# Requires: data.table, terra, tidyterra, ggplot2
# =============================================================================

library(data.table)
library(terra)
library(tidyterra)
library(ggplot2)

# ---- CONFIG -----------------------------------------------------------------

base_dir      <- "/path/to/global_v2"       # Parent directory containing band folders
out_dir       <- "."                         # Output directory for TIF and PNG files
variable      <- "cmass"                     # Output file stem (e.g. "cmass", "agpp", "srl")
variable_name <- "Biomass (kg C m⁻² )"      # Label for color bar
col_name      <- "Total"                     # Column to extract (e.g. "Total", a PFT name)
year_start    <- 2017
year_end      <- 2021

bands <- c("60S_30S", "30S_20S", "20S_10S", "10S_00N", "00N_10N",
           "10N_20N", "20N_30N", "30N_40N", "40N_50N", "50N_60N", "60N_90N")

# Whether to produce regional crop maps (set FALSE for a quick run)
make_regional_maps <- TRUE

# =============================================================================
# 1. READ PER-BAND DATA
# =============================================================================

cat("Reading per-band .out.gz files...\n")

all_data <- rbindlist(lapply(bands, function(b) {
  gz_file <- file.path(base_dir, b, paste0(variable, ".out.gz"))

  if (!file.exists(gz_file)) {
    cat("  SKIP:", b, "(not found)\n")
    return(NULL)
  }

  cat("  Reading:", b, "...")
  # Filter years with awk before reading into R (avoids long-vector limit)
  dt <- fread(cmd = paste("zcat", gz_file,
                          "| awk 'NR==1 || ($3 >=", year_start,
                          "&& $3 <=", year_end, ")'"))
  cat("", nrow(dt), "rows\n")
  dt[, .(Lon, Lat, Year, value = get(col_name))]
}))

cat("\nTotal rows:", nrow(all_data), "\n")
cat("Unique cells:", uniqueN(all_data[, .(Lon, Lat)]), "\n")

# =============================================================================
# 2. COMPUTE TEMPORAL MEAN
# =============================================================================

global_mean <- all_data[, .(value = mean(value, na.rm = TRUE)), by = .(Lon, Lat)]
rm(all_data); gc()

q99 <- quantile(global_mean$value, 0.99, na.rm = TRUE)
cat("Value range:", range(global_mean$value, na.rm = TRUE), "\n")
cat("99th percentile:", q99, "\n")
cat("Global mean:", mean(global_mean$value, na.rm = TRUE), "\n")
cat("Tropical (10S-10N):", mean(global_mean[Lat > -10 & Lat < 10]$value, na.rm = TRUE), "\n")
cat("Boreal (50N-70N):", mean(global_mean[Lat > 50 & Lat < 70]$value, na.rm = TRUE), "\n")

# =============================================================================
# 3. RASTERIZE AND SAVE GEOTIFF
# =============================================================================

cat("\nRasterizing to 30 arc-sec...\n")

template <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                 res = 1/120, crs = "EPSG:4326")
pts <- vect(global_mean, geom = c("Lon", "Lat"), crs = "EPSG:4326")
r   <- rasterize(pts, template, field = "value", fun = mean)

tif_base <- file.path(out_dir, paste0(variable, "_global_", year_start, "_", year_end, "_mean"))

# Uncompressed (fast access)
writeRaster(r, paste0(tif_base, ".tif"), overwrite = TRUE)
cat("Saved:", paste0(tif_base, ".tif"), "\n")

# Compressed (archival / sharing)
writeRaster(r, paste0(tif_base, "_compress.tif"),
            datatype = "FLT4S", overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "ZLEVEL=9"))
cat("Saved:", paste0(tif_base, "_compress.tif"), "\n")

# =============================================================================
# 4. GLOBAL MAP
# =============================================================================

cat("\nProducing global map...\n")

p_global <- ggplot() +
  geom_spatraster(data = r, maxcell = 5e6) +
  scale_fill_viridis_c(
    name   = variable_name,
    limits = c(0, q99),
    oob    = scales::squish,
    na.value = "transparent"
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 90)) +
  borders(colour = "grey30", linewidth = 0.1) +
  theme_minimal()

global_png <- file.path(out_dir,
                        paste0(variable, "_", col_name, "_global_",
                               year_start, "-", year_end, ".png"))
ggsave(global_png, p_global, width = 20, height = 8, dpi = 300)
cat("Saved:", global_png, "\n")


# =============================================================================
# 5. REGIONAL CROP MAPS (optional)
# =============================================================================

if (!make_regional_maps) {
  cat("\nSkipping regional maps (make_regional_maps = FALSE).\n")
  cat("Done.\n")
  q("no")
}

# ---- Helper: crop-and-plot one region ----
plot_region <- function(raster, region, var_name, out_prefix,
                        color_limits = NULL, maxcell = Inf, dpi = 300) {
  cropped <- crop(raster, ext(region$xmin, region$xmax, region$ymin, region$ymax))

  fill_args <- list(name = var_name, na.value = "transparent")
  if (!is.null(color_limits)) {
    fill_args$limits <- color_limits
    fill_args$oob    <- scales::squish
  }

  p <- ggplot() +
    geom_spatraster(data = cropped, maxcell = maxcell) +
    do.call(scale_fill_viridis_c, fill_args) +
    coord_sf(expand = FALSE) +
    labs(title = region$title) +
    theme_minimal()

  if (isTRUE(region$borders)) {
    p <- p + borders(colour = "grey30", linewidth = 0.1)
  }

  outfile <- file.path(out_dir, paste0(out_prefix, "_", region$id, ".png"))
  ggsave(outfile, p, width = region$width, height = region$height, dpi = dpi)
  cat("  Saved:", outfile, "\n")
}

# ---- Region definitions ----
# Each entry: id, title, extent (xmin, xmax, ymin, ymax), plot dimensions

regions_showcase <- list(
  # --- Global biomes ---
  list(id = "andes",       title = "Andes: Amazon to Páramo",
       xmin = -80, xmax = -75, ymin = -2,  ymax = 4,  width = 8,  height = 10),
  list(id = "east_africa", title = "East Africa: Rift Valley & Kilimanjaro",
       xmin = 28,  xmax = 38,  ymin = -5,  ymax = 3,  width = 8,  height = 8),
  list(id = "se_asia",     title = "SE Asia: Borneo, Sumatra, Java",
       xmin = 95,  xmax = 120, ymin = -10, ymax = 8,  width = 12, height = 8),
  list(id = "caribbean",   title = "Caribbean Islands",
       xmin = -85, xmax = -60, ymin = 10,  ymax = 25, width = 12, height = 8),
  list(id = "sahel",       title = "Sahel: Sahara to Guinea Forest",
       xmin = -15, xmax = 15,  ymin = 5,   ymax = 20, width = 12, height = 8),
  list(id = "patagonia",   title = "Patagonia: Valdivian Forest & Steppe",
       xmin = -76, xmax = -65, ymin = -55, ymax = -38, width = 8, height = 12),
  list(id = "alps",        title = "Central Europe: Alps & Rhine Valley",
       xmin = 5,   xmax = 18,  ymin = 45,  ymax = 55, width = 10, height = 8),
  list(id = "japan",       title = "Japan: Coast to Mountain Gradients",
       xmin = 128, xmax = 146, ymin = 30,  ymax = 46, width = 8,  height = 10),
  list(id = "scandinavia", title = "Scandinavia: Boreal to Tundra Treeline",
       xmin = 15,  xmax = 30,  ymin = 65,  ymax = 72, width = 10, height = 8),
  list(id = "canada",      title = "Canada: Taiga to Tundra",
       xmin = -120, xmax = -90, ymin = 55, ymax = 68, width = 12, height = 8),

  # --- Europe (ESA showcase) ---
  list(id = "europe",       title = "Europe",
       xmin = -24, xmax = 35, ymin = 32, ymax = 71, width = 10, height = 8,
       borders = TRUE, maxcell = 5e6, dpi = 500),
  list(id = "iberia",       title = "Iberian Peninsula",
       xmin = -10, xmax = 4,  ymin = 36, ymax = 44, width = 10, height = 8),
  list(id = "italy",        title = "Italy: Alps to Mediterranean",
       xmin = 6,   xmax = 19, ymin = 36, ymax = 47, width = 8,  height = 10),
  list(id = "british_isles", title = "British Isles",
       xmin = -11, xmax = 2,  ymin = 50, ymax = 59, width = 8,  height = 10),
  list(id = "balkans",      title = "Balkans & Greece",
       xmin = 13,  xmax = 30, ymin = 35, ymax = 47, width = 10, height = 10),
  list(id = "fennoscandia", title = "Fennoscandia: Fjords to Boreal",
       xmin = 4,   xmax = 32, ymin = 55, ymax = 72, width = 10, height = 12),
  list(id = "rhine_danube", title = "Central Europe: Rhine to Danube",
       xmin = 5,   xmax = 17, ymin = 46, ymax = 55, width = 10, height = 8),
  list(id = "pyrenees_alps", title = "Pyrenees to Alps",
       xmin = -2,  xmax = 8,  ymin = 42, ymax = 47, width = 12, height = 6),
  list(id = "fjords",       title = "Norwegian Fjords",
       xmin = 4,   xmax = 9,  ymin = 59, ymax = 63, width = 8,  height = 8),
  list(id = "carpathians",  title = "Carpathian Mountains",
       xmin = 17,  xmax = 27, ymin = 44, ymax = 50, width = 10, height = 8),

  # --- Amazon ---
  list(id = "amazon",      title = "Amazon Basin",
       xmin = -80, xmax = -45, ymin = -20, ymax = 5, width = 12, height = 10,
       maxcell = 5e6),
  list(id = "amazon_arc",  title = "Amazon: Deforestation Arc",
       xmin = -65, xmax = -50, ymin = -15, ymax = -5, width = 10, height = 8),
  list(id = "manaus",      title = "Amazon: Manaus & River Confluence",
       xmin = -62, xmax = -58, ymin = -4,  ymax = -1, width = 8,  height = 8),

  # --- Oceania ---
  list(id = "australia",    title = "Australia",
       xmin = 112, xmax = 155, ymin = -45, ymax = -10, width = 12, height = 10),
  list(id = "se_australia", title = "SE Australia: Great Dividing Range",
       xmin = 143, xmax = 154, ymin = -39, ymax = -28, width = 8,  height = 10),
  list(id = "nz",           title = "New Zealand: Alps to Coast",
       xmin = 165, xmax = 179, ymin = -48, ymax = -34, width = 8,  height = 12),
  list(id = "tasmania",     title = "Tasmania",
       xmin = 143, xmax = 149, ymin = -44, ymax = -40, width = 8,  height = 8)
)

# ---- Generate all regional maps ----
cat("\nGenerating regional maps...\n")

for (reg in regions_showcase) {
  plot_region(
    raster       = r,
    region       = reg,
    var_name     = variable_name,
    out_prefix   = variable,
    color_limits = c(0, q99),
    maxcell      = if (!is.null(reg$maxcell)) reg$maxcell else Inf,
    dpi          = if (!is.null(reg$dpi)) reg$dpi else 300
  )
}

cat("\nDone. All maps saved to:", out_dir, "\n")
