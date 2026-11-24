# =============================================================================
# extract_trait_parameters.R
#
# Extracts trait ranges and tradeoff parameters from the TRY and GRooT
# databases for use in LPJ-GUESS-NTD.
#
# This script produces two CSV outputs consumed by the .ins file:
#   1. trait_ranges_pfts_IQR.csv   — SLA, WSG and SRL min/max per PFT
#      (after IQR-based outlier removal)
#   2. tradeoffs_All_PFT.csv      — fitted h0/h1 tradeoff parameters per PFT
#      (plus a global "All" fit) for eight trait–trait relationships
#
# The tradeoff parameters map to .ins parameters as follows:
#
#   Tradeoff (this script)             .ins parameters
#   ─────────────────────────          ──────────────────────
#   SLA ↔ Leaf C:N                     k_cton_a / k_cton_b
#   SLA ↔ Leaf C:P (from N:P)         k_ctop_a / k_ctop_b
#   SLA ↔ Leaf turnover               k_turnover_leaf_a / k_turnover_leaf_b
#   Root C:N ↔ Root turnover          k_turnover_root_a / k_turnover_root_b
#   SRL ↔ Root diameter               k_droot_a / k_droot_b
#   Root diameter ↔ Mycorrhizal col.  k_mycocol_a / k_mycocol_b
#   Leaf C:N ↔ Leaf turnover          (used in extended configuration)
#
# Two additional WSG tradeoffs present in the .ins file are NOT derived
# from this script and come from literature/calibration:
#   - k_turnover_sap_a / k_turnover_sap_b   (WSG ↔ sapwood turnover)
#   - k_wsg_mort_a / k_wsg_mort_b           (WSG ↔ mortality)
#
# Functional forms:
#   - Linear in log-log space:  log(Y) = h0 + h1 * log(X)
#     i.e. Y = exp(h0 + h1 * log(X))
#   - Sigmoid (mycorrhizal colonization only):
#     mcol = 1 / (1 + exp(h0 * (log(RootD) + h1)))
#
# PFT borrowing (insufficient data):
#   - Needleleaved Deciduous uses Needleleaved Evergreen for:
#       SLA↔C:P, Root C:N↔Root turnover, Root D↔Mycorrhizal colonization
#   - Broadleaved Evergreen uses Broadleaved Deciduous for:
#       Root C:N↔Root turnover
#
# Data sources:
#   - TRY database (https://www.try-db.org):
#       Required TraitIDs:
#         SLA        — 4083, 3115, 3116, 3117, 125
#         WSG        — 4, 3877, 3473, 3452, 3472, 3764, 3736, 1204, 3064, 3996
#         Leaf C:N   — 146
#         Leaf N:P   — 56  (C:P is derived as C:N × N:P)
#         Leaf longevity — 12
#       These can be obtained in a single TRY data request. If split
#       across multiple requests, update the file paths in Section 3.
#   - GRooT database (https://groot-database.github.io):
#       GRooTAggregateSpeciesVersion.csv: SRL, root diameter, root C:N,
#       root lifespan, mycorrhizal colonization
#
# Unit conversions applied:
#   - WSG: density (g/cm³) × 500 → kgC/m³
#   - SRL: m/g × 2000 → m/kgC
#   - Leaf longevity (months) → leaf turnover (yr⁻¹): 12 / LL
#   - Root lifespan (days) → root turnover (yr⁻¹): 1 / (lifespan / 365)
#   - Root diameter: mm / 1000 → m
#   - Mycorrhizal colonization: percentage / 100 → fraction
#
# Requirements:
#   R packages: data.table, rtry, tidyr, dplyr, ggplot2, ggridges, cowplot
#
# Usage:
#   Update the file paths in Section 1 below, then source this script.
#
# Part of: lpj-guess-ntd preprocessing pipeline
# =============================================================================


# =============================================================================
# 1. Setup
# =============================================================================

library(data.table)
library(rtry)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(cowplot)

# --- Configurable paths (UPDATE THESE) ---
try_data_dir  <- "~/VESTA/trait_data/"           # Directory containing TRY data files
try_file_main <- "37325.txt"                     # TRY file with SLA, WSG, leaf longevity traits
try_file_cn_np <- "44739.txt"                    # TRY file with C:N (146) and N:P (56) traits
                                                 # (can be the same file as try_file_main)
groot_csv     <- "~/data/GRooT/GRooTAggregateSpeciesVersion.csv"  # GRooT aggregate file
output_dir    <- "~/VESTA/trait_data/"            # Where to write output CSVs
plot_dir      <- "~/VESTA/trait_data/"            # Where to write diagnostic plots

# --- Plot theme ---
std_size <- 20
textsizes <- theme(
  axis.title.x  = element_text(size = std_size),
  axis.text.x   = element_text(size = std_size),
  axis.title.y  = element_text(size = std_size),
  axis.text.y   = element_text(size = std_size),
  legend.text   = element_text(size = std_size),
  legend.title  = element_text(size = std_size),
  strip.text    = element_text(size = std_size)
)

# --- PFT names (consistent ordering) ---
pft_names <- c("Broadleaved Deciduous", "Broadleaved Evergreen",
               "C3 Grasses", "C4 Grasses",
               "Needleleaved Deciduous", "Needleleaved Evergreen")


# =============================================================================
# 2. PFT classification from TRY categorical traits
# =============================================================================

cat("Loading TRY categorical traits lookup table...\n")
cat_traits_lookup <- fread(
  file.path(try_data_dir, "TRY_Categorical_Traits_Lookup_Table.csv"),
  sep = ",", header = TRUE, fill = TRUE
)

# Classify species into PFTs based on growth form, leaf type, phenology, and
# photosynthetic pathway
cat_traits_lookup[, PFT := fcase(
  # Woody broadleaved evergreen (default phenology for broadleaved woody)
  (PlantGrowthForm %in% c("tree", "shrub/tree", "shrub") | Woodiness == "woody") &
    LeafType == "broadleaved" &
    (LeafPhenology == "evergreen" | LeafPhenology == ""),
  "Broadleaved Evergreen",

  # Woody needleleaved evergreen
  (PlantGrowthForm %in% c("tree", "shrub/tree", "shrub") | Woodiness == "woody") &
    LeafType %in% c("needleleaved", "scale-shaped", "scale-shaped/needleleaved") &
    LeafPhenology == "evergreen",
  "Needleleaved Evergreen",

  # Woody broadleaved deciduous
  (PlantGrowthForm %in% c("tree", "shrub/tree", "shrub") | Woodiness == "woody") &
    LeafType == "broadleaved" &
    LeafPhenology %in% c("deciduous", "deciduous/evergreen"),
  "Broadleaved Deciduous",

  # Woody needleleaved deciduous
  (PlantGrowthForm %in% c("tree", "shrub/tree", "shrub") | Woodiness == "woody") &
    LeafType %in% c("needleleaved", "scale-shaped", "scale-shaped/needleleaved") &
    LeafPhenology %in% c("deciduous", "deciduous/evergreen"),
  "Needleleaved Deciduous",

  # C3 grasses / herbs (default pathway for non-woody)
  (PlantGrowthForm %in% c("graminoid", "herb/shrub", "herb", "") | Woodiness == "non-woody") &
    (PhotosyntheticPathway == "C3" | PhotosyntheticPathway == ""),
  "C3 Grasses",

  # C4 grasses / herbs
  (PlantGrowthForm %in% c("graminoid", "herb/shrub", "herb", "") | Woodiness == "non-woody") &
    PhotosyntheticPathway == "C4",
  "C4 Grasses"
)]

# Lookup tables for joining with TRY data (by AccSpeciesID) and GRooT (by name)
cat_traits_lookup_pft <- cat_traits_lookup[!is.na(PFT), .(AccSpeciesID, PFT)]
cat_traits_lookup_pft_names <- cat_traits_lookup[!is.na(PFT), .(genusspecies = AccSpeciesName, PFT)]

# Save lookup tables
fwrite(cat_traits_lookup_pft, file.path(output_dir, "cat_traits_lookup_pft.csv"))
fwrite(cat_traits_lookup_pft_names, file.path(output_dir, "cat_traits_lookup_pft_sspnames.csv"))

cat("  PFT classification done:", nrow(cat_traits_lookup_pft), "species classified.\n")


# =============================================================================
# 3. Load trait databases
# =============================================================================

# TRY trait data — one or two files depending on how the data was requested.
# If all TraitIDs are in a single file, point both variables to the same file.
cat("Loading TRY data (SLA, WSG, leaf longevity)...\n")
TRYdataAll <- rtry_import(file.path(try_data_dir, try_file_main))

cat("Loading TRY data (C:N and N:P for derived C:P)...\n")
TRYdataCN_NP <- rtry_import(file.path(try_data_dir, try_file_cn_np))

cat("Loading GRooT aggregate species data...\n")
groot_data <- fread(groot_csv, header = TRUE, dec = ".")
groot_data[, genusspecies := paste(genusTNRS, speciesTNRS)]


# =============================================================================
# 4. Helper functions
# =============================================================================

#' Extract a TRY trait and compute species-level averages
#'
#' @param try_data  Full TRY dataset (from rtry_import)
#' @param trait_ids Numeric vector of TraitIDs to select
#' @param col_name  Name for the averaged trait column
#' @param scale     Multiplicative scaling factor (default 1)
#' @return data.table with AccSpeciesID and averaged trait value
extract_try_trait <- function(try_data, trait_ids, col_name, scale = 1) {
  trait_data <- rtry_select_row(try_data, TraitID %in% trait_ids)
  trait_explore <- rtry_explore(trait_data,
                                AccSpeciesID, TraitName, UnitName, StdValue,
                                sortBy = AccSpeciesID)
  trait_dt <- as.data.table(trait_explore)
  trait_dt[, .(value = mean(StdValue) * scale), by = "AccSpeciesID"]  |>
    setnames("value", col_name)
}

#' Extract a GRooT trait and compute species-level averages
#'
#' @param groot     GRooT data.table
#' @param trait     Trait name in GRooT (traitName column)
#' @param col_name  Name for the averaged trait column
#' @param scale     Multiplicative scaling factor (default 1)
#' @return data.table with genusspecies and averaged trait value
extract_groot_trait <- function(groot, trait, col_name, scale = 1) {
  trait_dt <- groot[traitName == trait]
  trait_dt[, .(value = mean(as.numeric(meanSpecies)) * scale), by = "genusspecies"] |>
    setnames("value", col_name)
}

#' IQR-based outlier removal and summary (min/max/mean) per PFT
#'
#' @param dt       data.table with PFT column and trait column
#' @param trait_col Name of the trait column
#' @return data.table with min, max, mean per PFT (after IQR filtering)
summarise_trait_ranges <- function(dt, trait_col) {
  dt[, {
    vals <- get(trait_col)
    Q1 <- quantile(vals, 0.25, na.rm = TRUE)
    Q3 <- quantile(vals, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower <- Q1 - 1.5 * IQR_val
    upper <- Q3 + 1.5 * IQR_val
    filtered <- vals[vals >= lower & vals <= upper]
    list(
      min = min(filtered, na.rm = TRUE),
      max = max(filtered, na.rm = TRUE),
      mean = mean(filtered, na.rm = TRUE)
    )
  }, by = PFT] |>
    setnames(c("PFT",
               paste0("min_", trait_col),
               paste0("max_", trait_col),
               paste0("mean_", trait_col)))
}

#' Fit linear tradeoff in log-log space per PFT
#'
#' Fits log(Y) = h0 + h1 * log(X) using NLS per PFT plus a global fit.
#' Both X and Y should already be log-transformed in the input data.
#'
#' @param dt           data.table with PFT, x_col, and y_col (all log-transformed)
#' @param x_col        Name of predictor column
#' @param y_col        Name of response column
#' @param tradeoff_id  Label for the tradeoff (e.g. "LeafCN-SLA")
#' @param equation     Equation string for documentation
#' @param borrowing    Named list of PFT borrowing (e.g. list(ND = "NE"))
#' @param starts       Named list of start values per PFT (optional overrides)
#' @return data.table with PFT, h0, h1, R2, tradeoff, equation
fit_linear_tradeoff <- function(dt, x_col, y_col, tradeoff_id, equation,
                                 borrowing = list(),
                                 starts = list()) {
  hdFunc <- function(x, h0, h1) { h0 + (x * h1) }

  # Short PFT abbreviations for borrowing
  pft_abbrev <- c("BD" = "Broadleaved Deciduous",
                   "BE" = "Broadleaved Evergreen",
                   "C3G" = "C3 Grasses",
                   "C4G" = "C4 Grasses",
                   "ND" = "Needleleaved Deciduous",
                   "NE" = "Needleleaved Evergreen")

  # Default starting values
  default_starts <- list(h0 = 10, h1 = 0)

  # Fit per PFT
  fits <- list()
  for (pft in pft_names) {
    pft_data <- dt[PFT == pft]
    start_vals <- if (pft %in% names(starts)) starts[[pft]] else default_starts
    fits[[pft]] <- tryCatch(
      nls(as.formula(paste(y_col, "~ hdFunc(", x_col, ", h0, h1)")),
          data = pft_data, start = start_vals),
      error = function(e) {
        warning(paste("NLS failed for", pft, "in", tradeoff_id, ":", e$message))
        NULL
      }
    )
  }

  # Apply PFT borrowing
  for (from_abbrev in names(borrowing)) {
    to_abbrev <- borrowing[[from_abbrev]]
    from_pft <- pft_abbrev[from_abbrev]
    to_pft <- pft_abbrev[to_abbrev]
    cat("  Borrowing:", from_pft, "<-", to_pft, "\n")
    fits[[from_pft]] <- fits[[to_pft]]
  }

  # Extract parameters and R²
  results <- rbindlist(lapply(pft_names, function(pft) {
    fit <- fits[[pft]]
    if (is.null(fit)) return(NULL)
    h0 <- summary(fit)$parameters[1]
    h1 <- summary(fit)$parameters[2]
    pft_data <- dt[PFT == pft]
    r2 <- summary(lm(
      as.formula(paste(y_col, "~ hdFunc(", x_col, ",", h0, ",", h1, ")")),
      data = pft_data
    ))$r.squared
    data.table(PFT = pft, h0 = h0, h1 = h1, R2 = r2)
  }))

  # Global fit ("All")
  fit_all <- nls(as.formula(paste(y_col, "~ hdFunc(", x_col, ", h0, h1)")),
                 data = dt, start = default_starts)
  h0_all <- summary(fit_all)$parameters[1]
  h1_all <- summary(fit_all)$parameters[2]
  r2_all <- summary(lm(
    as.formula(paste(y_col, "~ hdFunc(", x_col, ",", h0_all, ",", h1_all, ")")),
    data = dt
  ))$r.squared
  results <- rbind(results, data.table(PFT = "All", h0 = h0_all, h1 = h1_all, R2 = r2_all))

  results$tradeoff <- tradeoff_id
  results$equation <- equation
  return(results)
}

#' Fit sigmoid tradeoff for mycorrhizal colonization
#'
#' Fits mcol = 1 / (1 + exp(h0 * (log(x) + h1))) per PFT plus a global fit.
#' X (root diameter) is NOT log-transformed; Y (mycorrhizal colonisation) is a
#' fraction (0–1).
#'
#' @param dt           data.table with PFT, x_col, y_col
#' @param x_col        Root diameter column (metres, NOT log-transformed)
#' @param y_col        Mycorrhizal colonisation column (fraction)
#' @param borrowing    Named list of PFT borrowing
#' @return data.table with PFT, h0, h1, R2, tradeoff, equation
fit_sigmoid_tradeoff <- function(dt, x_col, y_col, borrowing = list()) {
  hdFunc_sig <- function(x, h0, h1) { 1.0 / (1.0 + exp(h0 * (log(x) + h1))) }

  pft_abbrev <- c("BD" = "Broadleaved Deciduous",
                   "BE" = "Broadleaved Evergreen",
                   "C3G" = "C3 Grasses",
                   "C4G" = "C4 Grasses",
                   "ND" = "Needleleaved Deciduous",
                   "NE" = "Needleleaved Evergreen")

  fits <- list()
  for (pft in pft_names) {
    pft_data <- dt[PFT == pft]
    fits[[pft]] <- tryCatch(
      nls(as.formula(paste(y_col, "~ hdFunc_sig(", x_col, ", h0, h1)")),
          data = pft_data, start = list(h0 = -3, h1 = 7.68)),
      error = function(e) {
        warning(paste("NLS sigmoid failed for", pft, ":", e$message))
        NULL
      }
    )
  }

  # Apply PFT borrowing
  for (from_abbrev in names(borrowing)) {
    to_abbrev <- borrowing[[from_abbrev]]
    from_pft <- pft_abbrev[from_abbrev]
    to_pft <- pft_abbrev[to_abbrev]
    cat("  Borrowing:", from_pft, "<-", to_pft, "\n")
    fits[[from_pft]] <- fits[[to_pft]]
  }

  results <- rbindlist(lapply(pft_names, function(pft) {
    fit <- fits[[pft]]
    if (is.null(fit)) return(NULL)
    h0 <- summary(fit)$parameters[1]
    h1 <- summary(fit)$parameters[2]
    pft_data <- dt[PFT == pft]
    r2 <- summary(lm(
      as.formula(paste(y_col, "~ hdFunc_sig(", x_col, ",", h0, ",", h1, ")")),
      data = pft_data
    ))$r.squared
    data.table(PFT = pft, h0 = h0, h1 = h1, R2 = r2)
  }))

  # Global fit
  fit_all <- nls(as.formula(paste(y_col, "~ hdFunc_sig(", x_col, ", h0, h1)")),
                 data = dt, start = list(h0 = -2, h1 = -1))
  h0_all <- summary(fit_all)$parameters[1]
  h1_all <- summary(fit_all)$parameters[2]
  r2_all <- summary(lm(
    as.formula(paste(y_col, "~ hdFunc_sig(", x_col, ",", h0_all, ",", h1_all, ")")),
    data = dt
  ))$r.squared
  results <- rbind(results, data.table(PFT = "All", h0 = h0_all, h1 = h1_all, R2 = r2_all))

  results$tradeoff <- "RootD-Mcol"
  results$equation <- "Mcol = 1.0 / (1.0 + exp(h0 * (log(RootD) + h1)))"
  return(results)
}

#' Generate per-PFT and global tradeoff plots
#'
#' @param dt           data.table with x, y, PFT, h0, h1, R2 columns
#' @param x_col        X axis column name (log-transformed)
#' @param y_col        Y axis column name (log-transformed)
#' @param x_lab        X axis label
#' @param y_lab        Y axis label
#' @param func         Prediction function (default: linear)
#' @return list with per_pft and global ggplot objects
plot_tradeoff <- function(dt, x_col, y_col, x_lab, y_lab,
                           func = function(x, h0, h1) h0 + x * h1) {
  # Per-PFT plot
  curve_data <- dt %>%
    group_by(PFT) %>%
    summarise(
      x_min = min(.data[[x_col]], na.rm = TRUE),
      x_max = max(.data[[x_col]], na.rm = TRUE),
      h0 = first(h0), h1 = first(h1),
      .groups = "drop"
    ) %>%
    filter(is.finite(x_min), is.finite(x_max)) %>%
    rowwise() %>%
    mutate(x_seq = list(seq(x_min, x_max, length.out = 200))) %>%
    unnest(cols = c(x_seq)) %>%
    mutate(y_pred = func(x_seq, h0, h1))

  annot <- dt %>%
    group_by(PFT, R2) %>%
    summarise(x = Inf, y = Inf, .groups = "drop") %>%
    mutate(label = paste("R² =", round(R2, 3)))

  p_pft <- ggplot(dt, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point() +
    geom_line(data = curve_data, aes(x = x_seq, y = y_pred), color = "blue", linewidth = 2) +
    textsizes +
    xlab(x_lab) + ylab(y_lab) +
    geom_text(data = annot, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1.1, vjust = 2, size = 5) +
    facet_wrap(~PFT)

  # Global plot (need global h0, h1, R2 — from "All" row passed separately)
  p_global <- NULL  # filled by caller if needed


  list(per_pft = p_pft, global = p_global)
}


# =============================================================================
# 5. Trait range extraction (SLA, WSG, SRL)
# =============================================================================

cat("\n--- Extracting trait ranges ---\n")

# --- SLA ---
cat("  SLA (TRY TraitIDs: 4083, 3115, 3116, 3117, 125)...\n")
SLA_avg <- extract_try_trait(TRYdataAll, c(4083, 3115, 3116, 3117, 125), "SLA_ssp_avg")
SLA_pft <- merge(SLA_avg, cat_traits_lookup_pft, by = "AccSpeciesID")
SLA_ranges <- summarise_trait_ranges(SLA_pft, "SLA_ssp_avg")
cat("    Species:", nrow(SLA_avg), " | With PFT:", nrow(SLA_pft), "\n")

# --- WSG ---
# WSG density (g/cm³) × 500 = kgC/m³
cat("  WSG (TRY TraitIDs: 4, 3877, 3473, 3452, 3472, 3764, 3736, 1204, 3064, 3996)...\n")
WSG_avg <- extract_try_trait(TRYdataAll,
                              c(4, 3877, 3473, 3452, 3472, 3764, 3736, 1204, 3064, 3996),
                              "WSG_ssp_avg", scale = 500)
WSG_pft <- merge(WSG_avg, cat_traits_lookup_pft, by = "AccSpeciesID")
WSG_ranges <- summarise_trait_ranges(WSG_pft, "WSG_ssp_avg")
cat("    Species:", nrow(WSG_avg), " | With PFT:", nrow(WSG_pft), "\n")

# --- SRL ---
# SRL from GRooT: m/g × 2000 = m/kgC
cat("  SRL (GRooT: Specific_root_length)...\n")
SRL_avg <- extract_groot_trait(groot_data, "Specific_root_length", "SRL", scale = 2000)
SRL_pft <- merge(SRL_avg, cat_traits_lookup_pft_names, by = "genusspecies")
SRL_ranges <- summarise_trait_ranges(SRL_pft, "SRL")
cat("    Species:", nrow(SRL_avg), " | With PFT:", nrow(SRL_pft), "\n")

# --- Combine and save ---
trait_ranges <- merge(SLA_ranges, WSG_ranges, by = "PFT")
trait_ranges <- merge(trait_ranges, SRL_ranges, by = "PFT")
fwrite(trait_ranges, file.path(output_dir, "trait_ranges_pfts_IQR.csv"))
cat("  Saved: trait_ranges_pfts_IQR.csv\n")


# =============================================================================
# 6. Tradeoff parameter fitting
# =============================================================================

cat("\n--- Fitting tradeoff relationships ---\n")

# ----- 6a. SLA ↔ Leaf C:N (TRY TraitID 146) -----
# .ins: k_cton_a, k_cton_b

cat("\n  [1/7] SLA <-> Leaf C:N\n")

CN_avg <- extract_try_trait(TRYdataAll, 146, "CN_ssp_avg")
SLA_CN <- merge(SLA_avg, CN_avg, by = "AccSpeciesID")
SLA_CN <- merge(SLA_CN, cat_traits_lookup_pft, by = "AccSpeciesID")

# Log-transform
SLA_CN[, SLA_ssp_avg := log(SLA_ssp_avg)]
SLA_CN[, CN_ssp_avg := log(CN_ssp_avg)]

tradeoffs_SLA_CN <- fit_linear_tradeoff(
  SLA_CN, "SLA_ssp_avg", "CN_ssp_avg",
  tradeoff_id = "LeafCN-SLA",
  equation = "LeafCN = exp(h0 + log(SLA) * h1)"
)
cat("    Rows:", nrow(SLA_CN), "\n")
print(tradeoffs_SLA_CN)


# ----- 6b. SLA ↔ Leaf C:P (derived from C:N × N:P) -----
# .ins: k_ctop_a, k_ctop_b
# C:P = C:N × N:P (TraitID 146 × TraitID 56)
# Note: This overwrites any direct C:P fit (TraitID 151) because the
# N:P-derived version uses a larger sample and is used in the .ins file.

cat("\n  [2/7] SLA <-> Leaf C:P (from N:P)\n")

# C:N (TraitID 146)
CN_NP_avg <- extract_try_trait(TRYdataCN_NP, 146, "CN_ssp_avg")

# N:P (TraitID 56)
NP_data <- rtry_select_row(TRYdataCN_NP, TraitID %in% 56)
NP_explore <- as.data.table(rtry_explore(NP_data,
                                          AccSpeciesID, TraitName, UnitName, StdValue,
                                          sortBy = AccSpeciesID))
NP_avg <- NP_explore[, .(NP_ssp_avg = mean(StdValue)), by = "AccSpeciesID"]

# Derive C:P = C:N × N:P
CN_NP <- merge(CN_NP_avg, NP_avg, by = "AccSpeciesID")
CN_NP[, CP_ssp_avg := CN_ssp_avg * NP_ssp_avg]

# Merge with SLA and PFT
SLA_CP <- merge(SLA_avg, CN_NP[, .(AccSpeciesID, CP_ssp_avg)], by = "AccSpeciesID")
SLA_CP <- merge(SLA_CP, cat_traits_lookup_pft, by = "AccSpeciesID")

# Log-transform
SLA_CP[, SLA_ssp_avg := log(SLA_ssp_avg)]
SLA_CP[, CP_ssp_avg := log(CP_ssp_avg)]

tradeoffs_SLA_CP <- fit_linear_tradeoff(
  SLA_CP, "SLA_ssp_avg", "CP_ssp_avg",
  tradeoff_id = "LeafCP-SLA",
  equation = "LeafCP = exp(h0 + log(SLA) * h1)",
  borrowing = list("ND" = "NE")  # ND uses NE fit (insufficient data)
)
cat("    Rows:", nrow(SLA_CP), "\n")
print(tradeoffs_SLA_CP)


# ----- 6c. SLA ↔ Leaf turnover (TRY TraitID 12: leaf longevity) -----
# .ins: k_turnover_leaf_a, k_turnover_leaf_b
# Leaf longevity (months) → turnover (yr⁻¹) = 12 / LL

cat("\n  [3/7] SLA <-> Leaf turnover\n")

LL_data <- rtry_select_row(TRYdataAll, TraitID %in% 12)
LL_explore <- as.data.table(rtry_explore(LL_data,
                                          AccSpeciesID, TraitName, UnitName, StdValue,
                                          sortBy = AccSpeciesID))
# Convert leaf longevity (months) to leaf turnover (yr⁻¹)
LT_avg <- LL_explore[, .(LT_ssp_avg = 12 / mean(StdValue)), by = "AccSpeciesID"]

SLA_LT <- merge(SLA_avg, LT_avg, by = "AccSpeciesID")
SLA_LT <- merge(SLA_LT, cat_traits_lookup_pft, by = "AccSpeciesID")
SLA_LT <- SLA_LT[PFT != ""]

# Log-transform
SLA_LT[, SLA_ssp_avg := log(SLA_ssp_avg)]
SLA_LT[, LT_ssp_avg := log(LT_ssp_avg)]

tradeoffs_SLA_LT <- fit_linear_tradeoff(
  SLA_LT, "SLA_ssp_avg", "LT_ssp_avg",
  tradeoff_id = "LT-SLA",
  equation = "LT = exp(h0 + log(SLA) * h1)"
)
cat("    Rows:", nrow(SLA_LT), "\n")
print(tradeoffs_SLA_LT)


# ----- 6d. SRL ↔ Root diameter (GRooT) -----
# .ins: k_droot_a, k_droot_b
# Root diameter: mm / 1000 → m

cat("\n  [4/7] SRL <-> Root diameter\n")

RootD_avg <- extract_groot_trait(groot_data, "Mean_Root_diameter", "Root_D", scale = 1 / 1000)
RootD_pft <- merge(RootD_avg, cat_traits_lookup_pft_names, by = "genusspecies")

SRL_RootD <- merge(SRL_pft, RootD_pft, by = c("genusspecies", "PFT"))

# Log-transform
SRL_RootD[, SRL := log(SRL)]
SRL_RootD[, Root_D := log(Root_D)]

tradeoffs_SRL_RootD <- fit_linear_tradeoff(
  SRL_RootD, "SRL", "Root_D",
  tradeoff_id = "RootD-SRL",
  equation = "RootD = exp(h0 + log(SRL) * h1)"
)
cat("    Rows:", nrow(SRL_RootD), "\n")
print(tradeoffs_SRL_RootD)


# ----- 6e. Root C:N ↔ Root turnover (GRooT) -----
# .ins: k_turnover_root_a, k_turnover_root_b
# Root lifespan (days) → root turnover (yr⁻¹) = 1 / (lifespan / 365)

cat("\n  [5/7] Root C:N <-> Root turnover\n")

RCN_avg <- extract_groot_trait(groot_data, "Root_C_N_ratio", "Root_CN")
RCN_pft <- merge(RCN_avg, cat_traits_lookup_pft_names, by = "genusspecies")

# Root lifespan → root turnover
root_lifespan <- groot_data[traitName == "Root_lifespan_mean"]
RT_avg <- root_lifespan[, .(turnroot = 1 / (mean(as.numeric(meanSpecies)) / 365)),
                         by = "genusspecies"]
RT_pft <- merge(RT_avg, cat_traits_lookup_pft_names, by = "genusspecies")

RCN_RT <- merge(RCN_pft, RT_pft, by = c("genusspecies", "PFT"))

# Log-transform
RCN_RT[, Root_CN := log(Root_CN)]
RCN_RT[, turnroot := log(turnroot)]

tradeoffs_RCN_RT <- fit_linear_tradeoff(
  RCN_RT, "Root_CN", "turnroot",
  tradeoff_id = "RT-RootCN",
  equation = "RT = exp(h0 + log(rootCN) * h1)",
  borrowing = list("ND" = "NE", "BE" = "BD")  # ND uses NE, BE uses BD
)
cat("    Rows:", nrow(RCN_RT), "\n")
print(tradeoffs_RCN_RT)


# ----- 6f. Root diameter ↔ Mycorrhizal colonisation (GRooT) -----
# .ins: k_mycocol_a, k_mycocol_b
# Sigmoid: mcol = 1 / (1 + exp(h0 * (log(RootD) + h1)))
# Root diameter in metres (NOT log-transformed for this fit)
# Colonisation percentage / 100 → fraction

cat("\n  [6/7] Root diameter <-> Mycorrhizal colonisation (sigmoid)\n")

MCOL_avg <- extract_groot_trait(groot_data, "Root_mycorrhizal colonization", "mcol",
                                 scale = 1 / 100)
MCOL_pft <- merge(MCOL_avg, cat_traits_lookup_pft_names, by = "genusspecies")

# Root diameter (already extracted above, re-use — but need non-log version)
RootD_avg_m <- extract_groot_trait(groot_data, "Mean_Root_diameter", "Root_D", scale = 1 / 1000)
RootD_pft_m <- merge(RootD_avg_m, cat_traits_lookup_pft_names, by = "genusspecies")

MRD_MCOL <- merge(MCOL_pft, RootD_pft_m, by = c("genusspecies", "PFT"))
MRD_MCOL <- MRD_MCOL[PFT != ""]

tradeoffs_MRD_MCOL <- fit_sigmoid_tradeoff(
  MRD_MCOL, "Root_D", "mcol",
  borrowing = list("ND" = "NE")
)
cat("    Rows:", nrow(MRD_MCOL), "\n")
print(tradeoffs_MRD_MCOL)


# ----- 6g. Leaf C:N ↔ Leaf turnover (TRY) -----
# Used in extended model configuration

cat("\n  [7/7] Leaf C:N <-> Leaf turnover\n")

LCN_LT <- merge(CN_avg, LT_avg, by = "AccSpeciesID")
LCN_LT <- merge(LCN_LT, cat_traits_lookup_pft, by = "AccSpeciesID")
LCN_LT <- LCN_LT[PFT != ""]

# Log-transform
LCN_LT[, CN_ssp_avg := log(CN_ssp_avg)]
LCN_LT[, LT_ssp_avg := log(LT_ssp_avg)]

tradeoffs_LCN_LT <- fit_linear_tradeoff(
  LCN_LT, "CN_ssp_avg", "LT_ssp_avg",
  tradeoff_id = "LT-LeafCN",
  equation = "LT = exp(h0 + log(leafCN) * h1)"
)
cat("    Rows:", nrow(LCN_LT), "\n")
print(tradeoffs_LCN_LT)


# =============================================================================
# 7. Combine and save all tradeoff parameters
# =============================================================================

cat("\n--- Combining tradeoff parameters ---\n")

tradeoffs_all <- rbind(
  tradeoffs_SLA_CN,
  tradeoffs_SLA_CP,
  tradeoffs_SLA_LT,
  tradeoffs_RCN_RT,
  tradeoffs_SRL_RootD,
  tradeoffs_MRD_MCOL,
  tradeoffs_LCN_LT
)

fwrite(tradeoffs_all, file.path(output_dir, "tradeoffs_All_PFT.csv"))
cat("  Saved: tradeoffs_All_PFT.csv\n")

cat("\n  Summary of fitted tradeoffs:\n")
print(tradeoffs_all[PFT == "All", .(tradeoff, h0, h1, R2)])


# =============================================================================
# 8. Diagnostic plots (optional)
# =============================================================================

cat("\n--- Generating diagnostic plots ---\n")

# Ridge density plots for trait distributions
theme_set(theme_light())

sla_ridge <- ggplot(SLA_pft[PFT %in% pft_names], aes(x = SLA_ssp_avg, y = PFT, fill = PFT)) +
  geom_density_ridges(alpha = 0.7, scale = 3, color = "black") +
  textsizes + xlab("SLA (m²/kgC)") + ylab("") +
  theme(legend.position = "none")

wsg_ridge <- ggplot(WSG_pft[PFT %in% pft_names], aes(x = WSG_ssp_avg, y = PFT, fill = PFT)) +
  geom_density_ridges(alpha = 0.7, scale = 3, color = "black") +
  textsizes + xlab("WSG (kgC/m³)") + ylab("") +
  theme(legend.position = "none")

srl_ridge <- ggplot(SRL_pft[PFT %in% pft_names], aes(x = SRL, y = PFT, fill = PFT)) +
  geom_density_ridges(alpha = 0.7, scale = 3, color = "black") +
  textsizes + xlab("SRL (m/kgC)") + ylab("") +
  theme(legend.position = "none")

# Tradeoff scatter plots with fitted curves
# (Using the hdFunc definitions from the fitting functions)
hdFunc_linear <- function(x, h0, h1) h0 + x * h1
hdFunc_sigmoid <- function(x, h0, h1) 1 / (1 + exp(h0 * (log(x) + h1)))

# Helper to make per-PFT + global plots for a linear tradeoff
make_tradeoff_plots <- function(dt, tradeoffs, x_col, y_col, x_lab, y_lab,
                                 func = hdFunc_linear) {
  # Merge fitted params into data
  dt_plot <- merge(dt, tradeoffs[PFT != "All", .(PFT, h0, h1, R2)], by = "PFT")
  all_params <- tradeoffs[PFT == "All"]

  # Per-PFT
  curve_data <- dt_plot %>%
    group_by(PFT) %>%
    summarise(x_min = min(.data[[x_col]], na.rm = TRUE),
              x_max = max(.data[[x_col]], na.rm = TRUE),
              h0 = first(h0), h1 = first(h1), .groups = "drop") %>%
    filter(is.finite(x_min), is.finite(x_max)) %>%
    rowwise() %>%
    mutate(x_seq = list(seq(x_min, x_max, length.out = 200))) %>%
    unnest(cols = c(x_seq)) %>%
    mutate(y_pred = func(x_seq, h0, h1))

  annot <- dt_plot %>%
    group_by(PFT, R2) %>%
    summarise(x = Inf, y = Inf, .groups = "drop") %>%
    mutate(label = paste("R² =", round(R2, 3)))

  p_pft <- ggplot(dt_plot, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point() +
    geom_line(data = curve_data, aes(x = x_seq, y = y_pred), color = "blue", linewidth = 2) +
    textsizes + xlab(x_lab) + ylab(y_lab) +
    geom_text(data = annot, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1.1, vjust = 2, size = 5) +
    facet_wrap(~PFT)

  # Global
  curve_global <- data.frame(
    x_seq = seq(min(dt[[x_col]], na.rm = TRUE), max(dt[[x_col]], na.rm = TRUE), length.out = 200)
  ) %>%
    mutate(y_pred = func(x_seq, all_params$h0, all_params$h1))

  p_global <- ggplot(dt, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point() +
    geom_line(data = curve_global, aes(x = x_seq, y = y_pred), color = "red", linewidth = 1) +
    textsizes +
    annotate("text", x = Inf, y = Inf,
             label = paste("R² =", round(all_params$R2, 4)),
             hjust = 1.1, vjust = 2, size = 6) +
    xlab(x_lab) + ylab(y_lab)

  list(per_pft = p_pft, global = p_global)
}

# Generate plots for each tradeoff
p_sla_cn <- make_tradeoff_plots(SLA_CN, tradeoffs_SLA_CN,
                                 "SLA_ssp_avg", "CN_ssp_avg",
                                 "log SLA (m²/kgC)", "log C:N ratio")

p_sla_cp <- make_tradeoff_plots(SLA_CP, tradeoffs_SLA_CP,
                                 "SLA_ssp_avg", "CP_ssp_avg",
                                 "log SLA (m²/kgC)", "log C:P ratio")

p_sla_lt <- make_tradeoff_plots(SLA_LT, tradeoffs_SLA_LT,
                                 "SLA_ssp_avg", "LT_ssp_avg",
                                 "log SLA (m²/kgC)", "log leaf turnover (yr⁻¹)")

p_srl_rd <- make_tradeoff_plots(SRL_RootD, tradeoffs_SRL_RootD,
                                 "SRL", "Root_D",
                                 "log SRL (m/kgC)", "log root diameter (m)")

p_rcn_rt <- make_tradeoff_plots(RCN_RT, tradeoffs_RCN_RT,
                                 "Root_CN", "turnroot",
                                 "log root C:N", "log root turnover (yr⁻¹)")

p_rd_mcol <- make_tradeoff_plots(MRD_MCOL, tradeoffs_MRD_MCOL,
                                  "Root_D", "mcol",
                                  "Root diameter (m)", "Mycorrhizal colonisation",
                                  func = hdFunc_sigmoid)

p_lcn_lt <- make_tradeoff_plots(LCN_LT, tradeoffs_LCN_LT,
                                 "CN_ssp_avg", "LT_ssp_avg",
                                 "log leaf C:N", "log leaf turnover (yr⁻¹)")

# Combined tradeoff panel
graphics.off()
png(file.path(plot_dir, "trait_range_trade-offs.png"), width = 1500, height = 1900)
theme_set(theme_cowplot(font_family = "Times"))
plot_grid(
  p_sla_cn$global,  p_sla_cn$per_pft,
  p_sla_cp$global,  p_sla_cp$per_pft,
  p_sla_lt$global,  p_sla_lt$per_pft,
  p_rcn_rt$global,  p_rcn_rt$per_pft,
  p_srl_rd$global,  p_srl_rd$per_pft,
  p_rd_mcol$global, p_rd_mcol$per_pft,
  labels = c("a", "", "b", "", "c", "", "d", "", "e", "", "f", ""),
  label_size = 25, nrow = 6, rel_widths = c(0.5, 1)
)
graphics.off()

# Alternative panel with Leaf C:N ↔ Leaf turnover instead of SLA ↔ Leaf turnover
graphics.off()
png(file.path(plot_dir, "trait_range_trade-offs_lcnturn.png"), width = 1500, height = 1900)
theme_set(theme_cowplot(font_family = "Times"))
plot_grid(
  p_sla_cn$global,  p_sla_cn$per_pft,
  p_sla_cp$global,  p_sla_cp$per_pft,
  p_lcn_lt$global,  p_lcn_lt$per_pft,
  p_rcn_rt$global,  p_rcn_rt$per_pft,
  p_srl_rd$global,  p_srl_rd$per_pft,
  p_rd_mcol$global, p_rd_mcol$per_pft,
  labels = c("a", "", "b", "", "c", "", "d", "", "e", "", "f", ""),
  label_size = 25, nrow = 6, rel_widths = c(0.5, 1)
)
graphics.off()

cat("\nDone. Output files:\n")
cat("  ", file.path(output_dir, "trait_ranges_pfts_IQR.csv"), "\n")
cat("  ", file.path(output_dir, "tradeoffs_All_PFT.csv"), "\n")
cat("  ", file.path(plot_dir, "trait_range_trade-offs.png"), "\n")
cat("  ", file.path(plot_dir, "trait_range_trade-offs_lcnturn.png"), "\n")
