# =============================================================================
# Parameter Optimization for LPJ-GUESS-NTD Turnover Tradeoffs
# =============================================================================
#
# Optimizes turnover tradeoff parameters (leaf, root, sapwood) for a given PFT
# against remote sensing GPP observations using Differential Evolution (DEoptim).
#
# Method:
#   1. DEoptim proposes candidate parameter vectors
#   2. Each evaluation writes the parameters into the .ins file
#   3. LPJ-GUESS-NTD runs in parallel across all grid points (batchtools)
#   4. Model GPP is compared to observed GPP (sum of squared errors)
#   5. DEoptim iterates until convergence or max iterations
#
# Parameters optimized (6 total, 3 tradeoff pairs):
#   par[1] = k_turnover_leaf_a   (SLA ↔ leaf turnover, intercept)
#   par[2] = k_turnover_leaf_b   (SLA ↔ leaf turnover, slope)
#   par[3] = k_turnover_root_a   (root C:N ↔ root turnover, intercept)
#   par[4] = k_turnover_root_b   (root C:N ↔ root turnover, slope)
#   par[5] = k_turnover_sap_a    (WSG ↔ sapwood turnover, intercept)
#   par[6] = k_turnover_sap_b    (WSG ↔ sapwood turnover, slope)
#
# All tradeoffs use the log-linear form: log(Y) = h0 + h1 * log(X)
#
# Observation data:
#   Remote sensing GPP (e.g. from MODIS or FLUXCOM) at the grid point locations,
#   stored as a CSV with columns including "gpp" (in gC/m²/yr × 1000, i.e.
#   kgC/m²/yr needs dividing by 1000).
#
# Requirements:
#   R packages: DEoptim, data.table, batchtools
#   A compiled LPJ-GUESS-NTD binary
#   .ins files for the target PFT (cf .ins + global .ins)
#   A gridlist file with the calibration grid points
#
# Note:
#   This script was used for exploratory optimization and was not applied to
#   the final parameter set. It is provided as a reference for future studies.
#
# Part of: lpj-guess-ntd optimization pipeline
# =============================================================================


# =============================================================================
# 1. Setup
# =============================================================================

library(DEoptim)
library(data.table)
library(batchtools)

# --- Configurable paths (UPDATE THESE) ---
pft           <- "BE"                                       # PFT to optimize
data_dir      <- "~/VESTA/eo_data/"                         # Remote sensing observations
base_dir      <- paste0("~/VESTA/optimization/parallel/", pft, "/")
reg_dir       <- "~/VESTA/optimization/parallel/registry"   # batchtools registry
guess_binary  <- "/home/mdantas/stormbringer/LPJ-GUESS/build/guess"

# --- .ins file names (relative to base_dir) ---
ins_cf_file   <- paste0(pft, "_random_vesta_cf.ins")        # CF driver .ins
ins_global_file <- "global_ntd_vesta.ins"                   # Global .ins (parameters here)

# --- Gridlist ---
gridlist_file <- paste0("global_gridlist_random_opt_", pft, "_test.txt")

# --- Observation data ---
obs_csv       <- paste0(data_dir, pft, "_points_cmass_gpp_obs_30points.csv")

# --- Optimization settings ---
n_gridpoints  <- 30       # Number of calibration grid points
n_population  <- 30       # DEoptim population size (NP)
n_iterations  <- 20       # Maximum DEoptim iterations
n_cores       <- n_gridpoints  # Parallel cores for LPJ-GUESS runs


# =============================================================================
# 2. Load observation data and gridlist
# =============================================================================

obs_data  <- fread(obs_csv, header = TRUE)
grid_list <- readLines(file.path(base_dir, "gridlist", gridlist_file))

cat("Loaded", nrow(obs_data), "observation points\n")
cat("Gridlist has", length(grid_list), "entries\n")


# =============================================================================
# 3. LPJ-GUESS parallel run wrapper
# =============================================================================

#' Run LPJ-GUESS for a single grid point
#'
#' Called by batchtools for each grid point in parallel. Reads the current .ins
#' files (which NTD_prediction() has already updated with candidate parameters),
#' writes per-run copies, and executes the model.
#'
#' @param i Grid point index (1 to n_gridpoints)
run_single_point <- function(i) {
  run_name <- paste0("run", i)
  run_dir  <- file.path(base_dir, "result", run_name)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  # Read current .ins files
  ins_cf_lines     <- readLines(file.path(base_dir, ins_cf_file))
  ins_global_lines <- readLines(file.path(base_dir, ins_global_file))

  # Per-run file names
  run_ins_cf     <- paste0(run_name, "_cf.ins")
  run_ins_global <- paste0(run_name, "_global.ins")

  # Write per-run gridlist (single grid point)
  run_gridlist <- paste0(run_name, "_", gridlist_file)
  writeLines(grid_list[i], file.path(base_dir, "gridlist", run_gridlist))

  # Update cf .ins: point to per-run gridlist and global .ins
  ins_cf_run <- sub(gridlist_file, run_gridlist, ins_cf_lines)
  ins_cf_run <- sub('import\\s+"[^"]*"', paste0('import "', run_ins_global, '"'), ins_cf_run)

  # Write per-run .ins files
  writeLines(ins_cf_run, file.path(run_dir, run_ins_cf))
  writeLines(ins_global_lines, file.path(run_dir, run_ins_global))

  # Run LPJ-GUESS
  orig_dir <- getwd()
  setwd(run_dir)
  system2(guess_binary, args = paste0("-input cf ", run_ins_cf),
          stdout = FALSE, stderr = FALSE)
  setwd(orig_dir)
}


# =============================================================================
# 4. Prediction function (called by DEoptim at each iteration)
# =============================================================================

#' Run LPJ-GUESS-NTD with candidate parameters and return predictions
#'
#' Writes candidate parameters into the global .ins file, runs the model in
#' parallel across all grid points, and collects the predicted GPP.
#'
#' @param par Numeric vector of 6 parameters (see header)
#' @return Numeric vector of predicted GPP values (kgC/m²/yr)
NTD_prediction <- function(par) {

  # --- Write candidate parameters into the global .ins file ---
  ins_lines <- readLines(file.path(base_dir, ins_global_file))

  ins_lines <- sub("k_turnover_leaf_a .*", paste0("k_turnover_leaf_a ", par[1]), ins_lines)
  ins_lines <- sub("k_turnover_leaf_b .*", paste0("k_turnover_leaf_b ", par[2]), ins_lines)
  ins_lines <- sub("k_turnover_root_a .*", paste0("k_turnover_root_a ", par[3]), ins_lines)
  ins_lines <- sub("k_turnover_root_b .*", paste0("k_turnover_root_b ", par[4]), ins_lines)
  ins_lines <- sub("k_turnover_sap_a .*",  paste0("k_turnover_sap_a ",  par[5]), ins_lines)
  ins_lines <- sub("k_turnover_sap_b .*",  paste0("k_turnover_sap_b ",  par[6]), ins_lines)

  writeLines(ins_lines, file.path(base_dir, ins_global_file))

  # --- Run model in parallel via batchtools ---
  if (dir.exists(reg_dir)) {
    unlink(reg_dir, recursive = TRUE)
  }
  reg <- makeRegistry(file.dir = reg_dir, make.default = TRUE)
  reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = n_cores)

  batchMap(fun = run_single_point, i = seq_len(n_gridpoints))
  submitJobs()
  waitForJobs()

  # --- Collect GPP predictions ---
  pred <- numeric(0)
  for (j in seq_len(n_gridpoints)) {
    result_file <- file.path(base_dir, "result", paste0("run", j), "agpp.out")
    if (!file.exists(result_file)) {
      warning("Missing output for run ", j)
      return(rep(NaN, n_gridpoints))
    }
    result_agpp <- fread(result_file)
    agpp_mean   <- result_agpp[, .(agpp_mean = mean(Total)), by = c("Lon", "Lat")]
    pred <- c(pred, agpp_mean$agpp_mean)
  }

  return(pred)
}


# =============================================================================
# 5. Objective function
# =============================================================================

#' Compute sum of squared errors between model and observed GPP
#'
#' @param par Numeric vector of 6 parameters
#' @return Scalar SSE value (lower is better)
objective_function <- function(par) {
  pred <- NTD_prediction(par)

  # Handle failed runs
  if (any(is.nan(pred))) {
    warning("NaN in predictions — returning large penalty")
    return(.Machine$double.xmax / 2)
  }

  # Observed GPP (convert from gC/m²/yr to kgC/m²/yr to match model output)
  obs <- obs_data$gpp[seq_len(n_gridpoints)] / 1000

  # Sum of squared errors
  sse <- sum((pred - obs)^2)
  return(sse)
}


# =============================================================================
# 6. Parameter bounds and optimization
# =============================================================================

# Parameter bounds: [intercept_min, intercept_max] for a, [slope_min, slope_max] for b
# All tradeoffs are in log-log space: log(Y) = a + b * log(X), so slopes are negative
#
# Reference values from trait extraction (extract_trait_parameters.R):
#   k_turnover_leaf_a ~ 3.34,  k_turnover_leaf_b ~ -0.94
#   k_turnover_root_a ~ 2.98,  k_turnover_root_b ~ -0.78
#   k_turnover_sap_a  ~ 608,   k_turnover_sap_b  ~ -1.25

lower_bounds <- c(
  1.5,    # k_turnover_leaf_a
  -0.75,  # k_turnover_leaf_b
  1.0,    # k_turnover_root_a
  -0.5,   # k_turnover_root_b
  400,    # k_turnover_sap_a
  -2.0    # k_turnover_sap_b
)

upper_bounds <- c(
  2.5,    # k_turnover_leaf_a
  -0.175, # k_turnover_leaf_b
  2.0,    # k_turnover_root_a
  -0.2,   # k_turnover_root_b
  800,    # k_turnover_sap_a
  0.0     # k_turnover_sap_b
)

cat("\n=== Starting DEoptim optimization for PFT:", pft, "===\n")
cat("Parameters: k_turnover_leaf_a/b, k_turnover_root_a/b, k_turnover_sap_a/b\n")
cat("Population:", n_population, " Max iterations:", n_iterations, "\n")
cat("Grid points:", n_gridpoints, " Parallel cores:", n_cores, "\n\n")

result <- DEoptim(
  fn      = objective_function,
  lower   = lower_bounds,
  upper   = upper_bounds,
  DEoptim.control(
    NP      = n_population,
    itermax = n_iterations,
    trace   = TRUE,
    reltol  = 1e-6
  )
)


# =============================================================================
# 7. Results
# =============================================================================

best_params <- result$optim$bestmem
names(best_params) <- c("k_turnover_leaf_a", "k_turnover_leaf_b",
                         "k_turnover_root_a", "k_turnover_root_b",
                         "k_turnover_sap_a",  "k_turnover_sap_b")

cat("\n=== Optimization complete ===\n")
cat("Best SSE:", result$optim$bestval, "\n")
cat("Best parameters:\n")
print(best_params)

# Save results
save(result, file = file.path(base_dir, "deoptim_result.RData"))
fwrite(as.data.table(t(best_params)),
       file.path(base_dir, "best_parameters.csv"))

cat("\nResults saved to:", base_dir, "\n")
