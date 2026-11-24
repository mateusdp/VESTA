# LPJ-GUESS-NTD: Global 1 km Simulation Pipeline for the VESTA (VEgetation Spatialization of Traits Algorithm) project

Preprocessing, HPC workflow, and postprocessing scripts for the first global
1 km resolution simulation with LPJ-GUESS-NTD (Nutrient-Trait Dynamics),
run on MareNostrum 5 at the Barcelona Supercomputing Center (BSC).

This repository accompanies the ESA Living Planet Fellowship project and
provides everything needed to reproduce the global campaign, from raw input
data to validated output maps.

## Model overview

LPJ-GUESS-NTD is a dynamic global vegetation model that represents
continuous trait variation within plant functional types (PFTs). Instead of
assigning fixed trait values, each cohort draws traits — specific leaf area
(SLA), wood specific gravity (WSG), specific root length (SRL), and
stoichiometric ratios — from PFT-specific distributions, capturing
within-PFT functional diversity.

The global campaign simulates ~239 million land cells at 30 arc-second
(~1 km) resolution, covering 60°S to 90°N, with 12 PFTs and 15 replicate
patches per cell. The model is forced by CHELSA v2.1 monthly climate
(1979–2021) and runs with coupled nitrogen and phosphorus cycling.

## Repository structure

```
├── README.md                          # This file
├── preprocessing/
│   ├── chelsa/
│   │   ├── download_chelsa.sh         # Download CHELSA v2.1 climate NetCDFs
│   │   └── process_chelsa_to_guess.R  # Convert to LPJ-GUESS cf-module format
│   ├── soil/
│   │   └── process_hwsd_to_guess.R    # HWSD v1.21 soil data → per-band .dat
│   └── traits/
│       └── extract_trait_parameters.R # Extract trait ranges from TRY database
├── hpc/
│   └── mn5/
│       ├── README.md                  # Detailed HPC workflow documentation
│       ├── setup_band.sh             # Prepare band-specific run directory
│       ├── submit_chelsa_template.sh # Build .ins, submit SLURM job
│       ├── guess.cmd.template        # SLURM batch script template
│       ├── recursive-ins-compilator.sh # Flatten .ins import directives
│       ├── gather_compressed.sh      # Download + concatenate gzipped output
│       ├── join_compress.sh          # On-cluster join (non-compressed builds)
│       ├── progress.sh               # Monitor simulation progress
│       └── template/                 # LPJ-GUESS instruction files (.ins)
│           ├── global_p-traits_daily_cf.ins
│           ├── global_p-traits-srlmyco_daily.ins
│           └── global_soiln.ins
├── optimization/
│   └── optimize_turnover_deoptim.R   # DEoptim trait turnover optimisation
├── LPJ-GUESS-NTD/
│   └── code/   # NTD model code
└── postprocessing/
    ├── plot_global.R                  # Read bands → GeoTIFF + global/regional PNGs
    ├── validate_gpp_gosif.R           # GPP validation at 0.5° vs GOSIF v4.0
    ├── validate_gpp_modis_regional.R  # GPP validation at 1 km vs MODIS MOD17
    ├── validate_traits_gbif.R         # Trait validation vs GBIF-TRY gapfilled maps
    └── validate_traits_splot.R        # Trait validation vs sPlotOpen CWMs
```

## Pipeline overview

The simulation pipeline has four stages, each with its own subdirectory.

### 1. Preprocessing

**Climate** (`preprocessing/chelsa/`): Downloads CHELSA v2.1 monthly
temperature, precipitation, and radiation NetCDFs (1979–2021) and converts
them to the LPJ-GUESS `cf` input module format, with per-band subsetting
and proper unit conversions. January 1979 is sourced from CHELSA v1.2 (v2.1
begins February 1979).

**Soil** (`preprocessing/soil/`): Processes the Harmonized World Soil
Database (HWSD v1.21) raster and attribute table into per-band soil
parameter files for LPJ-GUESS, mapping FAO texture classes to sand/clay/silt
fractions.

**Traits** (`preprocessing/traits/`): Extracts PFT-level trait distribution
parameters (mean, SD, min, max) from TRY database exports for SLA, WSG,
and SRL. These feed the `.ins` template files that define PFT trait ranges.

### 2. HPC workflow

**MareNostrum 5** (`hpc/mn5/`): The globe is split into 11 latitude bands
processed sequentially. Each band runs on 100 nodes (11,200 MPI processes)
for up to 48 hours. The workflow handles directory setup, `.ins` file
preparation, SLURM submission, progress monitoring, and output gathering.

See [`hpc/mn5/README.md`](hpc/mn5/README.md) for the full workflow
documentation, including build instructions, module configuration, and
guidance for adapting to other HPC systems.

### 3. Optimisation (optional)

**Parameter optimisation** (`optimization/`): Uses DEoptim to optimise trait
turnover parameters against remote sensing targets (e.g. biomass maps).
While this was not used to calibrate the production run, it provides a
reusable framework for future sensitivity analyses.

### 4. Postprocessing

**Map production** (`postprocessing/plot_global.R`): Reads per-band
`.out.gz` files, computes temporal means, and produces 30 arc-second
GeoTIFFs (plain and DEFLATE-compressed) plus publication-quality PNG maps at
global and regional scales (34 predefined regions).

**GPP validation at 0.5°** (`postprocessing/validate_gpp_gosif.R`):
Aggregates modelled annual GPP to 0.5° and compares against GOSIF SIF-based
GPP (Li & Xiao, 2019). Produces scatter plots, bias maps, and zonal mean
profiles with per-latitude-band metrics.

**GPP validation at 1 km** (`postprocessing/validate_gpp_modis_regional.R`):
Compares modelled GPP at native 1 km resolution against MODIS MOD17A3HGF
in selected regions with strong topographic or land-cover gradients,
demonstrating the added value of the 1 km resolution. Produces side-by-side
maps, density scatters, and E–W transects.

**Trait validation vs GBIF-TRY**
(`postprocessing/validate_traits_gbif.R`): Validates community-weighted
mean SLA, WSG, and SRL at 0.5° against gridded trait maps derived from
GBIF occurrences and TRY gapfilled trait values.

**Trait validation vs sPlotOpen**
(`postprocessing/validate_traits_splot.R`): Validates community-weighted
mean SLA, WSG, Leaf C:N, and Leaf C:P at plot locations (~95K plots)
against sPlotOpen vegetation-plot data.

## Simulation output

Each latitude band produces ~22 gzipped output files (`*.out.gz`). Key
output variables include:

| File | Description |
|------|-------------|
| `cmass.out.gz` | Vegetation carbon mass (kg C m⁻²) |
| `agpp.out.gz` | Annual gross primary production (kg C m⁻² yr⁻¹) |
| `cwm_traits.out.gz` | Community-weighted mean traits (SLA, WSG, CN, CP, …) |
| `sla.out.gz` | Community-weighted mean SLA per PFT |
| `wsg.out.gz` | Community-weighted mean WSG per PFT |
| `srl.out.gz` | Community-weighted mean SRL per PFT |

Total compressed output is approximately 4 TB across the 11 bands.

## Requirements

### Preprocessing and postprocessing (R)

- R ≥ 4.1
- Core packages: `terra`, `data.table`
- Plotting: `ggplot2`, `tidyterra`, `patchwork`, `scales`

### HPC

- LPJ-GUESS-NTD source code (compiled with MPI + NetCDF support)
- SLURM workload manager
- GNU coreutils (`split -n r/N` for round-robin gridlist splitting)

## Input data

| Dataset | Version | Resolution | Use |
|---------|---------|------------|-----|
| CHELSA | v2.1 (v1.2 for Jan 1979) | 30 arc-sec | Monthly temperature, precipitation, radiation |
| HWSD | v1.21 | 30 arc-sec | Soil texture and parameters |
| TRY | — | — | Trait distribution parameters |
| CO₂ concentrations | SRES A1B | — | Atmospheric CO₂ forcing |
| N deposition | LPJ-GUESS standard | 0.5° | Nitrogen deposition |
| P deposition / weathering | — | 0.5° | Phosphorus cycling |

## Validation data (not included)

| Dataset | Resolution | Script |
|---------|------------|--------|
| GOSIF GPP v4.0 | 0.05° | `validate_gpp_gosif.R` |
| MODIS MOD17A3HGF | 500 m | `validate_gpp_modis_regional.R` |
| GBIF-TRY gapfilled trait maps | 0.5° | `validate_traits_gbif.R` |
| sPlotOpen CWM/CWV | Plot-level | `validate_traits_splot.R` |

## Usage

All scripts use configurable paths at the top of each file — update the
`CONFIG` section before running. No script uses `setwd()` or
RStudio-specific calls, so they can be run non-interactively:

```bash
Rscript postprocessing/plot_global.R
Rscript postprocessing/validate_gpp_gosif.R
```

## References

- Smith, B. et al. (2014). Implications of incorporating N cycling and N
  limitations on primary production in an individual-based dynamic
  vegetation model. *Biogeosciences*, 11, 2027–2054.
- Karger, D.N. et al. (2017). Climatologies at high resolution for the
  earth's land surface areas. *Scientific Data*, 4, 170122.
- Li, X. & Xiao, J. (2019). A global, 0.05-degree product of
  solar-induced chlorophyll fluorescence derived from OCO-2, MODIS, and
  reanalysis data. *Remote Sensing*, 11, 517.
- Bruelheide, H. et al. (2019). sPlotOpen — An environmentally balanced,
  open-access, global dataset of vegetation plots. *Global Ecology and
  Biogeography*, 28, 1839–1855.

## Licence

LPJ-GUESS-NTD source code and these workflow scripts are provided for
reproducibility of the ESA Living Planet Fellowship project results.

## Contact

Mateus Dantas de Paula — Goethe-Universität Frankfurt / Senckenberg
Biodiversity and Climate Research Centre (SBiK-F)
