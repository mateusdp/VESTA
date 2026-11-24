# HPC Workflow — MareNostrum 5

Scripts for running the global 1 km LPJ-GUESS-NTD campaign on MareNostrum 5 (BSC). The globe is divided into 11 latitude bands processed sequentially, each using 100 nodes (11,200 MPI processes) for up to 48 hours.

## Latitude bands

```
60N_90N  50N_60N  40N_50N  30N_40N  20N_30N
10N_20N  00N_10N  10S_00N  20S_10S  30S_20S  60S_30S
```

Bands are processed one at a time due to storage constraints on `/gpfs/scratch/`. After each band finishes, output is gathered, compressed, and downloaded before the next band begins.

## Per-band workflow

```
1. bash setup_band.sh 60N_90N          # prepare .ins files + submit script
2. cd <run_dir> && bash submit.sh      # build guess.cmd, submit to SLURM
3. bash progress.sh                    # monitor (optional)
4. bash gather_compressed.sh -d <remote_run_dir>   # download + concatenate output
5. Repeat for the next band
```

Step 4 runs from the local machine (not on MN5). With `COMPRESS_OUTPUT` enabled, each MPI process writes `.out.gz` + `.out.hdr` files directly, and `gather_compressed.sh` concatenates them over sshfs. If `COMPRESS_OUTPUT` is not enabled, use `join_compress.sh` on-cluster instead.

## Scripts

| File | Description |
|------|-------------|
| `setup_band.sh` | **Step 1.** Creates a band-specific directory, copies the `.ins` template files, patches all band-dependent paths (gridlist, climate NetCDFs, soil map), and copies the submission script with the correct configuration. |
| `submit_chelsa_template.sh` | **Step 2.** Copied into each band's run directory as `submit.sh` by `setup_band.sh`. Builds a combined `.ins` file (flattening imports), fills `PLACEHOLDER_*` tokens in `guess.cmd.template`, and submits the job via `sbatch`. |
| `guess.cmd.template` | SLURM batch script template. At job start: builds the combined `.ins`, creates 11,200 run directories, splits the gridlist round-robin across processes (awk), splits the soil file (GNU `split -n r/N`), then launches LPJ-GUESS with `mpiexec`. |
| `recursive-ins-compilator.sh` | Recursively resolves `import` directives in LPJ-GUESS `.ins` files into a single combined file. Required because MPI processes start from different directories. Original by Johan Nord (2019). |
| `gather_compressed.sh` | **Post-job (production).** Runs from a local machine (e.g. mrpharmacist). Mounts the remote run directory via sshfs, concatenates per-process `.out.gz` streams into single gzipped files, and optionally cleans up. Used when LPJ-GUESS is compiled with `COMPRESS_OUTPUT` (the production configuration). |
| `join_compress.sh` | **Post-job (alternative).** Runs on-cluster. Appends per-process `.out` files and compresses with gzip. Used when `COMPRESS_OUTPUT` is *not* enabled at compile time. |
| `progress.sh` | Samples `guess.log` from evenly-spaced run directories to check simulation progress. Reports how many processes have finished. |

## Template `.ins` files

```
template/
├── global_p-traits_daily_cf.ins         # Main instruction file (climate paths, gridlist, soil)
├── global_p-traits-srlmyco_daily.ins    # PFT definitions, trait variation, run parameters
└── global_soiln.ins                     # Soil nitrogen cycling parameters
```

The `.ins` files use the default band `60S_30S` as a placeholder. `setup_band.sh` patches the band-specific references via `sed` before each run.

### Key simulation parameters (from `.ins`)

| Parameter | Value |
|-----------|-------|
| Vegetation mode | Cohort |
| Spinup | 500 years |
| Simulation period | 1979–2021 (43 years CHELSA v2.1 climate) |
| Patches per cell | 15 |
| Fire model | GLOBFIRM |
| Trait variation | SLA, wood specific gravity, specific root length |
| Nutrient limitation | Nitrogen + phosphorus |
| Input module | `cf` (climate forcing, NetCDF) |

### PFTs defined

BNE, BINE, BNS, TeNE, TeBS, IBS, TeBE, TrBE, TrIBE, TrBR, C3G, C4G

## MareNostrum 5 configuration

| Parameter | Value |
|-----------|-------|
| Nodes per band | 100 |
| Cores per node | 112 |
| Total MPI processes | 11,200 |
| Wall time | 48 hours |
| Partition | `gpp` (General Purpose) |
| Account | `ehpc336` (EuroHPC allocation) |
| QoS | `gp_ehpc` |

## Directory structure on MN5

```
/gpfs/scratch/ehpc336/LPJ_GUESS/
├── input/chelsa_bands/                     # Preprocessed inputs
│   ├── gridlist_{BAND}.txt                 # Per-band gridlists
│   ├── co2_1850-2100_sres_A1B.txt          # CO₂ concentrations
│   ├── GlobalNitrogenDeposition.bin        # N deposition
│   ├── GlobalTSPPO4depLandMonthly.csv      # P deposition
│   ├── lith_classes_shield_pdata_05.csv    # P weathering
│   ├── SimfireInput.bin                    # SIMFIRE data
│   └── {BAND}/
│       ├── temp_{BAND}.nc
│       ├── prec_{BAND}.nc
│       ├── rad_{BAND}.nc
│       └── soilmap_{BAND}.dat
├── ins/
│   ├── global_chelsa_v2_template/          # Template .ins files (this directory)
│   └── global_chelsa_v2_{BAND}/            # Created by setup_band.sh
└── scripts/
    ├── setup_band.sh
    ├── submit_chelsa_template.sh
    ├── guess.cmd.template
    └── recursive-ins-compilator.sh
```

## Ancillary input files

These files are **not** produced by the preprocessing pipeline and must be obtained separately:

| File | Description |
|------|-------------|
| `co2_1850-2100_sres_A1B.txt` | CO₂ concentrations (SRES A1B scenario, 1850–2100) |
| `GlobalNitrogenDeposition.bin` | Nitrogen deposition (LPJ-GUESS standard binary) |
| `GlobalTSPPO4depLandMonthly.csv` | Phosphorus deposition (monthly, global) |
| `lith_classes_shield_pdata_05.csv` | Phosphorus weathering by lithology class |
| `SimfireInput.bin` | SIMFIRE data (loaded but not used — fire model is GLOBFIRM) |

## Building LPJ-GUESS-NTD on MN5

The model is compiled with GCC + OpenMPI and linked against parallel NetCDF/HDF5. Load the following modules before running `cmake` and `make`:

```bash
module unload impi/2021.10.0
module load cmake
module load bsc/1.0
module unload intel

module load gcc/12.3.0
module load openmpi/4.1.5-gcc
module load hdf5/1.14.1-2-gcc-openmpi
module load pnetcdf/1.12.3-gcc-openmpi
module load netcdf/2023-06-14-gcc-openmpi
```

An Intel toolchain (`intel/2023.1` + `openmpi/4.1.5` + matching NetCDF/HDF5/pnetcdf) also works but is not used in the production runs.

Then configure and build:

```bash
mkdir build && cd build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-D COMPRESS_OUTPUT" \
    -DCMAKE_C_FLAGS="-I${HDF5_INCLUDE}" \
    -DCMAKE_EXE_LINKER_FLAGS="-L${HDF5_LIB}" \
    -DUNIT_TESTS=OFF \
    /path/to/LPJ-GUESS/source
make -j16
```

Key build options used in the production runs:

| Option | Value | Notes |
|--------|-------|-------|
| Build type | `Release` (`-O3`) | |
| `COMPRESS_OUTPUT` | Enabled (`-D COMPRESS_OUTPUT`) | Gzipped `.out` files from each process |
| MPI | Auto-detected via OpenMPI 4.1.5 (GCC) | `mpicxx` / `mpicc` from module |
| NetCDF | `c-4.9.2` / `cxx4-4.3.1` with HDF5 1.14.1-2 and pnetcdf 1.12.3 | GCC + OpenMPI variant |
| Unit tests | OFF | |

The HDF5 include/library paths (`-I` and `-L` flags) are needed because the MN5 module system does not always set them automatically. With the modules listed above, these resolve to `/apps/GPP/HDF5/1.14.1-2/GCC/OPENMPI/{include,lib}`.

The compiled binary is expected at the path set by `GUESS_BIN` in `submit_chelsa_template.sh`.

## Adapting for another HPC system

1. **Update paths** — all configurable paths are grouped at the top of `setup_band.sh` and `submit_chelsa_template.sh`. The `.ins` template files reference `/gpfs/scratch/ehpc336/` for ancillary inputs; update these or use `setup_band.sh` to patch them.

2. **Update SLURM parameters** — `NPROCESS`, `WALLTIME`, `PARTITION`, `ACCOUNT`, and `QOS` are set in `submit_chelsa_template.sh`. The `guess.cmd.template` uses `PLACEHOLDER_*` tokens that are filled at submission time.

3. **Update email** — `guess.cmd.template` contains a hardcoded `--mail-user` address.

4. **Update GUESS binary path** — set `GUESS_BIN` in `submit_chelsa_template.sh`.

5. **Check GNU coreutils** — the soil file split uses `split -n r/N` (round-robin by line), a GNU extension. This works on Linux but not on macOS/BSD.

## Notes

- The gridlist is split round-robin across processes (awk-based, ~600× faster than a bash loop for ~8.6M cells), ensuring even distribution of high- and low-productivity cells across processes.
- The combined `.ins` file is built inside the SLURM job. If the import flattening fails, the allocation is wasted. Consider moving this step to submission time for safety.
- SIMFIRE input is configured in the `.ins` file but the fire model is set to GLOBFIRM, so the SIMFIRE data is loaded but not used. This is harmless but potentially confusing.
- January 1979 climate data comes from CHELSA v1.2 (v2.1 starts February 1979). See the preprocessing documentation for details.
