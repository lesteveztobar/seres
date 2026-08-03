      # canopymicroenv

> **Microclimate niche modelling and vertical colonization simulation for epiphytic orchids — NW Ecuador**
> Master's thesis project — University of Bonn, 2026

---

## What this is

This repository contains the full analysis pipeline for my master's thesis on the **vertical stratification of epiphytic Maxillariinae orchids** and their microenvironmental correlates across seven cloud forest sites in the Chocó Andino of northwestern Ecuador (Maquipucuna, Mashpi, MindoTarabita, La Elenita, Mindo Mirador, Saloya, Yanayacu). The original five-site dataset included a site labelled "MiradorMindo," which was later found to bundle three physically distinct locations under one name; it was split into La Elenita, Mindo Mirador, and Saloya (see `scripts/00_data_conversion/rebuild_combined_csv.py`).

The core idea: instead of using coarse climate data to describe orchid habitat, this pipeline models the **exact microclimate at the height and location where each individual was observed** — at 0.25 m resolution across the full canopy vertical gradient (justified against finer/coarser alternatives; see [Notes](#notes)). These microclimate profiles characterise the realised niche of each species and feed a **3D spatially explicit colonization model** that simulates population dynamics (dispersal, establishment, survival, growth) across the canopy landscape, driven by microclimate-based vital rates.

> ⚠️ **Work in progress.** Microclimate pipeline complete for all 7 sites; colonization model and sensitivity experiments running. See [Status](#status).

---

## Study system

- **Focal group:** Maxillariinae orchids (tribe Maxillarieae, subtribe Maxillariinae)
- **Sites:** 7 cloud forest sites, ~800–1800 m elevation, NW Ecuador
- **Observations:** 141 individuals across all sites, with observed canopy height, elevation, genus, and photo references
- **Approach:** Field observations → microclimate modelling → niche characterisation → 3D vertical colonization simulation → parameter sensitivity experiments

---

## Pipeline overview

```
GeoJSON field exports
        ↓
  scripts/00_data_conversion/convert_observations.py     # parse field notes → structured CSV
        ↓
  geojson_to_csv/csv/                 # per-site CSVs (one per field site)
        ↓
  data/csv/combined_with_identification.csv  # merged dataset, all 7 sites (combinedv3.csv remains
                                      # available as a manual override for curation-quality work; see Notes)
        ↓
  scripts/01_microclimate/run_microclimate_site.R     # HPC: single site via Rscript / SLURM
    └── scripts/01_microclimate/lib.R                 # ERA5/DTM/LAI/albedo/vegetation/soil acquisition +
                                      #   niche extraction, canopy grid, climate lookups, runpointmodela()
                                      # height loop: production 0.25 m steps, per-height RDS to scratch
                                      # height ceiling: measured CanopyHeight_m → vhgt.tif p99 → hObs_max
  microenv_<site>[_h<step>].rds       # manifest → data/processed/ (~1.4 MB)
  /lustre/scratch/.../microenv_<site>[_h<step>]_heights/  # per-height spatial arrays
        ↓
  scripts/02_model/setup/characterize_niches.R  # pools each species' climate niche across every
                                      # site it was observed at → data/processed/species_niches.rds
        ↓
  scripts/02_model/run/run_colonization_onesite.R  # single-site colonization run (all-sites: batch_exp.sh, SLURM array)
    └── scripts/02_model/engine/get_colonization.R  # build_forest() · dispersal · establishment
                                      # survival/growth (IPM-style, equation primitives)
        ↓
  scripts/02_model/setup/make_params.R  # build parameter sweep RDS files for experiments
                                      # (incl. best_case.rds — persistence validation, see Notes)
  scripts/02_model/run/batch_exp.sh / run_colonization.sh  # SLURM: 7 sites × experiments in parallel
  scripts/simple_model/simple_colonization.R  # standalone 3D model without microclimate
  scripts/simple_model/simple_experiments.R   # sensitivity experiments for simple model
  scripts/simple_model/run_extinction_heatmap.R  # fine-scale extinction threshold scan
  scripts/02_model/plots/plot_all.R / plot_functions.R  # all project figures in one pass

Resolution justification (before committing to full experiment runs):
  scripts/02_model/resolution/height_res_array.sh              # generate coarser height-step microenv variants
  scripts/02_model/resolution/run_resolution_diagnostics.sh     # timing only: cache-build + short run per resolution
  scripts/02_model/resolution/run_height_resolution_experiment.sh  # outcomes: full best_case runs per height step
```

---

## Colonization model

The colonization model simulates a 3D canopy landscape as a voxel grid `[x, y, z]` where z is height in the canopy. Each run goes through:

1. **Forest structure** (`build_forest`): stochastic tree placement from forest inventory parameters (stem density, crown radius, height distribution). Each tree's trunk and crown are labelled with Johansson (1974) zones (JZ1–JZ5), and the maximum epiphyte carrying capacity of each voxel (`carCap_voxel`) is derived from available bark surface area.

2. **Dispersal** (pass 1): adults produce seeds; each seed travels an exponentially distributed distance in the downwind direction. Dispersal is vectorized using `tabulate()` over a padded array.

3. **Establishment** (pass 2): seeds become seedlings if the voxel is within the canopy, below carrying capacity, and a germination event occurs. Germination probability is scaled by local relative humidity and shortwave radiation from `get_clim()`.

4. **Survival and growth** (pass 3): three stage classes — seedlings (S, 0–1 cm pseudobulb), juveniles (J, 1–7 cm), adults (A, 7–20 cm). Survival is computed via `survival_logit()` — logistic in size and microclimate (temperature, RH, light). Stage transitions use `transition_logit()` driven by annual precipitation and RH. Pseudobulb growth uses `size_increment()` with a cost-of-reproduction penalty for fruiting individuals. Each function names its literature source directly in the code.

All inner loops run at height-slice level (not voxel level), replacing `rbinom(1,...)` per voxel with `rbinom(n,...)` per height slice. Experiments are parallelized across parameter values via `parallel::mclapply`.

---

## HPC usage (Marvin cluster, University of Bonn)

The microclimate model is computationally intensive. Each site is submitted as an independent SLURM job via `run_microenv.sh`, coordinated by `microenv_array.sh`:

```bash
# From /home/s38leste_hpc/canopymicroenv/
sbatch scripts/01_microclimate/microenv_array.sh 12 0.25          # all 7 sites × 12 months of ERA5, production 0.25 m
sbatch scripts/01_microclimate/run_microenv.sh Maquipucuna 12 0.25  # single site

# Check microenv regeneration progress (read-only)
Rscript scripts/01_microclimate/check_microenv_progress.R Maquipucuna 0.1,0.25,0.5,1.0

# Pool each species' climate niche across every site it was observed at
sbatch scripts/02_model/setup/characterize_niches.sh

# Resolution justification, before committing to full experiment runs
sbatch scripts/02_model/resolution/run_resolution_diagnostics.sh Maquipucuna        # timing only
sbatch scripts/02_model/resolution/run_height_resolution_experiment.sh Maquipucuna  # outcomes (best_case.rds)

# Sensitivity experiments (after microclimate is done)
sbatch scripts/02_model/run/batch_exp.sh    # 7 sites × experiments
```

Each microclimate job runs on the `lm_short` partition (large-memory nodes) with 4 CPUs and 500 GB RAM. It loads R/4.4.2, the Miniforge3 conda environment (`canopy_rgee`) for Earth Engine access, and allocates a Lustre scratch workspace for per-height temp files via `ws_allocate`. Logs are written to `logs/log_<jobid>.out`.

**Lustre scratch workspace:** per-height `.rds` files are written to `/lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_<site>[_h<step>]_heights/` during computation. The final `microenv_<site>[_h<step>].rds` in `data/processed/` is a small manifest (~1.4 MB) that records the height vector, scratch path, and baseline weather — not the full spatial arrays. **Do not release the scratch workspace** (`ws_release`) until downstream analysis is complete. The workspace expires in 90 days and can be extended up to 3 times.

---

## Repository structure

Scripts are organized by pipeline stage. Each of `00_data_conversion/`,
`01_microclimate/`, and `03_orchestration/` has a `run.sh` dispatcher — a
single documented entry point (`run.sh help`) that shells out to the
underlying scripts, which remain separate, independently runnable files
(some carry their own `#SBATCH` resource headers and can't be merged into a
shared script). `simple_model/` has the equivalent `run.R`.

```
canopymicroenv/
├── scripts/
│   │   # ── Stage 0: data conversion (GeoJSON → CSV) ──────────────────────
│   ├── 00_data_conversion/
│   │   ├── run.sh                     # dispatcher: convert / photos / rebuild
│   │   ├── convert_observations.py    # parse iNaturalist / field GeoJSON → structured CSV
│   │   ├── build_photo_lookup.py      # build photo-lookup table from field exports
│   │   ├── rebuild_combined_csv.py    # LaElenita/MindoMirador/Saloya migration + combined.csv rebuild
│   │   └── helper_functions.R         # COMBINED_COL_TYPES: shared read_csv() column spec
│   │
│   │   # ── Stage 1: microclimate pipeline ──────────────────────────────
│   ├── 01_microclimate/
│   │   ├── run_microclimate.sh        # dispatcher: site / progress / fix-dtm /
│   │   │                              #   submit-site / submit-array
│   │   ├── run_microclimate_site.R    # data acquisition + point model + grid model, single site — called by SLURM
│   │   ├── lib.R                      # ERA5/DTM/LAI/albedo/vegetation/soil acquisition +
│   │   │                              #   niche extraction, canopy grid, climate lookups, get_clim()
│   │   ├── check_microenv_progress.R  # read-only: per-site/height-step regen progress
│   │   ├── regenerate_missing_dtm.R   # repair utility: re-fetch a site's missing dtm.tif
│   │   ├── run_microenv.sh            # SLURM: single-site microclimate job (#SBATCH header)
│   │   └── microenv_array.sh          # SLURM: submits run_microenv.sh per site (#SBATCH header)
│   │
│   │   # ── Complex (microclimate-driven) colonization model ──────────────
│   ├── 02_model/
│   │   ├── run.R                      # dispatcher: params / isolation-files / niches / onesite /
│   │   │                              #   resolution-diagnostics / resolution-experiment /
│   │   │                              #   climate-variation / competition / summarize /
│   │   │                              #   check-niche / check-transitions / check-run / plots
│   │   ├── run.sh                     # dispatcher: submit-* (sbatch pass-through, each script
│   │   │                              #   keeps its own #SBATCH header) / progress / check-run
│   │   │
│   │   ├── engine/
│   │   │   └── get_colonization.R     # colonization functions: build_forest(),
│   │   │                              #   runcolonization(), pass1/2/3 sub-models,
│   │   │                              #   survival_logit(), transition_logit(), size_increment(),
│   │   │                              #   load_height(), run_experiment()/run_factorial_experiment()/
│   │   │                              #   run_replicated(), plot_abundance()/plot_3d_abundance_animated()
│   │   │                              #   (⚠️ animated plot not working yet — see Notes)
│   │   ├── config/
│   │   │   ├── paths.R                # path constants, load_observations()
│   │   │   └── patches.R              # runtime monkey-patches (ecmwfr/microclimdata/microclimf)
│   │   ├── setup/
│   │   │   ├── make_params.R          # build parameter sweep RDS files, incl. best_case.rds /
│   │   │   │                          #   realistic.rds (persistence validation — see Notes)
│   │   │   ├── make_isolation_species_files.R  # per-(site,species) isolation params + manifest
│   │   │   └── characterize_niches.R/.sh  # pools each species' climate niche across every site
│   │   │                              #   it was observed at → data/processed/species_niches.rds
│   │   ├── run/
│   │   │   ├── run_colonization_onesite.R  # colonization run: single site (interactive / SLURM)
│   │   │   ├── run_colonization.sh    # SLURM: one site x one experiment
│   │   │   └── batch_exp.sh           # SLURM: sites × experiments in parallel
│   │   ├── resolution/                # resolution justification, before committing to full runs
│   │   │   ├── resolution_diagnostics.R/run_resolution_diagnostics.sh  # timing-only:
│   │   │   │                          #   climate-cache build + short run cost
│   │   │   ├── height_resolution_experiment.R/run_height_resolution_experiment.sh  # outcomes:
│   │   │   │                          #   full best_case.rds runs at each height step
│   │   │   └── height_res_array.sh    # generate coarser height-step microenv variants
│   │   ├── diagnostics/               # read-only unless noted
│   │   │   ├── check_niche_suitability.R
│   │   │   ├── check_transition_rates.R
│   │   │   ├── check_colonization_run.R
│   │   │   ├── check_colonization_progress.sh
│   │   │   └── run_check_colonization_run.sh
│   │   ├── analysis/
│   │   │   ├── climate_variation_test.R  # vertical/elevational climate variation tests
│   │   │   │                          #   (incl. elevation_helpers.R, merged in here 2026-07-28)
│   │   │   ├── competition_analysis.R # isolation-vs-multi-species comparison
│   │   │   └── summarize_all_results.R
│   │   └── plots/
│   │       ├── plot_functions.R       # all plotting functions
│   │       ├── plot_all.R             # driver: every project figure in one pass
│   │       └── run_plots.sh           # SLURM wrapper for plot_all.R
│   │
│   │   # ── Simple / standalone model (no microclimate) ─────────────────────
│   ├── simple_model/
│   │   ├── run.R                      # dispatcher: baseline / experiments / heatmap
│   │   ├── simple_colonization.R      # standalone 3D colonization model (no microclimate)
│   │   ├── simple_experiments.R       # simple model sensitivity experiments + animations
│   │   └── run_extinction_heatmap.R   # extinction threshold scan (p_est × repro_rate)
│   │
│   │   # ── Stage 2: top-level pipeline orchestration ──────────────────────
│   ├── 03_orchestration/
│   │   ├── run.sh                     # dispatcher: pipeline / full-analysis /
│   │   │                              #   factorial-remaining / persistence / cleanup-scratch
│   │   ├── run_pipeline.sh            # menu-driven launcher: chains microenv/params/niche/
│   │   │                              #   experiments/plots via SLURM job dependencies
│   │   ├── run_full_analysis_pipeline.sh  # full automated run (resolution, factorial,
│   │   │                              #   competition, climate variation, plots)
│   │   ├── run_factorial_remaining_sites.sh  # catch-up: reproduction factorial for
│   │   │                              #   whichever sites are still missing it
│   │   ├── run_persistence_all_sites.sh  # best_case.rds / realistic.rds validation, all sites
│   │   └── cleanup_resolution_scratch.sh  # reclaim Lustre scratch (dry-run unless --yes)
│   │
│   │   # ── Legacy (not part of the active pipeline) ────────────────────────
│   ├── legacy/
│   │   └── allsites.R                 # earlier all-sites driver; predates the site split
│   │                                  #   and the 02_model API — kept for reference only
│   │
│   └── geojson-csv-sql-conversion-tools/  # Node.js + Python toolkit for
│                                          #   GeoJSON ↔ CSV ↔ SQL round-trips;
│                                          #   used for QA and manual corrections
│                                          #   during the observation cleaning step.
│                                          #   © Rudo Kemper, GPL-3.0
│                                          #   github.com/rudokemper/geojson-csv-sql-conversion-tools
├── geojson_to_csv/
│   ├── raw/                           # original GeoJSON exports from iNaturalist / OrganicMaps app
│   └── csv/                           # per-site CSVs produced by convert_observations.py
├── data/
│   ├── csv/                           # combined*.csv, Processed*.csv — observation datasets
│   ├── literature/                    # literature-sourced trait/reference data
│   ├── raw/                           # ERA5, DTM, LAI, albedo, etc. (not tracked)
│   ├── params/                        # parameter sweep RDS files built by make_params.R,
│   │                                  #   incl. best_case.rds (persistence validation)
│   └── processed/                     # microenv_<site>[_h<step>].rds (manifests, ~1.4 MB each),
│                                      #   species_niches.rds, pointmodel_*.rds, colonization
│                                      #   outputs (not tracked)
│                                      #   NB: full per-height arrays live in Lustre scratch
├── output/                            # figures, animations, resolution-diagnostic CSVs (not tracked)
└── logs/                              # timestamped run logs (not tracked)
```

---

## Status

| Step | Status |
|---|---|
| Field data collection (7 sites) | ✅ Complete |
| GeoJSON → CSV conversion | ✅ Complete |
| Observation cleaning + `hObs` extraction | ✅ Complete |
| ERA5 climate data download | ✅ Complete |
| DTM, landcover, vegetation, soil parameters | ✅ Complete |
| Point model height loop (`runpointmodela`, all 7 sites) | ✅ Complete |
| Grid microclimate model (`runmicro`, all 7 sites, production 0.25 m) | ✅ Complete (height files in Lustre scratch) |
| Resolution justification (timing + outcome comparison) | 🔄 In progress |
| Cross-site species niche characterization | ✅ Complete |
| Colonization model — single site | 🔄 In progress |
| Colonization model — all sites | 🔄 In progress |
| Forest structure (`build_forest`, Myster 2017 params) | ✅ Implemented |
| IPM equation primitives (`survival_logit` etc.) | ✅ Implemented |
| Parameter sensitivity experiments | ✅ Implemented |
| Extinction threshold heatmap | ✅ Complete |
| Species identification | 🔄 In progress |

---

## Key dependencies

- [`microclimf`](https://github.com/ilyamaclean/microclimf) — mechanistic microclimate model (Maclean 2026)
- [`microclimdata`](https://github.com/ilyamaclean/microclimdata) — automated input data acquisition
- [`mcera5`](https://github.com/dklinges9/mcera5) — ERA5 climate data download
- [`rgee`](https://github.com/r-spatial/rgee) — Google Earth Engine interface from R
- [`plotly`](https://plotly.com/r/) — interactive 3D visualisation
- [`gganimate`](https://gganimate.com/), [`gifski`](https://gif.ski/) — dispersal animations
- [`patchwork`](https://patchwork.data-imaginist.com/) — multi-panel plots
- [`sf`](https://r-spatial.github.io/sf/), `ggplot2`, `ggspatial`, `rnaturalearth` — spatial visualisation
- `terra`, `parallel`, `dplyr`, `readr`

---

## Notes

- **Johansson (1974) zones** (JZ1–JZ5) are used to assign habitat suitability and bark surface area within `build_forest()`. Zone boundaries are proportional canopy height: JZ1 < 10 %, JZ2 < 30 %, JZ3 < 50 %, JZ4 < 80 %, JZ5 = emergent crown. Carrying capacity per voxel is derived from trunk or effective crown surface area divided by mean epiphyte footprint.
- **Forest inventory baseline:** Myster (2017), Maquipucuna primary cloud forest, 1400 m: mean dsh 22.7 cm (trunk radius 0.114 m), 272–324 stems/ha (≥10 cm dsh). Used as default `forestparams` in one-site and all-sites runs.
- **Incremental microclimate saves:** `run_microclimate_site.R` saves each height to its own `.rds` in Lustre scratch before proceeding. If the job is killed (e.g. SLURM timeout or OOM), resubmitting resumes from the last completed height automatically.
- **microenv manifest format:** `microenv_<site>[_h<step>].rds` is a small list with `.heights` (numeric vector), `.height_dir` (path to scratch), and `.weather` (ERA5-cell baseline). The full per-height spatial arrays stay in scratch. Load a specific height with `readRDS(file.path(microenv$.height_dir, sprintf("h%.2f.rds", h)))`.
- **Observations CSV:** the pipeline reads `data/csv/combined_with_identification.csv` by default (all 7 sites; overridable via the `CANOPY_OBS_CSV` environment variable). `combinedv3.csv` — the earlier, more manually curated 5-site file — remains available as an override where curation quality matters more than the extra sites (`paths.R`); species IDs in the default file currently live in the `Identification` column pending manual promotion into `FinalID`.
- **Canopy height ceiling:** the top of the modelled canopy at each site prefers the field-measured `CanopyHeight_m` column in the observations CSV; where that's missing, falls back to the 99th percentile of a GEE canopy-height raster (`vhgt.tif`); and only as a last resort falls back to the tallest recorded epiphyte observation. Using the tallest *observed individual* alone would truncate the canopy below its true height at any site where the tallest recorded epiphyte happened to grow lower than the surrounding forest.
- **Species niche characterization:** `characterize_niches.R` pools each species' climate-niche observations across *every* site it was recorded at (not just the site being simulated), saving `data/processed/species_niches.rds`. Rerun it whenever the observations CSV gets new observations — every colonization run downstream picks up the refined niches automatically.
- **Persistence validation (`best_case.rds`, `realistic.rds`):** `best_case.rds` is a deliberately generous parameter set (every vital rate pushed to its most favourable tested value), used to confirm the model can sustain a population at all before interpreting non-persistence elsewhere as a genuine parameter effect rather than stochastic bad luck. `realistic.rds` runs the same check under literature-default values (no parameter pushed to an extreme), as the complementary lower/baseline bound. Together the two runs are used to decide how the sensitivity factorial's parameter ranges should be bracketed. Both reused as fixed parameter sets for the height-resolution outcome comparison (`height_resolution_experiment.R`).
- **Known issue — animated 3D abundance-over-time plot:** `plot_3d_abundance_animated()` (`get_colonization.R`) renders and saves, but the output currently looks visually wrong (not yet root-caused) and is not confirmed necessary for the thesis. Treat it as **not working** for now — use the static per-replicate abundance PNGs (`plot_abundance()`) and the static 3D snapshot (`plot_3d_abundance()`) as the reliable views instead.
- **Production resolution:** 10 m horizontal / 0.25 m vertical height-tier spacing. Chosen via `resolution_diagnostics.R` (timing: climate-cache build time scales linearly with height-tier count and dominates total cost, so height resolution — not horizontal resolution — governs compute budget at scale) and `height_resolution_experiment.R` (outcomes: full best-case runs at each candidate height step, to confirm the coarser spacing doesn't change results).
- **HPC path:** `paths.R` sets `BASE_DIR` to the HPC home directory. The `CANOPY_PYTHON` environment variable controls which Python interpreter is used for Earth Engine; it defaults to the `canopy_rgee` conda environment.
- **Runtime patches:** two monkey-patches are applied in `patches.R` to fix known bugs in `ecmwfr` and `microclimdata` without modifying package source.
- **Credentials:** `credentials.rds` (CDS API, NASA Earthdata, Google credentials) is excluded from version control. You will need your own.
- **Array dimensions:** the colonization model tracks `[x, y, z, timestep, species]` per life stage. Heights are scoped to the current site's observed range to keep the z-dimension manageable.

---

*Lizeth Estévez Tobar — Universität Bonn, 2026*  
*Supervisor: Juliano Sarmento Cabral*
