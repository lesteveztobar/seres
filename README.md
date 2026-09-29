# SERES — spatially explicit recruitment, establishment and survival

> A three-dimensional, individual-based population model of epiphyte
> colonisation driven by downscaled within-canopy microclimate.

---

## What it does

SERES represents a forest canopy as a three-dimensional grid of voxels,
each carrying its own microenvironmental state vector — temperature,
relative humidity and shortwave radiation — downscaled from reanalysis
climate data to the site's terrain and canopy structure. Population state
is tracked in three life stages (seedling, juvenile, adult) per taxon unit
per voxel. Each simulated year runs three sequential passes: **dispersal
and reproduction** (a wind-mediated seed kernel, gated by a fecundity
function chaining flowering, fruiting, pollination, mycorrhizal
germination and first-year survival probabilities), **establishment** (a
per-voxel binomial draw combining local climate suitability with a
taxon-specific climate niche score), and **survival, growth and stage
transitions** (monthly-compounded, climate- and size-dependent). Forest
structure — tree placement, Johansson-zone canopy geometry, and
bark-area-derived carrying capacity — is generated stochastically once per
site and held static within a run. The model is described in full in the
thesis (`report/methods.tex`, Section 2.2, "SERES: model description" --
part of the merged Section 2, Material and Methods).

## Status

Research code supporting a Master's thesis (University of Bonn, Plant
Sciences, 2026). Not a packaged release: expect active development,
in-flight parameter revisions, and scripts still being run interactively
rather than through a single pinned pipeline.

## Inputs

- **ERA5 reanalysis** (~0.25°, ~28 km resolution) — coarse climate forcing.
- **A global canopy height product** (Lang et al., Google Earth Engine
  raster) — fallback canopy-height ceiling where field measurement is
  unavailable. Tested and rejected as the basis for per-site mean tree
  height (see Discussion/Limitations in the thesis); used only for the
  fallback ceiling above.
- **Plot-based stand structure** — forest-inventory measurements from one
  site (Maquipucuna): stem density and crown radius are applied uniformly
  across sites, while mean tree height is generated per site from a
  ceiling-scaled rule (mean = 0.442 × that site's own canopy ceiling)
  whose coefficients are calibrated once against the Maquipucuna plot data
  (Section 2.2.6 of the thesis).
- **Field observation records** — per-individual species identity,
  coordinates and height above ground, expected as a CSV under
  `data/csv/` (see `scripts/02_model/config/paths.R` for the exact path
  and column names currently read).

## Repository layout

```
data/
  csv/            field observation CSVs
  processed/      microclimate manifests, niche caches, colonization run
                   output (.rds)
  params/         parameter-set .rds files consumed by run_colonization.R
docs/             restructuring and analysis notes (not part of the thesis text)
output/           figures and result CSVs referenced by report/
report/           the thesis itself (LaTeX, subfiles per section)
scripts/
  data_prep/      field-data parsing into the combined observation CSV
                   (formerly 00_data_conversion/)
  geojson-csv-sql-conversion-tools/   external GeoJSON-to-CSV helper,
                   called by data_prep/rebuild_combined_csv.py
  01_microclimate/  ERA5/microclimf downscaling per site (core pipeline,
                     unmoved by the reorg below)
  02_model/
    config/       shared path and helper definitions (core, unmoved)
    engine/       get_colonization.R -- the model itself (core, unmoved)
    run/          single-site and batch run entry points (core, unmoved)
    plots/        core figure-generation scripts, plot_all.R/plot_functions.R
                   (unmoved); per-experiment plotting lives under
                   experiments/A17_niche_suitability_plots/ instead
    analysis/, setup/, diagnostics/   small residue still pending manual
                   placement (see docs/orchestrator/2026-09-29/PHASE2_PLAN.csv,
                   category UNCLEAR) -- everything else that used to live
                   here has moved to diagnostics/ or experiments/ below
  diagnostics/    ad hoc checks and one-off investigations, not part of the
                   main pipeline (consolidated from what used to be spread
                   across 01_microclimate/, 02_model/{analysis,diagnostics,
                   setup}, and 03_orchestration/); includes the archived
                   archive_carcap_investigation_2026-09/ unit and
                   make_params.R (parameter-table generation)
  experiments/    one directory per named evaluation experiment (A02-A19;
                   see "Experiments" below), each holding that experiment's
                   driving script(s) and its own run_*.sh submission wrapper
  03_orchestration/  pipeline-level SLURM launchers:
                       run_overnight_20260929.sh, the current full re-run
                       launcher written against the new layout, plus two
                       earlier launchers still pending manual placement
                       (run_downstream_launcher.sh, run_phase_f_approved.sh)
  tests/          smoke tests and format checks
```

The repository was reorganized on 2026-09-29 (scripts regrouped from a
numbered-stage layout into `data_prep/`, `diagnostics/`, and
per-experiment `experiments/A**` directories); see
`docs/orchestrator/2026-09-29/SESSION_REPORT.md` and
`docs/orchestrator/2026-09-29/PHASE2_PLAN.csv` for the full rationale and
the handful of files still flagged for manual review.

`report/`, `docs/`, `output/`, and the contents of `data/` are excluded
from version control (see `.gitignore`), so the paths this README cites
under them exist in the author's working copy, not in a fresh clone.

## How to run

Requires R (developed under 4.5.2; the cluster job scripts load 4.4.2) with
`ggplot2`, `patchwork`, and the packages `microclimf` and `MASS` depend
on; scripts are run via `Rscript` or submitted to a SLURM cluster.

- **Single site, single configuration**: `Rscript
  scripts/02_model/run/run_colonization.R <site> [params_file]`.
- **All sites, one experiment, on a cluster**: `scripts/02_model/run/
  batch_exp.sh` (SLURM array over sites) or `scripts/02_model/run/
  run_colonization.sh` (single-job submission wrapper).
- **The full re-run** (carrying-capacity calibration gate, then
  every experiment below, the sensitivity design, results manifest, plots
  and thesis PDF, chained by SLURM dependencies): `bash
  scripts/03_orchestration/run_overnight_20260929.sh` from the repository
  root. It only submits jobs, so it is safe to start from a login node.
- **Outputs** land in `data/processed/` (`.rds` run output,
  `microenv_<site>*.rds` climate manifests, `species_niches*.rds` niche
  caches) and `output/` (figures, result CSVs).

## Experiments

Evaluation experiments carry names, not letter identifiers, and are
described in full in `report/methods.tex` (Section 2.3, "Design of the
evaluation experiments" -- merged from the former evaluation.tex, itself
now part of the merged Section 2, Material and Methods) and reported in
`report/results.tex` (Section 3) under the same names, in the same order:

| Experiment | Question | Driving script(s) |
|---|---|---|
| Vertical climate variation | Does the vertical microclimatic gradient retain its form across the elevational gradient? | `scripts/experiments/A08_climate_variation_between_sites/climate_variation_test.R`, `scripts/experiments/A09_elevation_climate_exchange/climate_variation_relative_height.R` |
| Climate covariation | Do radiation, temperature and humidity covary through the profile? | `scripts/experiments/A07_climate_decoupling/climate_decoupling.R` |
| Held-out validation | Is realised vertical position better explained by modelled microenvironment than by height or site climate, and is there an elevation–height exchange rate? | `scripts/experiments/A09b_held_out_validation/held_out_validation.R`, `scripts/02_model/analysis/elevation_canopy_exchange.R` (UNCLEAR-verdict, still pending manual placement) |
| Persistence, founder number and the reproduction–survival factorial | Which demographic constraint most limits establishment? | `scripts/02_model/run/run_colonization.R` (persistence/founder-number/factorial configurations), `scripts/diagnostics/make_params.R`; submission launchers in `scripts/experiments/A13_persistence/` and `scripts/experiments/A14_reproduction_survival_factorial/` |
| Competition and isolation | Does joint multi-taxon simulation produce competition or niche partitioning? | `scripts/experiments/A10_competition_isolation/make_isolation_species_files.R`, `scripts/experiments/A10_competition_isolation/competition_analysis.R` |
| Height resolution | Are results invariant to vertical discretisation? | `scripts/experiments/A11_height_resolution/resolution_diagnostics.R`, `height_resolution_experiment.R` |
| Horizontal resolution | Does horizontal microclimatic variation change where individuals establish? Run in both pooled and voxel modes (`CANOPY_CLIM_MODE`). | `get_colonization.R` (`CANOPY_CLIM_MODE=pooled\|voxel`), `scripts/experiments/A12_horizontal_resolution/run_horizontal_resolution.sh`, `horizontal_resolution_analysis.R` (paired pooled-vs-voxel test) |

Run matrices for Height resolution, Competition and isolation, and
Horizontal resolution had not been submitted at production scale as of
this writing (see `docs/methods_update_report.md`, "Phase F"); the table
above names the scripts that will produce their output, not necessarily a
completed result.

## Parameters

Full parameter tables are given in `report/appendix.tex`; the values
themselves are built by `scripts/diagnostics/make_params.R`. Broadly:

- **Literature-calibrated**: stage-transition and survival intercepts,
  growth rate, flowering/fruiting functions — sourced from epiphytic
  orchid demographic studies, mostly **outside Maxillariinae** (see
  `report/methods.tex`, Section 2.2.1, "Development goal and design
  principles", for this as a stated scope condition rather than a
  discovered limitation).
- **Assumed**: the stochastic noise magnitude and the mycorrhizal
  germination probability are declared modelling assumptions, not
  empirical estimates, and are swept across at least an order of
  magnitude in the sensitivity designs for that reason. The sensitivity
  design itself lives in `scripts/experiments/A05_sensitivity_v8/`: a
  Latin hypercube over the ranges in `data/params/sensitivity_ranges_<tag>.csv`
  (`build_lhs_design.R`, needs the `lhs` package), run through
  `sensitivity_run.R`, and analysed with PAWN indices, standardised
  regression coefficients and a regime classification tree
  (`pawn_analysis.R`, needs `rpart`).
- **Site-derived**: forest stem density and structural parameters (one
  site's plot data, applied uniformly); per-site canopy height and
  microclimate (site-specific).
- **Calibrated**: carrying capacity. `scripts/diagnostics/calibrate_k.R`
  scales the occupiable bark fraction so that Maquipucuna's carrying
  capacity matches a literature stand density (2,800 orchid stands/ha,
  Alzate-Q et al. 2019, scaled to Ecuadorian Maxillariinae richness: 562.9
  per ha) and writes `data/params/k_calibration.rds`. When that file
  exists, `default_forestparams()` reads it and overrides the two
  carrying-capacity terms for every run; set `CANOPY_K_CALIBRATION` to
  point at a different file, or to an empty string to ignore it.
- **Founder number** is capped at 500 in every parameter set built by
  `make_params.R`; the former 273-founder realistic variant is replaced by
  `realistic_75founders.rds`.

## Citation

> Estévez Tobar, L. (2026). *What sets an epiphyte's place in the canopy?
> Building and evaluating SERES, a three-dimensional model of
> recruitment, establishment and survival in Chocó Andino Maxillariinae.*
> Master's thesis, University of Bonn.

Repository: <https://github.com/lesteveztobar/seres>