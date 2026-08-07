#!/bin/bash
# run_persistence_all_sites.sh — niche re-characterization + best_case/realistic
# persistence runs (see make_params.R) for every site, one site at a time, then
# a single plotting pass over everything once every site is done.
#
# Chains everything with SLURM --dependency=afterok (same approach as
# run_pipeline.sh) so the whole thing runs unattended: submit it and walk
# away, no terminal/session needs to stay open. If any step fails, SLURM
# auto-cancels everything chained after it (DependencyNeverSatisfied)
# instead of quietly running on bad inputs.
#
# Sequencing:
#   1. characterize_niches.R + make_params.R -- once, shared by every site
#      (species_niches.rds isn't site-specific; neither is best_case.rds/
#      realistic.rds), so these two run in parallel with each other.
#   2. For each site in SITES, in order: best_case + realistic runs (in
#      parallel with each other, same as a single-site run would do), but
#      each site waits for the PREVIOUS site's pair to finish first -- so at
#      most one site's worth of heavy (256G/4cpu/8h) jobs runs at a time,
#      rather than all 5 sites x 2 configs = 10 heavy jobs competing for
#      cluster resources simultaneously.
#   3. One final plots job (run_plots.sh/plot_all.R) once the last site's
#      pair finishes -- it plots every site's results (and site maps, niche
#      suitability, etc.) in a single pass, so there's no need to replot
#      after every individual site.
#
# Every site is modeled with just the species actually observed there (no
# species_file arg passed to run_colonization.sh) -- see
# run_colonization.R's species_file/params$species_subset for how to
# override that.
#
# If a step fails partway through, whatever ran before it still has valid
# results -- rerun `sbatch scripts/02_model/plots/run_plots.sh` by hand any
# time to plot whatever's available so far, no need to rerun this whole
# script.
#
# Usage: sh scripts/03_orchestration/run_persistence_all_sites.sh
# Env overrides: HEIGHT_STEP=0.25 (must match an already-generated microenv)
# Run from: /home/s38leste_hpc/canopymicroenv/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/03_orchestration/lib.sh

# Master observations CSV -- mirrors paths.R's OBSERVATIONS_CSV (same
# CANOPY_OBS_CSV override, same default). SITES is derived from it rather
# than hardcoded, so a newly added site is picked up automatically.
CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combined_with_identification.csv}"
export CANOPY_OBS_CSV
mapfile -t SITES < <(awk -F, 'NR>1 && $2!="" {print $2}' "$CANOPY_OBS_CSV" | sort -u)
HEIGHT_STEP=${HEIGHT_STEP:-0.25}

mkdir -p logs
MANIFEST="logs/persistence_all_sites_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=()

echo "== Step 1: niche characterization + params regen (once, shared by every site) =="
NICHE_ID=$(sbatch --parsable scripts/02_model/setup/characterize_niches.sh "$HEIGHT_STEP")
PARAMS_ID=$(sbatch --parsable --partition=intelsr_short --account=ag_biob_scabral \
  --time=00:10:00 --ntasks=1 --output=logs/log_%j.out \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; Rscript scripts/02_model/setup/make_params.R")
echo "niche=$NICHE_ID  params=$PARAMS_ID" | tee -a "$MANIFEST"
ALL_IDS+=("$NICHE_ID" "$PARAMS_ID")

PREV_DEPS=("$NICHE_ID" "$PARAMS_ID")

echo
echo "== Step 2: best_case + realistic per site, one site at a time =="
for site in "${SITES[@]}"; do
  DEP_ARG="--dependency=$(join_dep "${PREV_DEPS[@]}")"

  BEST_ID=$(sbatch --parsable "$DEP_ARG" \
    scripts/02_model/run/run_colonization.sh "$site" data/params/best_case.rds best_case "$HEIGHT_STEP")
  REAL_ID=$(sbatch --parsable "$DEP_ARG" \
    scripts/02_model/run/run_colonization.sh "$site" data/params/realistic.rds realistic "$HEIGHT_STEP")

  echo "$site: best_case=$BEST_ID  realistic=$REAL_ID" | tee -a "$MANIFEST"
  ALL_IDS+=("$BEST_ID" "$REAL_ID")
  PREV_DEPS=("$BEST_ID" "$REAL_ID")
done

echo
echo "== Step 3: plot everything once every site is done =="
PLOT_ID=$(sbatch --parsable --dependency="$(join_dep "${PREV_DEPS[@]}")" scripts/02_model/plots/run_plots.sh)
echo "plots=$PLOT_ID" | tee -a "$MANIFEST"
ALL_IDS+=("$PLOT_ID")

ALL_IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")
echo
echo "All steps submitted -- job map saved to $MANIFEST"
echo "You can close this terminal now; SLURM will run the chain unattended."
echo "Check progress any time with: squeue -u \$USER"
echo "Or after the fact with:       sacct -j $ALL_IDS_CSV --format=JobID,JobName,State,Elapsed"
