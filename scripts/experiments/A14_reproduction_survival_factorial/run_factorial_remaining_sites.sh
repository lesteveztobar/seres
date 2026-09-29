#!/bin/bash
# run_factorial_remaining_sites.sh — reproduction_factorial_v3 sweep (p_poll x
# p_germ x p_s1 x n_founders, 625 combinations, see make_params.R) for Mashpi
# and MindoTarabita, one site at a time, then a single plotting pass.
#
# Scope note (2026-07-16): Maquipucuna already has a v3 factorial result
# (data/processed/colonization_Maquipucuna_reproduction_factorial_v3_h0.25.rds,
# run 2026-07-14, ~1h15m). MindoMirador and Yanayacu are deliberately excluded
# here -- their best_case runs are crashing with an unresolved "missing value
# where TRUE/FALSE needed" error (0/5 replicates succeeding), on hold until
# that's investigated; running their (much larger, 625-combination) factorial
# sweep first would just hit the same crash 625x instead of 5x, for no benefit.
# SITES is a hand-picked list, not derived from OBSERVATIONS_CSV -- it does
# NOT include LaElenita/Saloya (added 2026-07-24, see
# scripts/data_prep/rebuild_combined_csv.py); this script's whole point is catching up
# a specific historical gap for specific sites, not "run every site," so
# leave it as an explicit list rather than making it dynamic. All of best-
# case/every-site's results are stale anyway (code changed substantially
# since these were run -- size tracking, niche rewrite, timesteps=50), so
# add the two new sites here explicitly once you're re-running this for
# real rather than relying on this script picking them up automatically.
#
# No niche/params regen step here, unlike run_persistence_all_sites.sh --
# species_niches.rds and reproduction_factorial_v3.rds are already fresh from
# that run (nothing in characterize_niches.R/make_params.R has changed since),
# so re-running them would just reproduce byte-identical output.
#
# Chains everything with SLURM --dependency=afterok so the whole thing runs
# unattended: submit it and walk away. If a site's factorial run fails, SLURM
# auto-cancels everything chained after it instead of quietly running on bad
# inputs -- rerun `sbatch scripts/diagnostics/run_plots.sh` by hand any time
# to plot whatever's available so far.
#
# Sequencing: one site's factorial sweep at a time (each site's 625-combo x
# 4-core job already saturates the node it's on; running two sites' sweeps
# concurrently would just make both slower for no benefit), then one plots
# job once the last site finishes.
#
# Usage: sh scripts/experiments/A14_reproduction_survival_factorial/run_factorial_remaining_sites.sh
# Env overrides: HEIGHT_STEP=0.25 (must match an already-generated microenv)
# Run from: /home/s38leste_hpc/seres/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../../.."
source scripts/diagnostics/lib.sh

SITES=("Mashpi" "MindoTarabita")
HEIGHT_STEP=${HEIGHT_STEP:-0.25}
PARAMS_FILE="data/params/reproduction_factorial_v3.rds"

if [ ! -f "$PARAMS_FILE" ]; then
  echo "Missing $PARAMS_FILE -- run make_params.R first (e.g. sh scripts/experiments/A13_persistence/run_persistence_all_sites.sh," \
       "or just: sbatch --wrap=\"module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; Rscript scripts/diagnostics/make_params.R\")" >&2
  exit 1
fi

mkdir -p logs
MANIFEST="logs/factorial_remaining_sites_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=()

echo "== Reproduction factorial (v3) per site, one site at a time =="
PREV_DEPS=()
for site in "${SITES[@]}"; do
  if [ "${#PREV_DEPS[@]}" -gt 0 ]; then
    DEP_ARG="--dependency=$(join_dep "${PREV_DEPS[@]}")"
    JOB_ID=$(sbatch --parsable "$DEP_ARG" \
      scripts/02_model/run/run_colonization.sh "$site" "$PARAMS_FILE" reproduction_factorial_v3 "$HEIGHT_STEP")
  else
    JOB_ID=$(sbatch --parsable \
      scripts/02_model/run/run_colonization.sh "$site" "$PARAMS_FILE" reproduction_factorial_v3 "$HEIGHT_STEP")
  fi

  echo "$site: factorial=$JOB_ID" | tee -a "$MANIFEST"
  ALL_IDS+=("$JOB_ID")
  PREV_DEPS=("$JOB_ID")
done

echo
echo "== Plot everything once the last site is done =="
PLOT_ID=$(sbatch --parsable --dependency="$(join_dep "${PREV_DEPS[@]}")" scripts/diagnostics/run_plots.sh)
echo "plots=$PLOT_ID" | tee -a "$MANIFEST"
ALL_IDS+=("$PLOT_ID")

ALL_IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")
echo
echo "All steps submitted -- job map saved to $MANIFEST"
echo "You can close this terminal now; SLURM will run the chain unattended."
echo "Check progress any time with: squeue -u \$USER"
echo "Or after the fact with:       sacct -j $ALL_IDS_CSV --format=JobID,JobName,State,Elapsed"
