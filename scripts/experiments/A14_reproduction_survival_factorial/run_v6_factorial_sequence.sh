#!/bin/bash
# run_v6_factorial_sequence.sh — reproduction_factorial_v3, then _v4, then
# survival_factorial_v5 (see make_params.R), for every site, all three
# against combinedv6.csv -- once Saloya's still-in-flight
# combined_with_identification.csv (old dataset) sweep finishes.
#
# Design (per user, 2026-08-28): re-run v3 under the v6 (improved
# identification) dataset FIRST, alongside v4/v5, so v3's old-dataset vs.
# v6-dataset results are directly comparable using the identical parameter
# grid -- isolating the effect of the identification improvement itself --
# and so v4/v5's own results are internally consistent with that same v3(v6)
# baseline, all needed together to later choose the 3+3 parameters for the
# planned (not yet built) factorial v6 design.
#
# Sequencing:
#   0. Waits for Saloya's remaining old-dataset jobs (see SALOYA_OLD_JOBS
#      below) via --dependency=afterany -- this script can be submitted
#      immediately, it just won't start doing anything until that sweep
#      reaches a terminal state (matches this session's established
#      SLURM-native "chain it, don't poll for it" convention). Does NOT
#      cancel or otherwise touch those jobs.
#   1. reproduction_factorial_v3.rds, one site at a time (same design as
#      run_factorial_remaining_sites.sh -- each site's 625-combo job already
#      saturates a full node, so sites are chained afterok rather than run
#      concurrently), tagged reproduction_factorial_v3_v6.
#   2. Once EVERY site's v3(v6) is done: reproduction_factorial_v4.rds
#      (also 625 combos), same one-site-at-a-time chain, tagged
#      reproduction_factorial_v4_v6.
#   3. Once EVERY site's v4(v6) is done: survival_factorial_v5.rds (216
#      combos), same pattern, tagged survival_factorial_v5_v6.
# No auto-plotting step at the end (unlike run_persistence_all_sites.sh) --
# plot_new_figures.R's factorial call is hardcoded to the _v3_ tag; extend it
# for _v6 tags once real results are in, rather than guessing the right
# call shape now.
#
# Per-site resourcing: default run_colonization.sh SBATCH directives
# (lm_short, 1600G, 8h) work for every site except Saloya, which needed
# vlm_long/3500G/3-day limits this session (tallest canopy of any site --
# 49.7m, 125 height tiers -- see SALOYA_* overrides below, applied the same
# way for all three of its factorial-family jobs here).
#
# Usage: sh scripts/experiments/A14_reproduction_survival_factorial/run_v6_factorial_sequence.sh
# Env overrides: HEIGHT_STEP=0.4
# Run from: /home/s38leste_hpc/seres/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../../.."
source scripts/diagnostics/lib.sh

export CANOPY_OBS_CSV=data/csv/combinedv6.csv
export CANOPY_NICHE_CACHE=data/processed/species_niches_v6.rds
export CANOPY_NICHE_BACKGROUND=data/processed/niche_background_density_v6.rds

mapfile -t SITES < <(awk -F, 'NR>1 && $2!="" {print $2}' "$CANOPY_OBS_CSV" | sort -u)
HEIGHT_STEP=${HEIGHT_STEP:-0.4}

# Saloya's still-running/still-pending OLD-dataset (combined_with_
# identification.csv) jobs -- see this session's transcript. afterany (not
# afterok) so the v6 sequence starts regardless of whether any individual
# one of these succeeds or fails, matching this session's established
# convention for the OAT chain itself.
SALOYA_OLD_JOBS=(27169759 27178947 27169761 27169762 27169763 27169764 27169765 27169766)
START_DEP="--dependency=$(join_dep "${SALOYA_OLD_JOBS[@]}" | sed 's/^afterok/afterany/')"

# Saloya-specific resource override, reused for all three of its factorial
# jobs below (same values used earlier this session for its own sweep).
saloya_extra_args() {
  local site="$1"
  if [ "$site" = "Saloya" ]; then
    echo "--partition=vlm_long --time=3-00:00:00 --mem=3500G"
  fi
}

mkdir -p logs
MANIFEST="logs/v6_factorial_sequence_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=()

submit_stage() {
  local params_file="$1" exp_tag="$2" first_dep_arg="$3"
  local prev_deps=()
  local dep_arg="$first_dep_arg"

  for site in "${SITES[@]}"; do
    local extra
    extra=$(saloya_extra_args "$site")
    # shellcheck disable=SC2086
    local job_id
    job_id=$(sbatch --parsable $dep_arg $extra \
      scripts/02_model/run/run_colonization.sh "$site" "$params_file" "$exp_tag" "$HEIGHT_STEP")
    # >&2, not stdout: this function's stdout is captured via $(...) by the
    # caller (LAST_V3=$(submit_stage ...)) to get the final job id -- tee's
    # own stdout copy would otherwise corrupt that capture with every
    # progress line, not just the last one (caught in dry-run testing,
    # 2026-08-28).
    echo "$site: $exp_tag=$job_id" | tee -a "$MANIFEST" >&2
    ALL_IDS+=("$job_id")
    prev_deps=("$job_id")
    dep_arg="--dependency=$(join_dep "${prev_deps[@]}")"
  done
  # Return the LAST job id of this stage, for the next stage to depend on.
  echo "${prev_deps[0]}"
}

echo "== Stage 1: reproduction_factorial_v3 (v6 dataset), one site at a time =="
echo "   (waits for Saloya's old-dataset sweep to reach a terminal state first)"
LAST_V3=$(submit_stage "$(pwd)/data/params/reproduction_factorial_v3.rds" \
  "reproduction_factorial_v3_v6" "$START_DEP")

echo
echo "== Stage 2: reproduction_factorial_v4 (v6 dataset), one site at a time =="
LAST_V4=$(submit_stage "$(pwd)/data/params/reproduction_factorial_v4.rds" \
  "reproduction_factorial_v4_v6" "--dependency=$(join_dep "$LAST_V3")")

echo
echo "== Stage 3: survival_factorial_v5 (v6 dataset), one site at a time =="
submit_stage "$(pwd)/data/params/survival_factorial_v5.rds" \
  "survival_factorial_v5_v6" "--dependency=$(join_dep "$LAST_V4")" > /dev/null

ALL_IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")
echo
echo "All steps submitted -- job map saved to $MANIFEST"
echo "Nothing runs until Saloya's old-dataset jobs (${SALOYA_OLD_JOBS[*]}) finish."
echo "You can close this terminal now; SLURM will run the chain unattended."
echo "Check progress any time with: squeue -u \$USER"
echo "Or after the fact with:       sacct -j $ALL_IDS_CSV --format=JobID,JobName,State,Elapsed"
