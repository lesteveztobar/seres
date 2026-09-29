#!/bin/bash
# run_v6_factorial_sequence_concurrent.sh — concurrent variant of
# run_v6_factorial_sequence.sh (2026-08-28): reproduction_factorial_v3,
# then _v4, then survival_factorial_v5 (see make_params.R), against
# combinedv6.csv, for every remaining site (Saloya dropped from production
# 2026-08-28 -- tallest canopy/slowest site, not worth the wall-clock
# budget with under a month left; its old-dataset sweep was cancelled, not
# waited on).
#
# Differs from run_v6_factorial_sequence.sh in exactly one way: within a
# stage, every site is submitted CONCURRENTLY (no --dependency between
# sites) instead of chained one-at-a-time. The sequential version's own
# header explains why it originally chained sites: "each site's 625-combo
# job already saturates a full node, so sites are chained afterok rather
# than run concurrently" -- but that's about one job's mclapply workers
# saturating ONE node, not about running multiple separate SLURM jobs (each
# its own node reservation) at once. sinfo confirms plenty of headroom for
# this (2026-08-28: lm_short/lm_medium share ~22 nodes, ~1300+ idle CPUs,
# 2TB/node -- 6 concurrent 1600G jobs fit across separate nodes with room
# to spare). Per user request: try concurrent, monitor for trouble (OOM
# kills, node-health drains, /tmp collisions -- the specific failure modes
# already seen and fixed once this session), fall back to the sequential
# script if it recurs.
#
# Between-STAGE dependencies are kept exactly as the sequential version
# (v4 waits for every site's v3 to succeed; v5 waits for every site's v4)
# -- that's about factorial-design interpretability (see the sequential
# script's header), not resource contention, so concurrency doesn't change
# it. Only the reason to skip is dropped: no more waiting on Saloya's
# now-cancelled old-dataset jobs.
#
# Usage: sh scripts/experiments/A14_reproduction_survival_factorial/run_v6_factorial_sequence_concurrent.sh
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
# Drop Saloya -- see batch_exp.sh's matching 2026-08-28 exclusion.
SITES=("${SITES[@]/Saloya}")
SITES=($(printf '%s\n' "${SITES[@]}" | grep -v '^$'))
HEIGHT_STEP=${HEIGHT_STEP:-0.4}

mkdir -p logs
MANIFEST="logs/v6_factorial_sequence_concurrent_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"

# Submits every site in $SITES concurrently (same $dep_arg for all of
# them -- no chaining between sites), returns the list of job ids on
# stdout (space-separated) for the caller to build the next stage's
# dependency from via join_dep.
submit_stage_concurrent() {
  local params_file="$1" exp_tag="$2" dep_arg="$3"
  local stage_ids=()
  for site in "${SITES[@]}"; do
    # shellcheck disable=SC2086
    local job_id
    job_id=$(sbatch --parsable $dep_arg \
      scripts/02_model/run/run_colonization.sh "$site" "$params_file" "$exp_tag" "$HEIGHT_STEP")
    # >&2, not stdout: this function's stdout is captured via $(...) to get
    # the stage's job ids -- see submit_stage()'s matching note in
    # run_v6_factorial_sequence.sh.
    echo "$site: $exp_tag=$job_id" | tee -a "$MANIFEST" >&2
    stage_ids+=("$job_id")
  done
  echo "${stage_ids[@]}"
}

echo "== Stage 1: reproduction_factorial_v3 (v6 dataset), all sites concurrently =="
STAGE1_IDS=($(submit_stage_concurrent "$(pwd)/data/params/reproduction_factorial_v3.rds" \
  "reproduction_factorial_v3_v6" ""))

echo
echo "== Stage 2: reproduction_factorial_v4 (v6 dataset), all sites concurrently =="
echo "   (waits for every site's stage 1 job to succeed)"
STAGE2_IDS=($(submit_stage_concurrent "$(pwd)/data/params/reproduction_factorial_v4.rds" \
  "reproduction_factorial_v4_v6" "--dependency=$(join_dep "${STAGE1_IDS[@]}")"))

echo
echo "== Stage 3: survival_factorial_v5 (v6 dataset), all sites concurrently =="
echo "   (waits for every site's stage 2 job to succeed)"
STAGE3_IDS=($(submit_stage_concurrent "$(pwd)/data/params/survival_factorial_v5.rds" \
  "survival_factorial_v5_v6" "--dependency=$(join_dep "${STAGE2_IDS[@]}")"))

ALL_IDS=("${STAGE1_IDS[@]}" "${STAGE2_IDS[@]}" "${STAGE3_IDS[@]}")
ALL_IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")
echo
echo "All steps submitted -- job map saved to $MANIFEST"
echo "Stage 1 starts immediately (no dependency)."
echo "You can close this terminal now; SLURM will run the sequence unattended."
echo "Check progress any time with: squeue -u \$USER"
echo "Or after the fact with:       sacct -j $ALL_IDS_CSV --format=JobID,JobName,State,Elapsed"
