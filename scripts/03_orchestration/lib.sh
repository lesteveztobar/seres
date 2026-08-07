#!/bin/bash
# lib.sh — shared shell helpers for the 03_orchestration SLURM launcher
# scripts. Consolidates join_dep(), previously copy-pasted byte-identical
# into run_pipeline.sh, run_factorial_remaining_sites.sh,
# run_full_analysis_pipeline.sh, and run_persistence_all_sites.sh.
#
# Source this after `cd`-ing to the repo root (all four callers already do
# `cd "$(dirname "$0")/../.."` before sourcing), e.g.:
#   source scripts/03_orchestration/lib.sh

# Formats a SLURM --dependency=afterok:<id>:<id>:... string from the given
# job IDs.
join_dep() {
  local out="afterok"
  for id in "$@"; do out="$out:$id"; done
  echo "$out"
}
