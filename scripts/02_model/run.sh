#!/bin/bash
# run.sh — single entry point for the complex_model SLURM/shell scripts.
#
# Seven of these carry their own #SBATCH resource headers (partition/time/
# mem/cpus tuned per job type) and are dispatched here as pass-throughs to
# `sbatch` rather than merged -- SLURM requires those directives at the top
# of the exact file you submit, so folding them into a shared script would
# lose per-job resource control. The two that don't carry a header
# (progress/check-run) are dispatched directly.
#
# For the plain R driver/diagnostic/plotting scripts, see run.R in this
# directory.
#
# Run from: /home/s38leste_hpc/canopymicroenv/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."
HERE="scripts/02_model"

usage() {
  cat <<'EOF'
Usage: scripts/02_model/run.sh <subcommand> [args...]

SLURM submission (pass-through to sbatch; each script keeps its own
#SBATCH resource header):
  submit-colonization SITE [PARAMS] [EXPTAG] [HEIGHT_STEP] [SPECIES_FILE]
      sbatch run_colonization.sh -- one site x one experiment.

  submit-batch-exp
      sbatch batch_exp.sh -- regenerates params, then submits the full
      sites x experiments matrix (many submit-colonization jobs).

  submit-niches [HEIGHT_STEP]
      sbatch characterize_niches.sh (default height_step=0.4).

  submit-resolution-diagnostics SITE [HEIGHT_STEPS] [HORIZ_RES]
      sbatch run_resolution_diagnostics.sh -- timing only.

  submit-resolution-experiment SITE [HEIGHT_STEPS] [PARAMS_FILE]
      sbatch run_height_resolution_experiment.sh -- full outcome comparison.

  submit-height-res-array SITE
      sbatch height_res_array.sh -- generate coarser height-step microenv
      variants for one site (chains 3 sequential jobs internally).

  submit-plots
      sbatch run_plots.sh.

Direct (no SBATCH header, runs immediately on the login node):
  progress SITE EXPTAG
      check_colonization_progress.sh -- read-only run-progress check.

  check-run [SITE] [EXPTAG]
      run_check_colonization_run.sh -- thin wrapper for check_colonization_run.R.

  help
      Show this message.
EOF
}

[ $# -ge 1 ] || { usage; exit 1; }
CMD="$1"; shift

case "$CMD" in
  submit-colonization)            sbatch "$HERE/run/run_colonization.sh" "$@" ;;
  submit-batch-exp)               sbatch "$HERE/run/batch_exp.sh" "$@" ;;
  submit-niches)                  sbatch "$HERE/setup/characterize_niches.sh" "$@" ;;
  submit-resolution-diagnostics)  sbatch "$HERE/resolution/run_resolution_diagnostics.sh" "$@" ;;
  submit-resolution-experiment)   sbatch "$HERE/resolution/run_height_resolution_experiment.sh" "$@" ;;
  submit-height-res-array)        sbatch "$HERE/resolution/height_res_array.sh" "$@" ;;
  submit-plots)                   sbatch "$HERE/plots/run_plots.sh" "$@" ;;
  progress)                       sh "$HERE/diagnostics/check_colonization_progress.sh" "$@" ;;
  check-run)                      sh "$HERE/diagnostics/run_check_colonization_run.sh" "$@" ;;
  help|-h|--help) usage ;;
  *) echo "Unknown subcommand: $CMD" >&2; usage; exit 1 ;;
esac
