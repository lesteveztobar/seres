#!/bin/bash
# run_microclimate.sh — single entry point for the microclimate stage.
#
# Dispatches to the underlying scripts. microenv_array.sh carries its own
# #SBATCH resource directives (partition/time/mem/cpus), which SLURM requires
# at the top of the exact file passed to `sbatch` -- it stays a separate file
# for that reason, and generates a throwaway per-site job script (with its
# own #SBATCH header) at submission time rather than delegating to another
# persisted template file. `submit-site` below reuses that same mechanism via
# microenv_array.sh's single-site filter, so there's only one place the
# per-site SBATCH directives live. Everything else here (the plain R entry
# points) has no SBATCH header and can be dispatched directly.
#
# Run from: /home/s38leste_hpc/canopymicroenv/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."
HERE="scripts/01_microclimate"

usage() {
  cat <<'EOF'
Usage: scripts/01_microclimate/run_microclimate.sh <subcommand> [args...]

Direct (no SLURM submission -- assumes the run environment, e.g. an
interactive SLURM allocation or a node with modules/conda already loaded):
  site SITE [N_MONTHS] [HEIGHT_STEP]
      Run run_microclimate_site.R directly for one site (default
      N_MONTHS=12, HEIGHT_STEP=0.1). This is what the generated per-site job
      invokes under SLURM -- use `submit-site` below for production runs.

  progress [SITES] [HEIGHT_STEPS]
      Run check_microenv_progress.R (read-only). SITES and HEIGHT_STEPS are
      comma-separated, e.g. `progress Maquipucuna,Mashpi 0.1,0.25,0.5,1.0`.

  fix-dtm
      Run regenerate_missing_dtm.R: re-fetch any site's missing dtm.tif.

SLURM submission (pass-through to the SBATCH-bearing script):
  submit-site SITE [N_MONTHS] [HEIGHT_STEP]
      sbatch microenv_array.sh restricted to one site -- single-site
      production run.

  submit-array [N_MONTHS] [HEIGHT_STEP]
      sbatch microenv_array.sh -- submits one per-site job for every site.

  help
      Show this message.
EOF
}

[ $# -ge 1 ] || { usage; exit 1; }
CMD="$1"; shift

case "$CMD" in
  site)          Rscript "$HERE/run_microclimate_site.R" "$@" ;;
  progress)      Rscript "$HERE/check_microenv_progress.R" "$@" ;;
  fix-dtm)       Rscript "$HERE/regenerate_missing_dtm.R" "$@" ;;
  submit-site)   sbatch "$HERE/microenv_array.sh" "${2:-12}" "${3:-0.1}" "$1" ;;
  submit-array)  sbatch "$HERE/microenv_array.sh" "$@" ;;
  help|-h|--help) usage ;;
  *) echo "Unknown subcommand: $CMD" >&2; usage; exit 1 ;;
esac
