#!/bin/bash
# run.sh — single entry point for top-level pipeline orchestration.
# None of the three launchers below carry their own #SBATCH header (they run
# directly on the login node and submit the real compute jobs themselves via
# `sbatch`), so unlike scripts/01_microclimate/run_microclimate.sh, they can be dispatched
# here without any SLURM-header constraint. Each remains a separate,
# independently runnable file underneath.
#
# Run from: /home/s38leste_hpc/seres/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."

usage() {
  cat <<'EOF'
Usage: scripts/diagnostics/run.sh <subcommand> [args...]

Subcommands:
  factorial-remaining
      run_factorial_remaining_sites.sh -- catch-up run of the reproduction
      factorial for whichever sites are still missing it.

  persistence
      run_persistence_all_sites.sh -- best_case.rds / realistic.rds
      persistence validation for every site.

  cleanup-scratch [--yes]
      cleanup_resolution_scratch.sh -- reclaim Lustre scratch from finished
      non-production height-resolution runs. Dry-run unless --yes is given.

  help
      Show this message.
EOF
}

[ $# -ge 1 ] || { usage; exit 1; }
CMD="$1"; shift

case "$CMD" in
  factorial-remaining)  sh "scripts/experiments/A14_reproduction_survival_factorial/run_factorial_remaining_sites.sh" "$@" ;;
  persistence)          sh "scripts/experiments/A13_persistence/run_persistence_all_sites.sh" "$@" ;;
  cleanup-scratch)      sh "scripts/diagnostics/cleanup_resolution_scratch.sh" "$@" ;;
  help|-h|--help) usage ;;
  *) echo "Unknown subcommand: $CMD" >&2; usage; exit 1 ;;
esac
