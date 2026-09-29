#!/bin/bash
# cleanup_resolution_scratch.sh — reclaim scratch space from non-production
# height-tier data once it's no longer needed.
#
# Production colonization runs (realistic_273founders baseline, isolation
# runs, reproduction factorial -- steps 4/5/8, run_full_analysis_pipeline.sh)
# only ever read the 0.25m microenv. The 0.1/0.5/1.0m variants exist solely
# for the height-resolution sensitivity comparison (step 7,
# height_resolution_experiment.R) -- once that comparison has produced its
# output file for a site, those variants' per-height scratch directories
# (CANOPY_SCRATCH/microenv_<site>_heights [[0.1m, unsuffixed]],
# microenv_<site>_h0.50_heights, microenv_<site>_h1.00_heights) are dead
# weight. The 0.25m directory (microenv_<site>_h0.25_heights) is NEVER
# touched here -- it stays needed indefinitely by ongoing/future production
# runs for that site.
#
# Safety: dry-run by default (just prints what it would delete and the
# space it would reclaim). Pass --yes to actually delete.
#
# Usage:
#   sh scripts/diagnostics/cleanup_resolution_scratch.sh          # dry run
#   sh scripts/diagnostics/cleanup_resolution_scratch.sh --yes     # actually delete
# Run from: /home/s38leste_hpc/seres/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."

DO_DELETE=0
[ "${1:-}" = "--yes" ] && DO_DELETE=1

CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combined_with_identification.csv}"
mapfile -t SITES < <(awk -F, 'NR>1 && $2!="" {print $2}' "$CANOPY_OBS_CSV" | sort -u)

RESOLUTION_PARAMS_TAG="realistic_75founders"
SCRATCH="${CANOPY_SCRATCH:-}"
if [ -z "$SCRATCH" ]; then
  SCRATCH=$(ws_find seres 2>/dev/null || true)
fi
if [ -z "$SCRATCH" ] || [ ! -d "$SCRATCH" ]; then
  echo "ERROR: could not resolve the scratch workspace (ws_find seres failed and CANOPY_SCRATCH is unset)." >&2
  exit 1
fi
echo "Scratch workspace: $SCRATCH"
[ "$DO_DELETE" -eq 0 ] && echo "(DRY RUN -- pass --yes to actually delete)"
echo

TOTAL_BYTES=0
for site in "${SITES[@]}"; do
  resolution_out="data/processed/height_resolution_${site}_${RESOLUTION_PARAMS_TAG}.rds"
  if [ ! -f "$resolution_out" ]; then
    echo "$site: resolution test not done yet ($resolution_out missing) -- skipping."
    continue
  fi
  echo "$site: resolution test done -- checking non-0.25m scratch dirs..."
  for dirname in "microenv_${site}_heights" "microenv_${site}_h0.50_heights" "microenv_${site}_h1.00_heights"; do
    d="$SCRATCH/$dirname"
    [ -d "$d" ] || continue
    sz_bytes=$(du -sb "$d" 2>/dev/null | cut -f1)
    sz_human=$(du -sh "$d" 2>/dev/null | cut -f1)
    TOTAL_BYTES=$((TOTAL_BYTES + sz_bytes))
    if [ "$DO_DELETE" -eq 1 ]; then
      echo "  DELETING $d ($sz_human)"
      rm -rf "$d"
    else
      echo "  would delete: $d ($sz_human)"
    fi
  done
done

echo
TOTAL_HUMAN=$(numfmt --to=iec-i --suffix=B "$TOTAL_BYTES" 2>/dev/null || echo "${TOTAL_BYTES} bytes")
if [ "$DO_DELETE" -eq 1 ]; then
  echo "Reclaimed: $TOTAL_HUMAN"
else
  echo "Would reclaim: $TOTAL_HUMAN -- rerun with --yes to actually delete."
fi
