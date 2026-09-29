#!/bin/bash
# organize_output.sh -- sorts top-level files in output/ into
#   output/figures/<topic>/, output/tables/<topic>/, output/runs/
# by filename prefix, then rewrites \includegraphics{../output/<file>} in
# report/*.tex to the new location. Moves only (never deletes); files that
# match no rule stay where they are and are listed. Dry run unless --apply.
# Usage: bash scripts/diagnostics/organize_output.sh [--apply]
set -euo pipefail
cd "$(dirname "$0")/../.."
APPLY=0; [ "${1:-}" = "--apply" ] && APPLY=1
OUT=output

topic_for() {
  case "$1" in
    site_map*|tree_diagram*)                                   echo sites ;;
    held_out*|elevation_canopy_crossing*)                      echo validation ;;
    competition*)                                              echo competition ;;
    height_resolution*|resolution_*|horizontal_resolution*)    echo resolution ;;
    sensitivity*|sens_*|pawn*)                                 echo sensitivity ;;
    niche_*)                                                   echo niche ;;
    climate_*|elevation_*|vpd*|footprint_*|raster_*|vhgt_*)    echo climate ;;
    *_audit*|valid_cell*|grid_origin*|home_cell*|monthly_validity*|pooled_fallback*|horizon_margin*|gridsnap*|carcap_*) echo audits ;;
    factorial_*|persistence*|bounded_*|colonization_*|abundance_*|results_manifest*|founder*) echo model ;;
    *) echo "" ;;
  esac
}

move() { # src dest_dir
  if [ "$APPLY" = 1 ]; then mkdir -p "$2"; if [ -d "$1" ]; then mv -n "$1" "$2/"; else mv -f "$1" "$2/"; fi; fi
  echo "  $1 -> $2/"
}

unmatched=()
for f in "$OUT"/*; do
  name=$(basename "$f")
  case "$name" in archive|figures|tables|runs|provisional_*) continue ;; esac
  # Sensitivity arrays may still be running when this fires -- leave their live
  # run directories and the analysis folder in place (any design tag).
  case "$name" in sensitivity_runs_*|sensitivity) continue ;; esac
  if [ -d "$f" ]; then move "$f" "$OUT/runs"; continue; fi
  topic=$(topic_for "$name")
  if [ -z "$topic" ]; then unmatched+=("$name"); continue; fi
  case "$name" in
    *.png|*.pdf|*.svg|*.html|*.jpg) dest="$OUT/figures/$topic" ;;
    *) dest="$OUT/tables/$topic" ;;
  esac
  move "$f" "$dest"
  if [ "$APPLY" = 1 ] && [ "$dest" != "${dest#$OUT/figures}" ]; then
    rel="${dest#$OUT/}"
    for tex in report/*.tex; do
      grep -q "../output/$name}" "$tex" && sed -i "s#../output/$name}#../output/$rel/$name}#g" "$tex" && echo "    updated $tex"
    done
  fi
done
[ ${#unmatched[@]} -gt 0 ] && { echo "Left in place (no rule):"; printf '  %s\n' "${unmatched[@]}"; }
[ "$APPLY" = 1 ] || echo "(dry run -- rerun with --apply to move)"
