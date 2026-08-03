#!/bin/bash
# run.sh — single entry point for the data-conversion stage (GeoJSON -> structured CSV).
# Dispatches to the underlying scripts, which remain separate, independently
# runnable files (see each script's own --help / docstring for full options).
# Run from: /home/s38leste_hpc/canopymicroenv/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."
HERE="scripts/00_data_conversion"

usage() {
  cat <<'EOF'
Usage: scripts/00_data_conversion/run.sh <subcommand> [args...]

Subcommands:
  convert --input RAW.csv --output OUT.csv [--source NAME] [--site NAME]
      Parse a raw GeoJSON-derived CSV into the structured EpiphytesDatabase
      format (convert_observations.py).

  photos --input CSV [CSV...] --sites SITE [SITE...] --output OUT.csv
         [--offset N] [--sd-card PATH]
      Build the photo-lookup table from raw field-export CSVs
      (build_photo_lookup.py).

  rebuild
      Run the LaElenita/MindoMirador/Saloya migration end to end and rebuild
      data/csv/combined.csv from every site's Processed*.csv
      (rebuild_combined_csv.py). No arguments.

  help
      Show this message.
EOF
}

[ $# -ge 1 ] || { usage; exit 1; }
CMD="$1"; shift

case "$CMD" in
  convert)  python3 "$HERE/convert_observations.py" "$@" ;;
  photos)   python3 "$HERE/build_photo_lookup.py" "$@" ;;
  rebuild)  python3 "$HERE/rebuild_combined_csv.py" "$@" ;;
  help|-h|--help) usage ;;
  *) echo "Unknown subcommand: $CMD" >&2; usage; exit 1 ;;
esac
