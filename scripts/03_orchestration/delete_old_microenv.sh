#!/bin/bash
#SBATCH --job-name=del_%x_%j
#SBATCH --output=logs/delete_%x_%j.log
#SBATCH --error=logs/delete_%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=lm_short

set -euo pipefail

SITE="${1:?Usage: sbatch delete_old_microenv.sh <SiteName>}"
BASE=/lustre/scratch/data/s38leste_hpc-canopymicroenv
MANIFEST_DIR=/home/s38leste_hpc/canopymicroenv/data/processed

mkdir -p logs

echo "=== Deleting old-format microenv data for site: $SITE ==="
echo "Started: $(date)"

for dir in "${BASE}/microenv_${SITE}_heights" \
           "${BASE}/microenv_${SITE}_h0.25_heights" \
           "${BASE}/microenv_${SITE}_h0.50_heights" \
           "${BASE}/microenv_${SITE}_h1.00_heights"; do
  if [ -d "$dir" ]; then
    echo "-- Deleting: $dir"
    du -sh "$dir"
    find "$dir" -type f -delete
    rmdir "$dir"
    echo "   Done."
  else
    echo "-- SKIP (not found): $dir"
  fi
done

for manifest in "${MANIFEST_DIR}/microenv_${SITE}.rds" \
                "${MANIFEST_DIR}/microenv_${SITE}_h0.25.rds" \
                "${MANIFEST_DIR}/microenv_${SITE}_h0.50.rds" \
                "${MANIFEST_DIR}/microenv_${SITE}_h1.00.rds"; do
  if [ -f "$manifest" ]; then
    echo "-- Deleting stale manifest: $manifest"
    rm "$manifest"
  else
    echo "-- SKIP (not found): $manifest"
  fi
done

echo "Finished: $(date)"
echo "=== Done with $SITE ==="
