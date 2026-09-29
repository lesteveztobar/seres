#!/bin/bash
# run_best_combo_pipeline.sh — chains: (already-running) targeted search ->
# pick_best_combo.R -> per-site best_combo re-run (5 reps, real
# run_replicated() output with obs_val -- the factorial's own summary-only
# output doesn't carry that) -> held_out_validation.R (Phase E rewrite).
set -euo pipefail
cd "$(dirname "$0")/../../.."
source scripts/diagnostics/lib.sh

SEARCH_JOB_IDS=("$@")   # pass the 6 targeted-search job IDs as args
SITES=("Maquipucuna" "Mashpi" "Yanayacu" "MindoMirador" "MindoTarabita" "Saloya")

mkdir -p logs
MANIFEST="logs/best_combo_pipeline_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=()
log_job() { echo "$1: $2" | tee -a "$MANIFEST"; ALL_IDS+=("$2"); }
join_dep_any() { local out="afterany"; for id in "$@"; do out="$out:$id"; done; echo "$out"; }

PICK_ID=$(sbatch --parsable --dependency="$(join_dep_any "${SEARCH_JOB_IDS[@]}")" \
  --partition=intelsr_short --account=ag_biob_scabral --time=00:15:00 --ntasks=1 \
  --output=logs/log_%j.out --error=logs/log_%j.err \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; export CANOPY_OBS_CSV=data/csv/combinedv6.csv; Rscript scripts/experiments/A15_best_combo_search/pick_best_combo.R")
log_job "pick_best_combo" "$PICK_ID"

RERUN_IDS=()
for site in "${SITES[@]}"; do
  id=$(sbatch --parsable --dependency=afterok:"$PICK_ID" \
    --export=ALL,CANOPY_CLIM_MODE=voxel,CANOPY_OBS_CSV=data/csv/combinedv6.csv \
    scripts/02_model/run/run_colonization.sh "$site" "data/params/best_combo_${site}.rds" best_combo 0.4)
  log_job "best_combo_rerun_${site}" "$id"
  RERUN_IDS+=("$id")
done

VALIDATE_ID=$(sbatch --parsable --dependency="$(join_dep_any "${RERUN_IDS[@]}")" \
  --partition=lm_short --account=ag_biob_scabral --time=02:00:00 --ntasks=1 --cpus-per-task=4 --mem=900G \
  --output=logs/log_%j.out --error=logs/log_%j.err \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a Miniforge3/24.1.2-0 UDUNITS/2.2.28-GCCcore-13.2.0; export LD_LIBRARY_PATH=\"/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:\$LD_LIBRARY_PATH\"; export CANOPY_OBS_CSV=data/csv/combinedv6.csv; Rscript scripts/experiments/A09b_held_out_validation/held_out_validation.R")
log_job "held_out_validation" "$VALIDATE_ID"

SUMMARY_ID=$(sbatch --parsable --dependency="$(join_dep_any "${ALL_IDS[@]}")" \
  --job-name=best_combo_pipeline_summary --mail-type=END --mail-user=s38leste@uni-bonn.de \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  scripts/diagnostics/batch_summary_mailer.sh "best_combo_pipeline" "${ALL_IDS[@]}")
log_job "summary" "$SUMMARY_ID"

echo "Best-combo pipeline chained: $(echo "${ALL_IDS[@]}" | wc -w) jobs -- job map: $MANIFEST"
