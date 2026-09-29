#!/bin/bash
# run_horizontal_resolution.sh — the new paired experiment (2026-09-08):
# same configuration, same seeds-per-replicate, run in `pooled` and `voxel`
# CANOPY_CLIM_MODE at every site, 10 replicates each (raised from 3 per
# the author's own instruction -- this is the headline result for the 3D
# rebuild, worth the extra ~cheap insurance). The only experiment run in
# both modes -- everything else in Phase F is voxel-only.
set -euo pipefail
cd "$(dirname "$0")/../../.."
source scripts/diagnostics/lib.sh

SITES=("Maquipucuna" "Mashpi" "Yanayacu" "MindoMirador" "MindoTarabita" "Saloya")
# SKIP_SITES (space-separated, default empty): drop sites from this run, e.g. SKIP_SITES="Saloya"
if [ -n "${SKIP_SITES:-}" ]; then
  _keep=(); for _s in "${SITES[@]}"; do case " $SKIP_SITES " in *" $_s "*) echo "SKIPPING $_s (SKIP_SITES)";; *) _keep+=("$_s");; esac; done
  SITES=("${_keep[@]}")
fi
# 2026-09-24: first pass at 5 reps (not 10) for a quicker initial read.
# horizontal_resolution.rds (10 reps) is kept on disk as the manual,
# NOT-auto-scheduled follow-up if this first pass's pooled/voxel effect
# turns out close to the noise floor (median 3.0%/p90 5.2% CV in the
# bounded regime, report/evaluation.tex) and worth resolving more precisely.
PARAMS="data/params/realistic_75founders.rds"

mkdir -p logs
MANIFEST="logs/horizontal_resolution_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=()
log_job() { echo "$1: $2" | tee -a "$MANIFEST"; ALL_IDS+=("$2"); }
join_dep_any() { local out="afterany"; for id in "$@"; do out="$out:$id"; done; echo "$out"; }

for site in "${SITES[@]}"; do
  id=$(sbatch --parsable --export=ALL,CANOPY_CLIM_MODE=pooled,CANOPY_OBS_CSV=data/csv/combinedv6.csv \
    scripts/02_model/run/run_colonization.sh "$site" "$PARAMS" horizontal_resolution_pooled 0.4)
  log_job "horizontal_resolution_pooled_${site}" "$id"
  id=$(sbatch --parsable --export=ALL,CANOPY_CLIM_MODE=voxel,CANOPY_OBS_CSV=data/csv/combinedv6.csv \
    scripts/02_model/run/run_colonization.sh "$site" "$PARAMS" horizontal_resolution_voxel 0.4)
  log_job "horizontal_resolution_voxel_${site}" "$id"
done

SUMMARY_ID=$(sbatch --parsable --dependency="$(join_dep_any "${ALL_IDS[@]}")" \
  --job-name=horizontal_resolution_summary --mail-type=END --mail-user=s38leste@uni-bonn.de \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  scripts/diagnostics/batch_summary_mailer.sh "horizontal_resolution" "${ALL_IDS[@]}")
log_job "summary" "$SUMMARY_ID"

ALL_IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")
echo "Horizontal resolution submitted: $(echo "${ALL_IDS[@]}" | wc -w) jobs -- job map: $MANIFEST"
