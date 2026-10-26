#!/bin/bash
# run_height_res_remaining.sh -- finish A11 (height resolution, 0.4
# vs 1.0 m) at the five sites other than Maquipucuna. The full re-run (run_full_rerun.sh)
# chain never ran them: their h1.00 microclimates failed or were flagged
# FAILED, so --kill-on-invalid-dep cancelled the A11 jobs.
#   MindoMirador, Yanayacu, Saloya: h1.00 manifests were saved (the jobs
#     only errored after "Manifest saved") -> verify, then A11.
#   Mashpi: fresh DTM download fell ~0.5 km short of the domain (fixed in
#     lib.R get_dtm(), 0.01deg pad) -> rebuild h1.00, then A11.
#   MindoTarabita: OOM at 32 cores (11/30 heights) -> resume at 8 cores
#     (existing height files are skipped), then A11.
# Usage (repo root, login node -- only submits):
#   bash scripts/03_orchestration/run_height_res_remaining.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

OBS=data/csv/combinedv6.csv
NICHE=data/processed/species_niches.rds
BG=data/processed/niche_background_density.rds
ENV="--export=ALL,CANOPY_CLIM_MODE=voxel,CANOPY_OBS_CSV=$OBS,CANOPY_NICHE_CACHE=$NICHE,CANOPY_NICHE_BACKGROUND=$BG"
ACCT="--account=ag_biob_scabral"
R_PRE='module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a Miniforge3/24.1.2-0 UDUNITS/2.2.28-GCCcore-13.2.0; export LD_LIBRARY_PATH=/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH; unset PYTHONPATH; cd /home/s38leste_hpc/seres'
MICRO_PRE='export CANOPY_PYTHON=/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python; export LD_PRELOAD=/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3; export CANOPY_SCRATCH=$(ws_allocate seres 90); mkdir -p $CANOPY_SCRATCH/tmp; export TMPDIR=$CANOPY_SCRATCH/tmp'
A11=scripts/experiments/A11_height_resolution/run_height_resolution_experiment.sh
PARAMS=data/params/realistic_75founders.rds

mkdir -p logs
MANIFEST="logs/height_res_remaining_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
log_job() { echo "$1: $2" | tee -a "$MANIFEST"; }
sub() { sbatch --parsable $ACCT "$@"; }
mem_multi() { case "$1" in Saloya|MindoMirador) echo 1600G ;; *) echo 600G ;; esac; }

# Manifest check: fails (so A11 never starts) unless the h1.00 manifest
# exists and lists at least one valid height tier.
check_cmd() {
  echo "$R_PRE; Rscript -e 'source(\"scripts/02_model/config/paths.R\"); source(\"scripts/02_model/config/shared_helpers.R\"); source(\"scripts/02_model/engine/get_colonization.R\"); f <- file.path(PROCESSED_DIR, \"microenv_$1_h1.00.rds\"); h <- microenv_heights(readRDS(f)); cat(\"$1 h1.00:\", length(h), \"tiers, max\", max(h), \"m\\n\"); stopifnot(length(h) > 1)'"
}

for site in MindoMirador Yanayacu Saloya; do
  cid=$(sub $ENV --job-name=check_h1_$site --partition=intelsr_short --time=00:20:00 --mem=16G \
    --output=logs/%x_%j.out --wrap="$(check_cmd $site)")
  log_job "check_h1_$site" "$cid"
  id=$(sub $ENV --dependency=afterok:$cid --kill-on-invalid-dep=yes --job-name=A11_hres_$site \
    --time=16:00:00 --mem=$(mem_multi $site) $A11 "$site" "0.4,1.0" $PARAMS)
  log_job "A11_hres_$site" "$id"
done

for site in Mashpi MindoTarabita; do
  cpus=32; [ "$site" = MindoTarabita ] && cpus=8
  mid=$(sub --job-name=microenv_${site}_h1.00 --partition=lm_medium --time=24:00:00 --cpus-per-task=$cpus --mem=1500G \
    --output=logs/%x_%j.out --export=ALL,CANOPY_OBS_CSV=$OBS \
    --wrap="$R_PRE; $MICRO_PRE; Rscript scripts/01_microclimate/run_microclimate_site.R $site 12 1.0")
  log_job "microenv_${site}_h1.00" "$mid"
  # afterany + manifest check: the microclimate script has errored *after*
  # saving a complete manifest before, which afterok would treat as failure.
  cid=$(sub $ENV --dependency=afterany:$mid --kill-on-invalid-dep=yes --job-name=check_h1_$site \
    --partition=intelsr_short --time=00:20:00 --mem=16G --output=logs/%x_%j.out --wrap="$(check_cmd $site)")
  log_job "check_h1_$site" "$cid"
  id=$(sub $ENV --dependency=afterok:$cid --kill-on-invalid-dep=yes --job-name=A11_hres_$site \
    --time=16:00:00 --mem=$(mem_multi $site) $A11 "$site" "0.4,1.0" $PARAMS)
  log_job "A11_hres_$site" "$id"
done
echo "Manifest: $MANIFEST"
