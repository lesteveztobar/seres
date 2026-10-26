#!/bin/bash
# submit_supervisor_updates.sh -- launch fast post-processing (v10 PAWN + Figure 6)
# and parallel constant-physical-dispersal (maxDispZ_m = 2.0 m) height-resolution runs.
set -euo pipefail
cd /home/s38leste_hpc/seres

OBS=data/csv/combinedv6.csv
NICHE=data/processed/species_niches.rds
BG=data/processed/niche_background_density.rds
ENV="--export=ALL,CANOPY_CLIM_MODE=voxel,CANOPY_OBS_CSV=$OBS,CANOPY_NICHE_CACHE=$NICHE,CANOPY_NICHE_BACKGROUND=$BG"
ACCT="--account=ag_biob_scabral"
PARAMS=data/params/realistic_75founders.rds
A11=scripts/experiments/A11_height_resolution/run_height_resolution_experiment.sh

mkdir -p logs/results_figs

# 1. Fast job: v10 PAWN analysis + Figure 6 regeneration
R_ENV='module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a UDUNITS/2.2.28-GCCcore-13.2.0; export LD_LIBRARY_PATH=/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH; unset PYTHONPATH; export CANOPY_OBS_CSV=data/csv/combinedv6.csv; cd /home/s38leste_hpc/seres'

jid_fast=$(sbatch --parsable $ACCT --job-name=v10_pawn_fig6 --partition=intelsr_short --time=00:45:00 --ntasks=1 --cpus-per-task=2 --mem=32G \
  --output=logs/v10_pawn_fig6_%j.out --error=logs/v10_pawn_fig6_%j.err \
  --wrap="$R_ENV; Rscript scripts/experiments/A05_sensitivity/pawn_analysis_v10.R; Rscript scripts/04_results_figures/fig_niche_taxon_curves.R")
echo "v10_pawn_fig6: $jid_fast"

# 2. Parallel height-resolution jobs for Maquipucuna (0.1, 0.25, 0.5, 0.75, 1.0 m; 0.40 m is unchanged)
for step in 0.1 0.25 0.5 0.75 1.0; do
  tag=$(printf "_step%.2f" "$step")
  jid=$(sbatch --parsable $ACCT $ENV --job-name=hres_Maq_${step} --partition=lm_short --time=02:00:00 --cpus-per-task=8 --mem=400G \
    --output=logs/hres_Maq_${step}_%j.out $A11 Maquipucuna "$step" "$PARAMS" "$tag")
  echo "hres_Maq_${step}: $jid"
done

# 3. Parallel height-resolution jobs for 1.0 m at MindoMirador, MindoTarabita, Saloya, Yanayacu (0.40 m is unchanged)
for site in MindoMirador MindoTarabita Saloya Yanayacu; do
  mem=400G
  if [ "$site" = "Saloya" ] || [ "$site" = "MindoMirador" ]; then mem=800G; fi
  jid=$(sbatch --parsable $ACCT $ENV --job-name=hres_${site}_1.0 --partition=lm_short --time=02:00:00 --cpus-per-task=8 --mem=$mem \
    --output=logs/hres_${site}_1.0_%j.out $A11 "$site" "1.0" "$PARAMS" "_step1.00")
  echo "hres_${site}_1.0: $jid"
done
