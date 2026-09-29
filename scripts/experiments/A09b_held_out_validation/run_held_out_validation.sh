#!/bin/bash
# run_held_out_validation.sh -- held-out validation (A09b) + per-site-median
# crossing height, post-processing only. Usage: sbatch run_held_out_validation.sh [run_tag]
#SBATCH --job-name=held_out_validation
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
module purge
module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a Miniforge3/24.1.2-0 UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
cd /home/$USER/seres
export CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combinedv6.csv}"
export CANOPY_NICHE_CACHE="${CANOPY_NICHE_CACHE:-data/processed/species_niches.rds}"
export CANOPY_NICHE_BACKGROUND="${CANOPY_NICHE_BACKGROUND:-data/processed/niche_background_density.rds}"
RUN_TAG="${1:-realistic}"
status=0
Rscript scripts/experiments/A09b_held_out_validation/held_out_validation.R "$RUN_TAG" || status=1
Rscript scripts/02_model/analysis/elevation_canopy_exchange.R || status=1
exit $status
