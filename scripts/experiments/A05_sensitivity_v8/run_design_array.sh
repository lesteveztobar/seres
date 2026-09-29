#!/bin/bash
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=1200G
#SBATCH --mail-type=FAIL,TIME_LIMIT,TIME_LIMIT_80
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%A_%a.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%A_%a.err
# Generic design-array runner (2026-09-24). Required via --export:
#   CANOPY_SENS_DESIGN, CANOPY_SENS_OUTDIR, CANOPY_ARRAY_N ; optional CANOPY_SENS_TIMESTEPS.
# Submit with e.g.: sbatch --job-name=fact_v8 --array=1-75%12 --export=ALL,CANOPY_SENS_DESIGN=...,CANOPY_SENS_OUTDIR=...,CANOPY_ARRAY_N=75 run_design_array.sh
module purge
module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a Miniforge3/24.1.2-0 UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
cd /home/$USER/seres
export CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combinedv6.csv}"
export CANOPY_CLIM_MODE=voxel
export CANOPY_ARRAY_IDX="${SLURM_ARRAY_TASK_ID}"
: "${CANOPY_SENS_DESIGN:?not set}" "${CANOPY_SENS_OUTDIR:?not set}" "${CANOPY_ARRAY_N:?not set}"
Rscript scripts/experiments/A05_sensitivity_v8/sensitivity_run.R
