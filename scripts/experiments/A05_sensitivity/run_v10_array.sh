#!/bin/bash
#SBATCH --partition=intelsr_short,lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=480G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%A_%a.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%A_%a.err
# v10 campaign array runner (2026-10-02). Required via --export:
#   CANOPY_SENS_DESIGN, CANOPY_SENS_OUTDIR ; optional CANOPY_SENS_TIMESTEPS, CANOPY_RUN_TIMEOUT_S.
# Array index = design$task (interleaved by arm x site, see build_v10_design.R).
module purge
module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a Miniforge3/24.1.2-0 UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
cd /home/$USER/seres
export CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combinedv6.csv}"
export CANOPY_CLIM_MODE=voxel
export CANOPY_ARRAY_IDX="${SLURM_ARRAY_TASK_ID}"
: "${CANOPY_SENS_DESIGN:?not set}" "${CANOPY_SENS_OUTDIR:?not set}"
Rscript scripts/experiments/A05_sensitivity/sensitivity_run_v10.R
