#!/bin/bash
#SBATCH --job-name=literature_diagnostic
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=1600G
#SBATCH --mail-type=FAIL,TIME_LIMIT,TIME_LIMIT_80
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
# Literature-parameter diagnostic against the RECALIBRATED carrying
# capacity (2026-09-10). 30 founders/species, p_poll=0.3, p_germ=0.001,
# p_s1=0.45, S=1.76e6, 50 timesteps, 5 spin-up generations, 5 replicates,
# voxel mode. One site per argument; run both.
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
cd /home/$USER/seres
export CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combinedv6.csv}"
export CANOPY_CLIM_MODE=voxel
SITE="${1:?usage: sbatch run_literature_diagnostic.sh <SITE>}"
Rscript scripts/02_model/run/run_colonization.R "$SITE" data/params/literature_diagnostic.rds literature_diag 0.4
