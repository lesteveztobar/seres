#!/bin/bash
# diag_subset_fidelity_job.sh — Stage A verification gate for the
# every-other-day temporal subsample. Not part of the pipeline; submitted
# manually. Same resource profile as diag_single_height_job.sh (1 cpu, no
# mclapply contention, single runmicro() call) since diag_subset_fidelity.R
# does exactly one runmicro() call at one height.
#
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=900G
#SBATCH --job-name=diag_subset_fidelity
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/diag_subset_fidelity_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
export CANOPY_PYTHON="/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python"

cd /home/s38leste_hpc/canopymicroenv
Rscript scripts/01_microclimate/diag_subset_fidelity.R Maquipucuna 2.10 2 \
  /lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_Maquipucuna_h0.25_heights/h2.10.rds
