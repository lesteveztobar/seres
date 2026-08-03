#!/bin/bash
# diag_single_height_job.sh — isolated single-height diagnostic run.
# Not part of the pipeline; submitted manually for one-off resourcing diagnostics.
# 1 cpu, no mclapply contention -- isolates the true per-height cost.
# /usr/bin/time is unavailable on this cluster, so memory/time are read from
# sacct (cgroup-accounted MaxRSS/Elapsed) after the job completes.
#
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=900G
#SBATCH --job-name=diag_single_height
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/diag_single_height_%j.out

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
Rscript scripts/01_microclimate/diag_single_height.R LaElenita 0.25 5.10
