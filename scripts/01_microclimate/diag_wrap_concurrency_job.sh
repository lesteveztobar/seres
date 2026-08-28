#!/bin/bash
# diag_wrap_concurrency_job.sh — Step 3 concurrency test, wrap()/unwrap()+mclapply.
# Not part of the pipeline; submitted manually.
# Usage: sbatch diag_wrap_concurrency_job.sh <site> <mode> <n_cores> <n_heights> [height_start] [height_step] [method]

#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1500G
#SBATCH --job-name=diag_wrap_conc
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/diag_wrap_conc_%j.out

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
SITE="${1:-LaElenita}"
MODE="${2:-parallel}"
NCORES="${3:-4}"
NHEIGHTS="${4:-8}"
HSTART="${5:-5.10}"
HSTEP="${6:-0.10}"
METHOD="${7:-Cpp}"
Rscript scripts/01_microclimate/diag_wrap_concurrency.R "$SITE" "$MODE" "$NCORES" "$NHEIGHTS" "$HSTART" "$HSTEP" "$METHOD"
