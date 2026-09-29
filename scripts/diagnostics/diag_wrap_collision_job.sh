#!/bin/bash
# diag_wrap_collision_job.sh — isolated wrap()/unwrap()+mclapply fork diagnostic.
# Not part of the pipeline; submitted manually.
# Usage: sbatch diag_wrap_collision_job.sh <site> <mode:serial|parallel> <n_cores> <n_heights> [height_start] [height_step]

#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=900G
#SBATCH --job-name=diag_wrap
#SBATCH --output=/home/s38leste_hpc/seres/logs/diag_wrap_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
export CANOPY_PYTHON="/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python"

cd /home/s38leste_hpc/seres
SITE="${1:-LaElenita}"
MODE="${2:-serial}"
NCORES="${3:-1}"
NHEIGHTS="${4:-1}"
HSTART="${5:-5.10}"
HSTEP="${6:-0.10}"
Rscript scripts/diagnostics/diag_wrap_collision.R "$SITE" "$MODE" "$NCORES" "$NHEIGHTS" "$HSTART" "$HSTEP"
