#!/bin/bash
# diag_inspect_runmicro.sh — one-off: dump microclimf::runmicro() source to a log.
# Not part of the pipeline.
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=diag_inspect_runmicro
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/diag_inspect_runmicro_%j.out

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
Rscript -e '
source("scripts/02_model/config/patches.R")
library(rgee)
library(readr)
library(mcera5)
library(microclimf)
library(microclimdata)
library(terra)
library(luna)
cat("=== formals(.runmicronosnow) ===\n")
print(args(microclimf:::.runmicronosnow))
cat("\n=== body(.runmicronosnow) ===\n")
b <- deparse(body(microclimf:::.runmicronosnow))
cat(b, sep="\n")
cat("\n=== body(.runmicronosnow) length (lines) ===", length(b), "\n")
' 2>&1
