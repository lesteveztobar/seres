#!/bin/bash
# diag_compare_r_cpp_job.sh — Step 2 correctness comparison, R vs Cpp output.
# Not part of the pipeline; submitted manually.
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=900G
#SBATCH --job-name=diag_compare_r_cpp
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/diag_compare_r_cpp_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/HDF5/1.14.3-gompi-2023b/lib:/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"

cd /home/s38leste_hpc/canopymicroenv
Rscript scripts/01_microclimate/diag_compare_r_cpp.R
