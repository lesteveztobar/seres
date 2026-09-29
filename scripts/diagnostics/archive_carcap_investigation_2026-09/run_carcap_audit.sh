#!/bin/bash
#SBATCH --job-name=carcap_audit
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --mail-type=FAIL,TIME_LIMIT,TIME_LIMIT_80
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
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
Rscript scripts/diagnostics/archive_carcap_investigation_2026-09/carcap_audit.R
