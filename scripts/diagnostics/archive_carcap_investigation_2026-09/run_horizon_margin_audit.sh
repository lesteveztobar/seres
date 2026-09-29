#!/bin/bash
#SBATCH --job-name=horizon_margin_audit
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
cd /home/$USER/seres
Rscript scripts/diagnostics/archive_carcap_investigation_2026-09/horizon_margin_audit.R
