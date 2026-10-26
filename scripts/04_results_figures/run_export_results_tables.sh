#!/bin/bash
# run_export_results_tables.sh -- Phase 1 exports (export_results_tables.R).
#   sbatch scripts/04_results_figures/run_export_results_tables.sh
#SBATCH --job-name=export_results_tables
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --output=/home/s38leste_hpc/seres/logs/results_figs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/results_figs/%x_%j.err
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
cd /home/s38leste_hpc/seres
# The realistic runs' own inputs (scripts/03_orchestration/run_full_rerun.sh).
export CANOPY_OBS_CSV=data/csv/combinedv6.csv
export CANOPY_NICHE_CACHE=data/processed/species_niches.rds
export CANOPY_NICHE_BACKGROUND=data/processed/niche_background_density.rds
Rscript scripts/04_results_figures/export_results_tables.R
