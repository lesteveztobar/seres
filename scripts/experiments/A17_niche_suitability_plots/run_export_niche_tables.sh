#!/bin/bash
# run_export_niche_tables.sh — CSV exports of the niche inputs and the
# per-timestep $summary of every Results run (export_niche_tables.R).
# Pinned to the cache the Results runs used (species_niches.rds, v8 rebuild
# of 2026-09-25) and checks n_train against that build's log.
#
# Usage (from /home/s38leste_hpc/seres/):
#   sbatch scripts/experiments/A17_niche_suitability_plots/run_export_niche_tables.sh
#
#SBATCH --job-name=export_niche_tables
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
# ─────────────────────────────────────────────────────────────────────────────
set -e
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"

cd /home/$USER/seres
export CANOPY_OBS_CSV=data/csv/combinedv6.csv
export CANOPY_NICHE_CACHE=data/processed/species_niches.rds
export CANOPY_NICHE_BACKGROUND=data/processed/niche_background_density.rds
export CANOPY_NICHE_BUILD_LOG=logs/niche_rebuild_v8_27876673.err
Rscript scripts/experiments/A17_niche_suitability_plots/export_niche_tables.R
