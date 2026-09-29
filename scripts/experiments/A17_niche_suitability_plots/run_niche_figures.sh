#!/bin/bash
# run_niche_figures.sh — regenerate ONLY the combined all-sites
# niche-suitability figure family (plot_niche_all_sites.R): Fig A
# (niche_suitability_by_site.png), Fig B (niche_forest_all_sites.png), and
# the three supplements. Separate from run_plots.sh so this ~1 h
# climate-cache-bound job can be resubmitted on its own.
#
# Usage:
#   # default dataset (combined_with_identification.csv + species_niches.rds):
#   sbatch scripts/experiments/A17_niche_suitability_plots/run_niche_figures.sh
#
#   # v6 dataset (combinedv6.csv + the _v6 niche caches):
#   sbatch --export=ALL,\
# CANOPY_OBS_CSV=data/csv/combinedv6.csv,\
# CANOPY_NICHE_CACHE=data/processed/species_niches_v6.rds,\
# CANOPY_NICHE_BACKGROUND=data/processed/niche_background_density_v6.rds \
#     scripts/experiments/A17_niche_suitability_plots/run_niche_figures.sh
#
# Run from: /home/s38leste_hpc/seres/
#
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --output=/home/s38leste_hpc/seres/logs/log_%j.out
# ─────────────────────────────────────────────────────────────────────────────

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"

cd /home/$USER/seres
echo "CANOPY_OBS_CSV=${CANOPY_OBS_CSV:-<paths.R default>}"
echo "CANOPY_NICHE_CACHE=${CANOPY_NICHE_CACHE:-<paths.R default>}"
Rscript scripts/experiments/A17_niche_suitability_plots/plot_niche_all_sites.R
