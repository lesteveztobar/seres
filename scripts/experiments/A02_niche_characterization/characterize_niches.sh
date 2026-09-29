#!/bin/bash
# characterize_niches.sh — precompute each species' realized climate niche
# pooled across every site it was observed at, saving
# data/processed/species_niches_v7.rds / niche_background_density_v7.rds for
# init_colonization() to use.
#
# Rerun this whenever you add observations to the master CSV, or (2026-09-06,
# Phase C of the v7 3D re-run) whenever the background-weighting logic
# changes -- see characterize_niches.R for the corrected 1/pixel_count
# per-site weighting this run applies.
#
# Usage: sbatch scripts/experiments/A02_niche_characterization/characterize_niches.sh [height_step]
#   height_step defaults to 0.4 (the production resolution).
# Run from: /home/s38leste_hpc/seres/
#
#SBATCH --job-name=characterize_niches_v7c
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# 1400G (2026-09-06): the two-voxel-cache-per-site version of this script
# OOM'd at 400G; this run builds only ONE voxel cache per site (cc_voxel,
# reused for both presence and the now-per-pixel-weighted background), so
# peak memory should be lower than that run, but keeping the known-good
# figure rather than re-testing a lower one under time pressure.
#SBATCH --mem=1400G
#SBATCH --mail-type=FAIL,TIME_LIMIT,TIME_LIMIT_80
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

cd /home/$USER/seres
export CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combinedv6.csv}"
Rscript scripts/experiments/A02_niche_characterization/characterize_niches.R "${1:-0.4}"
