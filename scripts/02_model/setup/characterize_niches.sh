#!/bin/bash
# characterize_niches.sh — precompute each species' realized climate niche
# pooled across every site it was observed at, saving
# data/processed/species_niches.rds for init_colonization() to use.
#
# Rerun this whenever you add observations to the master CSV (paths.R's
# OBSERVATIONS_CSV -- data/csv/combined_with_identification.csv by default,
# override by exporting CANOPY_OBS_CSV before submitting) — every
# colonization run downstream picks up the refined niches automatically.
# Requires microenv_<site>[_h<step>].rds for every site with observations
# already generated (run_microenv.sh / microenv_array.sh).
#
# Usage: sbatch scripts/02_model/setup/characterize_niches.sh [height_step]
#   height_step defaults to 0.4 (the production resolution as of
#   2026-08-17, was 0.25).
# Run from: /home/s38leste_hpc/canopymicroenv/
#
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# --mem=900G (was 128G, 2026-07-27): each per-height climate file is ~6GB
# loaded (run_microclimate_site.R's own comment: 10 spatial arrays x
# 374x372x288) -- with build_clim_cache() now parallelized across 32 forked
# workers (get_colonization.R), peak concurrent memory is ~32 x 6GB = 192GB,
# which OOM-killed this job at 128G. lm_medium nodes have 2TB available.
#SBATCH --mem=900G
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

cd /home/$USER/canopymicroenv
Rscript scripts/02_model/setup/characterize_niches.R "${1:-0.4}"
