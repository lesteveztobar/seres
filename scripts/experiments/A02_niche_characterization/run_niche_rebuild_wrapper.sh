#!/bin/bash
#SBATCH --job-name=niche_rebuild_v8
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1400G
#SBATCH --mail-type=FAIL,TIME_LIMIT,TIME_LIMIT_80,END
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
set -e
module purge
module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a
cd /home/$USER/seres
export CANOPY_OBS_CSV=data/csv/combinedv6.csv
OLD=data/processed/archive_niche_pre_v8_$(date +%Y%m%d)
mkdir -p $OLD; cp data/processed/species_niches.rds data/processed/niche_background_density.rds $OLD/
Rscript scripts/experiments/A02_niche_characterization/characterize_niches.R 0.4
Rscript scripts/experiments/A02_niche_characterization/compare_niche_caches.R $OLD
