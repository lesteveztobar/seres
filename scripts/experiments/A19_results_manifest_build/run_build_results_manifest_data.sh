#!/bin/bash
#SBATCH --job-name=build_results_manifest_data
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=400G
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
cd /home/$USER/seres
export CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combinedv6.csv}"
Rscript scripts/experiments/A19_results_manifest_build/build_results_manifest_data.R
