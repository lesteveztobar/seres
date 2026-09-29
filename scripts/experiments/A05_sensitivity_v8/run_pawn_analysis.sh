#!/bin/bash
# run_pawn_analysis.sh -- generic launcher for pawn_analysis.R.
# Usage: sbatch run_pawn_analysis.sh <tag> [ranges_csv]
#SBATCH --job-name=pawn_analysis
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
module purge
module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a
cd /home/$USER/seres
export CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combinedv6.csv}"
Rscript scripts/experiments/A05_sensitivity_v8/pawn_analysis.R "$@"
