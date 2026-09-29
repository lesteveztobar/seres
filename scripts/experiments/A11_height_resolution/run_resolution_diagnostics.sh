#!/bin/bash
# run_resolution_diagnostics.sh — time climate-cache build (height resolution)
# and a short colonization run (horizontal resolution) across several
# resolution combinations, to find the coarsest resolution that's still fine
# enough before committing to full experiment runs.
#
# Requires microenv_<site>.rds (0.1m) and microenv_<site>_h<step>.rds for any
# coarser steps to already exist — run height_res_array.sh <site> first for
# the coarser ones.
#
# Usage: sbatch run_resolution_diagnostics.sh <site> [height_steps] [horiz_resolutions]
#   e.g.: sbatch run_resolution_diagnostics.sh Maquipucuna "0.1,0.25,0.5,1.0" "5,10,20,40"
# Run from: /home/s38leste_hpc/seres/
#
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --output=/home/s38leste_hpc/seres/logs/log_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

cd /home/$USER/seres

SITE=${1:?"Usage: sbatch run_resolution_diagnostics.sh <site> [height_steps] [horiz_resolutions]"}
HEIGHT_STEPS=$2
HORIZ_RES=$3

Rscript scripts/experiments/A11_height_resolution/resolution_diagnostics.R "$SITE" "$HEIGHT_STEPS" "$HORIZ_RES"
