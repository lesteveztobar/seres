#!/bin/bash
# run_height_resolution_experiment.sh — run the full colonization model (real
# timesteps/spinup, several best_case.rds replicates) at each microenv height
# resolution, to check whether height resolution changes RESULTS (not just
# runtime -- see resolution_diagnostics.R for the timing-only comparison).
#
# Requires microenv_<site>.rds (0.1m) and microenv_<site>_h<step>.rds for any
# coarser steps to already exist — run height_res_array.sh <site> first.
#
# Usage: sbatch run_height_resolution_experiment.sh <site> [height_steps] [params_file]
#   e.g.: sbatch run_height_resolution_experiment.sh Maquipucuna "0.1,0.25,0.5,1.0" data/params/best_case.rds
# Run from: /home/s38leste_hpc/seres/
#
#SBATCH --partition=lm_medium
#SBATCH --account=ag_biob_scabral
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
# --mem=1600G (was 300G -> 900G -> 1600G, 2026-07-27): MindoMirador's
# resolution job still showed replicate-level dropout at 900G (1/5, then
# 0/5 succeeded -- consistent with individual mclapply workers getting
# silently OOM-killed) -- same size-tracking memory growth as
# run_colonization.sh, see that script's own comment.
#SBATCH --mem=1600G
#SBATCH --output=/home/s38leste_hpc/seres/logs/log_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

cd /home/$USER/seres

SITE=${1:?"Usage: sbatch run_height_resolution_experiment.sh <site> [height_steps] [params_file]"}
HEIGHT_STEPS=$2
PARAMS_FILE=$3

Rscript scripts/experiments/A11_height_resolution/height_resolution_experiment.R "$SITE" "$HEIGHT_STEPS" "$PARAMS_FILE"
