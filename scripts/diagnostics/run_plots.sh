#!/bin/bash
# run_plots.sh — generate every project plot from whatever results exist in
# data/processed/ and geojson_to_csv/. Safe to re-run any time.
# Usage: sh scripts/diagnostics/run_plots.sh          (local/interactive node)
#        sbatch scripts/diagnostics/run_plots.sh      (SLURM)
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
Rscript scripts/02_model/plots/plot_all.R
