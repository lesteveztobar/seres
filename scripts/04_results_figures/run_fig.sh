#!/bin/bash
# run_fig.sh -- run one or more results-figure scripts on a compute node.
#   sbatch scripts/04_results_figures/run_fig.sh fig_competition.R [fig_x.R ...]
#SBATCH --job-name=results_fig
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/home/s38leste_hpc/seres/logs/results_figs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/results_figs/%x_%j.err
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export CANOPY_OBS_CSV=data/csv/combinedv6.csv
cd /home/s38leste_hpc/seres
status=0
for f in "$@"; do
  echo "=== $f ==="
  Rscript "scripts/04_results_figures/$f" || { echo "FAILED: $f"; status=1; }
done
exit $status
