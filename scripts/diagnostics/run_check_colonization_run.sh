#!/bin/bash
# run_check_colonization_run.sh — check_colonization_run.R with the
# module/library setup plot_functions.R needs (it loads sf, which needs
# GDAL/PROJ/GEOS/UDUNITS).
# Usage: sh scripts/diagnostics/run_check_colonization_run.sh [site] [exp_tag]
#   Defaults to Maquipucuna / best_case_h0.25 if no args given.
# Run from: /home/s38leste_hpc/seres/
# ─────────────────────────────────────────────────────────────────────────────

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"

cd /home/$USER/seres
Rscript scripts/diagnostics/check_colonization_run.R "$@"
