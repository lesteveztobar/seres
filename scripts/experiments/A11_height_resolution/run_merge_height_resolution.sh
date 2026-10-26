#!/bin/bash
set -e
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
cd /home/s38leste_hpc/seres
Rscript scripts/experiments/A11_height_resolution/merge_and_analyze_height_resolution.R

