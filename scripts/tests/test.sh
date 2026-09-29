#!/bin/bash
#SBATCH --job-name=la_elenita_climtest
#SBATCH --output=/home/s38leste_hpc/seres/logs/la_elenita_spatial_compare.%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/la_elenita_spatial_compare.%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=lm_medium
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s38leste@uni-bonn.de

module purge
module load R

cd /home/s38leste_hpc/seres
Rscript scripts/tests/test.R --compare