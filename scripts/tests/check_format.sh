#!/bin/bash
#SBATCH --job-name=check_format
#SBATCH --output=logs/check_format_%j.log
#SBATCH --error=logs/check_format_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=lm_short

mkdir -p logs

Rscript -e '
x <- readRDS("/lustre/scratch/data/s38leste_hpc-seres/microenv_LaElenita_heights/h0.10.rds")
cat("Top-level names:\n")
print(names(x))
cat("\nStructure (names only, 2 levels deep):\n")
str(x, max.level = 2, list.len = 10)
'