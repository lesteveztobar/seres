#!/bin/bash
#SBATCH --job-name=du_inventory
#SBATCH --output=logs/du_inventory_%j.log
#SBATCH --error=logs/du_inventory_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=lm_short

mkdir -p logs

for site in LaElenita Maquipucuna Mashpi MindoMirador MindoTarabita Saloya Yanayacu; do
  for dir in /lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_${site}_heights \
             /lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_${site}_h0.25_heights \
             /lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_${site}_h0.50_heights \
             /lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_${site}_h1.00_heights; do
    if [ -d "$dir" ]; then
      echo "$(du -sh "$dir")"
    else
      echo "MISSING: $dir"
    fi
  done
done