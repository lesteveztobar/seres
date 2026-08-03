#!/bin/bash
# run_colonization.sh — run colonization model for a single site × experiment
# Called by batch_exp.sh; do not submit directly.
#
#SBATCH --partition=lm_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# --mem=1600G (was 256G -> 900G -> 1600G, 2026-07-27): Saloya's
# realistic_273founders run OOM-killed even at 900G (881GB used, right at
# the ceiling) -- explicit per-voxel size_S/size_J/size_A tracking roughly
# triples each replicate's state array footprint versus abundance-only
# tracking, and Saloya/MindoMirador's tall canopies (many height tiers)
# make this worse than at other sites. lm_short and lm_medium share the
# same physical nodes (2TB each, sinfo) -- 1600G leaves ~400G headroom for
# the OS/other jobs sharing the node rather than claiming the whole thing.
#SBATCH --mem=1600G
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"

cd /home/$USER/canopymicroenv

SITE=$1      # e.g. "Maquipucuna"

PARAMS=$2    # path to params RDS (e.g. data/params/p_poll.rds), or omit for literature defaults

EXPTAG=$3

HEIGHT_STEP=$4  # microenv height-tier spacing to run at (must already exist); default 0.25 if omitted

SPECIES_FILE=$5  # optional: path to an RDS of FinalID values (params$species_subset) -- see
                 # run_colonization.R's species_file arg. Omit for every species observed at the site.

Rscript scripts/02_model/run/run_colonization.R "$SITE" "$PARAMS" "$EXPTAG" "$HEIGHT_STEP" "$SPECIES_FILE"


