#!/bin/bash
# batch_exp.sh — submit one sensitivity-experiment job per site to the Marvin cluster
# Each job sweeps one parameter (5 values, handled inside run_colonization.R)
# across all timesteps for one site. SITES x EXP = 5 x 9 = 45 jobs.
#
# RECOMMENDED ORDER: run founder_number for one site alone first (e.g.
#   sbatch scripts/02_model/run/run_colonization.sh Maquipucuna "$(pwd)/data/params/n_founders.rds" founder_number
# ) to find an n_founders value that actually persists — the fecundity math
# means every other sweep is likely to come back extinct until that's
# fixed. Once found, update N_FOUNDERS_DEFAULT in make_params.R, re-run it,
# and only then submit the rest (including via this script).
#
# run_colonization.R also takes an optional 4th arg, height_step
# (microenv resolution to run at — must already exist), e.g.
#   sbatch scripts/02_model/run/run_colonization.sh Maquipucuna "$(pwd)/data/params/n_founders.rds" founder_number 0.25
# Defaults to 0.25 (production resolution) if omitted — every job below runs
# at that default unless this script is edited to pass a 4th argument.
#
# CAUTION: reproduction_factorial is a 125-combo factorial (p_poll x p_germ x
# p_s1), much heavier than the other 5-value one-at-a-time sweeps — consider
# running it for one site directly via run_colonization.sh before including
# it in a full SITES x EXP batch_exp.sh submission.
#
# Usage: sbatch batch_exp.sh
# Run from: /home/s38leste_hpc/canopymicroenv/
#
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --output=/home/s38leste_hpc/canopymicroenv/logs/log_array_%j.out

module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a

Rscript scripts/02_model/setup/make_params.R

# Master observations CSV -- mirrors paths.R's OBSERVATIONS_CSV (same
# CANOPY_OBS_CSV override, same default). SITES is derived from it rather
# than hardcoded, so a newly added site is picked up automatically.
CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combined_with_identification.csv}"
export CANOPY_OBS_CSV
mapfile -t SITES < <(awk -F, 'NR>1 && $2!="" {print $2}' "$CANOPY_OBS_CSV" | sort -u)

# EXP and PARAMS are paired by index — PARAMS[k] is the RDS make_params.R
# built for the experiment EXP[k]. Keep these two arrays in sync.
EXP=("pollination_success" "adult_survival_intercept" "germination_probability" "reproduction_cost" "climate_sensitivity_rh" "precipitation_sensitivity" "founder_number" "reproduction_factorial")
PARAMS=("p_poll.rds" "beta0A.rds" "p_germ.rds" "cost_repro.rds" "beta_rh.rds" "beta_precip.rds" "n_founders.rds" "reproduction_factorial.rds")

PARAMS_DIR="/home/$USER/canopymicroenv/data/params"

for i in "${!SITES[@]}"; do
    for k in "${!EXP[@]}"; do
        sbatch /home/$USER/canopymicroenv/scripts/02_model/run/run_colonization.sh \
            "${SITES[$i]}" \
            "$PARAMS_DIR/${PARAMS[$k]}" \
            "${EXP[$k]}"
    done
done
