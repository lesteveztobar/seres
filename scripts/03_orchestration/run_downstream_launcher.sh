#!/bin/bash
#SBATCH --job-name=downstream_launcher
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=s38leste@uni-bonn.de
#SBATCH --output=/home/s38leste_hpc/seres/logs/%x_%j.out
#SBATCH --error=/home/s38leste_hpc/seres/logs/%x_%j.err
# Runs ONLY after the sanity gate passed (afterok chain): submits everything
# downstream of MindoMirador + the rebuilt niche cache (2026-09-24).
set -euo pipefail
cd /home/$USER/seres
export CANOPY_OBS_CSV=data/csv/combinedv6.csv
# Saloya left out for now (few observations: 17 raw / 5 identified); set SKIP_SITES="" to include it.
export SKIP_SITES="${SKIP_SITES-Saloya}"
export CANOPY_SKIP_SITES="$SKIP_SITES"
D=scripts/experiments/A05_sensitivity_v8
MAP=logs/downstream_$(date +%Y%m%d_%H%M%S).jobs; : > $MAP
sub() { local name=$1; shift; local id; id=$(sbatch --parsable "$@"); echo "$name: $id" | tee -a $MAP; }
arr() { # name design outdir ntasks
  sub "$1" --job-name="$1" --array=1-$4%12 --export=ALL,CANOPY_SENS_DESIGN=$2,CANOPY_SENS_OUTDIR=$3,CANOPY_ARRAY_N=$4 $D/run_design_array.sh; }
# --- Phase F (5 approved + horizontal resolution) ---
bash scripts/03_orchestration/run_phase_f_approved.sh   | tee -a $MAP
bash scripts/experiments/A12_horizontal_resolution/run_horizontal_resolution.sh | tee -a $MAP
# --- re-centred factorial, K sweep, clean sensitivity re-run at final inputs ---
arr factorial_v8      data/params/factorial_v8_design.rds              factorial_v8_runs      75
arr ksweep            data/params/ksweep_design.rds                    ksweep_runs            14
arr sens_rerun_s1     data/params/sensitivity_design.rds               sens_rerun_stage1      60
arr sens_rerun_s2     data/params/sensitivity_design_stage2_full.rds   sens_rerun_stage2      60
# --- Phase D climate analyses (2026-09-24): these read only microclimate
# manifests, not the niche cache or any colonization output, so they don't
# need to wait for the niche rebuild/sanity gate this launcher itself is
# gated on -- only for MindoMirador + LaElenita's Phase A reruns. Defaults
# to the actual in-flight job IDs (this launcher, 27876675, was already
# submitted before this edit landed, so it can't pick up new --export vars
# at submission time); CLIM_DEP1/CLIM_DEP2 can override for a later manual
# resubmission once both sites are confirmed fresh.
CLIM_DEP1="${CLIM_DEP1:-27876466}"
CLIM_DEP2="${CLIM_DEP2:-27876667}"
CLIM_GATE="--dependency=afterok:${CLIM_DEP1}:${CLIM_DEP2}"
sub climate_decoupling $CLIM_GATE --export=ALL,EXCLUDE_SALOYA=1 scripts/experiments/A07_climate_decoupling/run_climate_decoupling.sh
sub climvar_relheight  $CLIM_GATE --export=ALL,CANOPY_SKIP_SITES=Saloya scripts/experiments/A09_elevation_climate_exchange/run_climate_variation_relative_height.sh
echo "job map: $MAP"
