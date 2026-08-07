#!/bin/bash
# run_pipeline.sh — menu-driven launcher for the full canopymicroenv HPC
# pipeline (microenv -> params -> niche -> experiments -> plots).
#
# Chains stages with SLURM job dependencies (--dependency=afterok) instead of
# polling, so the WHOLE CHAIN is submitted in a few seconds and then runs
# unattended overnight: you don't need to keep a terminal/process open, and
# you don't need to come back and resubmit the next step by hand. If a step
# fails, everything chained after it is auto-cancelled by SLURM
# (DependencyNeverSatisfied) instead of quietly running on bad inputs.
#
# Usage:
#   sh scripts/03_orchestration/run_pipeline.sh                                  # interactive prompt
#   sh scripts/03_orchestration/run_pipeline.sh microenv params niche experiments plots
#   sh scripts/03_orchestration/run_pipeline.sh params, niche, plots             # commas ok too
#   sh scripts/03_orchestration/run_pipeline.sh -y params niche plots            # skip confirmation
#
# Steps run in exactly the order you give them — this script does not
# reorder them for you, so make sure the order makes sense (e.g. microenv
# before niche/experiments, since both read microenv_<site>.rds).
#
# Env overrides (defaults match the rest of the pipeline's production settings):
#   N_MONTHS=12       months of ERA5 for microenv (see run_microenv.sh)
#   HEIGHT_STEP=0.25  microenv height-tier spacing used by every stage below
#                     (must already exist for niche/experiments if microenv
#                     isn't part of this run)
#
# Run from: /home/s38leste_hpc/canopymicroenv/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/03_orchestration/lib.sh

# Master observations CSV -- mirrors paths.R's OBSERVATIONS_CSV (same
# CANOPY_OBS_CSV override, same default). Keep the two in sync if either
# changes. SITES is derived from it rather than hardcoded, so a newly added
# site is picked up automatically.
CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combined_with_identification.csv}"
export CANOPY_OBS_CSV
mapfile -t SITES < <(awk -F, 'NR>1 && $2!="" {print $2}' "$CANOPY_OBS_CSV" | sort -u)
# EXP and PARAMS are paired by index, same pairing as batch_exp.sh
EXP=("pollination_success" "adult_survival_intercept" "germination_probability" "reproduction_cost" "climate_sensitivity_rh" "precipitation_sensitivity" "founder_number" "reproduction_factorial")
PARAMS=("p_poll.rds" "beta0A.rds" "p_germ.rds" "cost_repro.rds" "beta_rh.rds" "beta_precip.rds" "n_founders.rds" "reproduction_factorial.rds")
PARAMS_DIR="$(pwd)/data/params"

N_MONTHS=${N_MONTHS:-12}
HEIGHT_STEP=${HEIGHT_STEP:-0.25}

VALID_STEPS=("microenv" "params" "niche" "experiments" "plots")

# ── parse args ─────────────────────────────────────────────────────────────
AUTO_YES=0
ARGS=()
for arg in "$@"; do
  case "$arg" in
    -y|--yes) AUTO_YES=1 ;;
    *) ARGS+=("$arg") ;;
  esac
done

if [ "${#ARGS[@]}" -eq 0 ]; then
  echo "canopymicroenv pipeline — available steps: ${VALID_STEPS[*]}"
  echo "Enter the steps you want to run, in order (space or comma separated):"
  read -r -p "> " REPLY
  ARGS=("$REPLY")
fi

RAW="${ARGS[*]//,/ }"
read -r -a STEPS <<< "$RAW"

if [ "${#STEPS[@]}" -eq 0 ]; then
  echo "No steps given, exiting."
  exit 1
fi

for s in "${STEPS[@]}"; do
  ok=0
  for v in "${VALID_STEPS[@]}"; do [ "$s" = "$v" ] && ok=1; done
  if [ "$ok" -eq 0 ]; then
    echo "Unknown step: '$s' — valid steps are: ${VALID_STEPS[*]}"
    exit 1
  fi
done

echo
echo "Plan (in order): ${STEPS[*]}"
echo "  height_step=$HEIGHT_STEP  n_months=$N_MONTHS"
echo "  microenv would submit ${#SITES[@]} jobs; experiments would submit $(( ${#SITES[@]} * ${#EXP[@]} + 1 )) jobs (incl. its own params regen)"
if [ "$AUTO_YES" -ne 1 ]; then
  read -r -p "Submit this chain to SLURM now? [y/N] " CONFIRM
  [[ "$CONFIRM" =~ ^[Yy]$ ]] || { echo "Aborted, nothing submitted."; exit 0; }
fi

mkdir -p logs
MANIFEST="logs/pipeline_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
DEP=()   # job IDs the *next* step must wait on

for step in "${STEPS[@]}"; do
  DEP_ARG=()
  if [ "${#DEP[@]}" -gt 0 ]; then
    DEP_ARG=(--dependency="$(join_dep "${DEP[@]}")")
  fi
  NEW_IDS=()

  case "$step" in

    microenv)
      if ls data/processed/microenv_*.rds >/dev/null 2>&1; then
        echo "NOTE: existing microenv_*.rds manifests found — per microenv_array.sh," \
             "if these predate the current canopy-height ceiling, move or delete them" \
             "first or this step will see them and exit without regenerating anything."
      fi
      for site in "${SITES[@]}"; do
        id=$(sbatch --parsable "${DEP_ARG[@]}" scripts/01_microclimate/run_microenv.sh "$site" "$N_MONTHS" "$HEIGHT_STEP")
        echo "  microenv/$site -> job $id"
        NEW_IDS+=("$id")
      done
      ;;

    params)
      id=$(sbatch --parsable "${DEP_ARG[@]}" --partition=intelsr_short \
             --account=ag_biob_scabral --time=00:10:00 --ntasks=1 \
             --output=logs/log_%j.out \
             --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; Rscript scripts/02_model/setup/make_params.R")
      echo "  params -> job $id"
      NEW_IDS+=("$id")
      ;;

    niche)
      id=$(sbatch --parsable "${DEP_ARG[@]}" scripts/02_model/setup/characterize_niches.sh "$HEIGHT_STEP")
      echo "  niche -> job $id"
      NEW_IDS+=("$id")
      ;;

    experiments)
      # mirrors batch_exp.sh's own SITES x EXP loop (which regenerates params
      # first), but submitted job-by-job here so we get real per-job IDs to
      # chain "plots" on instead of the array-submitter's own short-lived job id
      pid=$(sbatch --parsable "${DEP_ARG[@]}" --partition=intelsr_short \
             --account=ag_biob_scabral --time=00:10:00 --ntasks=1 \
             --output=logs/log_%j.out \
             --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; Rscript scripts/02_model/setup/make_params.R")
      echo "  experiments/params -> job $pid"
      for i in "${!SITES[@]}"; do
        for k in "${!EXP[@]}"; do
          id=$(sbatch --parsable --dependency="afterok:$pid" scripts/02_model/run/run_colonization.sh \
                 "${SITES[$i]}" "$PARAMS_DIR/${PARAMS[$k]}" "${EXP[$k]}" "$HEIGHT_STEP")
          echo "  experiments/${SITES[$i]}/${EXP[$k]} -> job $id"
          NEW_IDS+=("$id")
        done
      done
      ;;

    plots)
      id=$(sbatch --parsable "${DEP_ARG[@]}" scripts/02_model/plots/run_plots.sh)
      echo "  plots -> job $id"
      NEW_IDS+=("$id")
      ;;
  esac

  echo "$step: ${NEW_IDS[*]}" >> "$MANIFEST"
  DEP=("${NEW_IDS[@]}")
done

echo
echo "All steps submitted — job map saved to $MANIFEST"
echo "You can close this terminal now; SLURM will run the chain overnight."
echo "Check progress any time with: squeue -u \$USER"
echo "Or after the fact with:       sacct -j <jobid1>,<jobid2>,... --format=JobID,JobName,State,Elapsed"
