#!/bin/bash
# run_full_resolution_pipeline.sh — runs the full sequence unattended:
#   1. height_res_array.sh Maquipucuna       (0.1/0.25/0.5/1.0m, sequential)
#   2. run_height_resolution_experiment.sh    (real colonization comparison,
#                                               depends on step 1 finishing)
#   3. microenv_array.sh for the other 6 sites (excludes Maquipucuna, since
#                                                 step 1 already covers it;
#                                                 depends on step 2 finishing
#                                                 so nothing writes to scratch
#                                                 concurrently with step 1/2)
#
# Each stage is chained via --dependency=afterok against the previous
# stage's LAST job ID, so this whole pipeline can be submitted once and left
# running overnight unattended -- nothing here blocks or waits interactively.
#
# Usage: bash run_full_resolution_pipeline.sh
# Run from: /home/s38leste_hpc/seres/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"

echo "== Stage 1: height_res_array.sh Maquipucuna =="
stage1_out=$(bash /home/s38leste_hpc/seres/scripts/experiments/A11_height_resolution/height_res_array.sh Maquipucuna 2>&1 | tee /dev/stderr)
stage1_id=$(echo "$stage1_out" | grep "LAST_JOB_ID=" | tail -1 | cut -d= -f2)
if [ -z "$stage1_id" ]; then
  echo "ERROR: could not capture stage 1's last job ID -- aborting pipeline."
  exit 1
fi
echo "Stage 1 last job: $stage1_id"

echo "== Stage 2: run_height_resolution_experiment.sh (depends on $stage1_id) =="
stage2_id=$(sbatch --parsable --dependency=afterok:${stage1_id} \
  /home/s38leste_hpc/seres/scripts/experiments/A11_height_resolution/run_height_resolution_experiment.sh \
  Maquipucuna "0.1,0.25,0.5,1.0" data/params/realistic.rds)
echo "Stage 2 job: $stage2_id"

echo "== Stage 3: microenv_array.sh, all sites except Maquipucuna (depends on $stage2_id) =="
# microenv_array.sh is itself a short dispatcher that submits its own chain
# of real jobs -- we can't directly sbatch --dependency on ITS internal
# chain from here without it already having run. So stage 3 is submitted as
# a dependent dispatcher: it won't start executing (and therefore won't
# submit any real site jobs) until stage 2 completes.
stage3_id=$(sbatch --parsable --dependency=afterok:${stage2_id} \
  /home/s38leste_hpc/seres/scripts/01_microclimate/microenv_array.sh 12 0.25 "" Maquipucuna)
echo "Stage 3 dispatcher job: $stage3_id (will submit the real 6-site chain once it runs)"

echo ""
echo "Pipeline submitted. Check progress with: squeue --me"
echo "Stage 1 (Maquipucuna resolution) -> Stage 2 (colonization comparison) -> Stage 3 (other 6 sites, 0.25m)"