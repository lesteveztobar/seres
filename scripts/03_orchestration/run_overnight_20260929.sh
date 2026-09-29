#!/bin/bash
# run_overnight_20260929.sh -- full re-run on the engine with K calibrated to
# Alzate-Q et al. 2019 (562.9 Maxillariinae/ha at Maquipucuna) and n_founders
# capped at 500. Priority via --nice (lower = sooner); real dependencies only
# where one job consumes another's output. Every colonisation job waits
# afterok on the K calibration gate (G2).
#   1 A11 height resolution   2 A10 competition/isolation   3 A12 horizontal
#   4 A08/A09 elevation        5 A07 covariation             6 PAWN v9
#   then A09b held-out, results manifest, plots, output organisation, PDF build.
# Usage (repo root, login node -- only submits): bash scripts/03_orchestration/run_overnight_20260929.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

OBS=data/csv/combinedv6.csv
NICHE=data/processed/species_niches.rds
BG=data/processed/niche_background_density.rds
SITES=("Maquipucuna" "Mashpi" "Yanayacu" "MindoMirador" "MindoTarabita" "Saloya")
if [ -n "${SKIP_SITES:-}" ]; then
  _keep=(); for _s in "${SITES[@]}"; do case " $SKIP_SITES " in *" $_s "*) ;; *) _keep+=("$_s");; esac; done
  SITES=("${_keep[@]}")
fi
ENV="--export=ALL,CANOPY_CLIM_MODE=voxel,CANOPY_OBS_CSV=$OBS,CANOPY_NICHE_CACHE=$NICHE,CANOPY_NICHE_BACKGROUND=$BG"
ACCT="--account=ag_biob_scabral"
R_PRE='module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a Miniforge3/24.1.2-0 UDUNITS/2.2.28-GCCcore-13.2.0; export LD_LIBRARY_PATH=/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH; unset PYTHONPATH; cd /home/s38leste_hpc/seres'

mkdir -p logs
MANIFEST="logs/overnight_20260929_$(date +%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=(); MAIN_IDS=()
log_job() { echo "$1: $2" | tee -a "$MANIFEST"; ALL_IDS+=("$2"); }
dep_any() { local o="afterany"; for i in "$@"; do o="$o:$i"; done; echo "$o"; }
sub() { sbatch --parsable $ACCT "$@"; }
# Memory per colonisation job (2026-09-29): past PAWN tasks running 16
# simulations peaked at ~375G; tall-canopy sites (Saloya, MindoMirador) have
# OOM'd near 880G, so they keep the full request. Single-species isolation
# arrays are ~n_species times smaller.
mem_multi() { case "$1" in Saloya|MindoMirador) echo 1600G ;; *) echo 600G ;; esac; }
mem_iso()   { case "$1" in Saloya|MindoMirador) echo 900G ;;  *) echo 300G ;; esac; }


echo "== G2: K calibration gate =="
G2=$(sub $ENV --job-name=k_calibrate --partition=lm_short --time=03:00:00 --cpus-per-task=2 --mem=200G \
  --output=logs/%x_%j.out --wrap="$R_PRE; Rscript scripts/diagnostics/calibrate_k.R")
log_job k_calibrate "$G2"; MAIN_IDS+=("$G2")
GATE="--dependency=afterok:$G2 --kill-on-invalid-dep=yes"

echo "== P1: A11 height resolution =="
for site in "${SITES[@]}"; do
  if [ "$site" = "Maquipucuna" ]; then
    id=$(sub $ENV $GATE --nice=0 --job-name=A11_hres_$site --time=16:00:00 \
      scripts/experiments/A11_height_resolution/run_height_resolution_experiment.sh "$site" "0.1,0.25,0.4,0.5,0.75,1.0" data/params/realistic_75founders.rds)
  else
    mid=$(sub --nice=0 --job-name=microenv_${site}_h1.00 --partition=lm_medium --time=24:00:00 --cpus-per-task=32 --mem=900G \
      --output=logs/%x_%j.out --export=ALL,CANOPY_OBS_CSV=$OBS \
      --wrap="$R_PRE; export CANOPY_PYTHON=/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python; export LD_PRELOAD=/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3; export CANOPY_SCRATCH=\$(ws_allocate seres 90); mkdir -p \$CANOPY_SCRATCH/tmp; export TMPDIR=\$CANOPY_SCRATCH/tmp; Rscript scripts/01_microclimate/run_microclimate_site.R $site 12 1.0")
    log_job "microenv_${site}_h1.00" "$mid"; MAIN_IDS+=("$mid")
    id=$(sub $ENV --dependency=afterok:$G2:$mid --kill-on-invalid-dep=yes --nice=0 --job-name=A11_hres_$site --time=16:00:00 \
      scripts/experiments/A11_height_resolution/run_height_resolution_experiment.sh "$site" "0.4,1.0" data/params/realistic_75founders.rds)
  fi
  log_job "A11_hres_$site" "$id"; MAIN_IDS+=("$id")
done

echo "== P2: A10 realistic baselines (also the A09b input) + isolation + analysis =="
BASE_IDS=(); ISO_IDS=()
for site in "${SITES[@]}"; do
  id=$(sub $ENV $GATE --nice=100 --mem=$(mem_multi $site) --job-name=A10_base_$site scripts/02_model/run/run_colonization.sh "$site" data/params/realistic.rds realistic 0.4)
  log_job "A10_base_$site" "$id"; BASE_IDS+=("$id"); MAIN_IDS+=("$id")
done
{
  read -r _header
  while IFS=, read -r site species species_file exp_tag _gen; do
    site=${site//\"/}; species_file=${species_file//\"/}; exp_tag=${exp_tag//\"/}
    [ "$site" = "Yanayacu" ] && continue   # single taxon unit: no competitors to isolate from
    keep=0; for s in "${SITES[@]}"; do [ "$s" = "$site" ] && keep=1; done
    [ "$keep" = 1 ] || continue
    id=$(sub $ENV $GATE --nice=100 --mem=$(mem_iso $site) --job-name=A10_iso_${site}_${exp_tag#isolation_} \
      scripts/02_model/run/run_colonization.sh "$site" data/params/realistic.rds "$exp_tag" 0.4 "$species_file")
    log_job "A10_iso_${site}_${exp_tag}" "$id"; ISO_IDS+=("$id"); MAIN_IDS+=("$id")
  done
} < data/params/isolation_manifest.csv
id=$(sub $ENV --dependency="$(dep_any "${BASE_IDS[@]}" "${ISO_IDS[@]}")" --nice=100 --job-name=A10_competition_analysis \
  --partition=lm_short --time=01:00:00 --mem=200G --output=logs/%x_%j.out \
  --wrap="$R_PRE; Rscript scripts/experiments/A10_competition_isolation/competition_analysis.R realistic")
log_job A10_competition_analysis "$id"; MAIN_IDS+=("$id")

echo "== P3: A12 horizontal resolution (pooled vs voxel) =="
for site in "${SITES[@]}"; do
  for mode in pooled voxel; do
    id=$(sub --export=ALL,CANOPY_CLIM_MODE=$mode,CANOPY_OBS_CSV=$OBS,CANOPY_NICHE_CACHE=$NICHE,CANOPY_NICHE_BACKGROUND=$BG \
      $GATE --nice=200 --mem=$(mem_multi $site) --job-name=A12_${mode}_$site \
      scripts/02_model/run/run_colonization.sh "$site" data/params/realistic_75founders.rds horizontal_resolution_$mode 0.4)
    log_job "A12_${mode}_$site" "$id"; MAIN_IDS+=("$id")
  done
done

echo "== P4: A08 / A09 elevation gradient (climate only, no gate) =="
id=$(sub $ENV --nice=300 --job-name=A08_climate_variation --partition=lm_medium --time=08:00:00 --cpus-per-task=32 --mem=900G \
  --output=logs/%x_%j.out \
  --wrap="$R_PRE; Rscript scripts/experiments/A08_climate_variation_between_sites/climate_variation_test.R && Rscript scripts/experiments/A08_climate_variation_between_sites/climate_variation_between_sites.R")
log_job A08_climate_variation "$id"; MAIN_IDS+=("$id")
id=$(sub $ENV --nice=300 scripts/experiments/A09_elevation_climate_exchange/run_climate_variation_relative_height.sh)
log_job A09_relative_height "$id"; MAIN_IDS+=("$id")
id=$(sub $ENV --nice=300 --job-name=A09_elevation_bins --partition=lm_medium --time=08:00:00 --cpus-per-task=16 --mem=800G \
  --output=logs/%x_%j.out \
  --wrap="$R_PRE; Rscript scripts/experiments/A09_elevation_climate_exchange/step1_step2_v2.R && D=\$(ls -td output/microclimate_reanalysis_* | head -1) && Rscript scripts/experiments/A09_elevation_climate_exchange/step2_elevation_bins_v2.R \"\$D\"")
log_job A09_elevation_bins "$id"; MAIN_IDS+=("$id")

echo "== P5: A07 covariation =="
id=$(sub $ENV --nice=400 scripts/experiments/A07_climate_decoupling/run_climate_decoupling.sh)
log_job A07_climate_decoupling "$id"; MAIN_IDS+=("$id")

echo "== P6: PAWN global sensitivity (lowest priority; analysed with whatever completes) =="
# Design is fully described by these settings plus the ranges file; the same
# parameter sets and seeds are run at every site.
PAWN_TAG="${PAWN_TAG:-v9}"; PAWN_POINTS=500; PAWN_REPS=3; PAWN_SEED=20260929; PAWN_TASKS=50
PAWN_RANGES="data/params/sensitivity_ranges_${PAWN_TAG}.csv"
PAWN_DIR=scripts/experiments/A05_sensitivity_v8
pawn_mem() { case "$1" in Mashpi) echo 700G ;; Saloya|MindoMirador) echo 450G ;; *) echo 350G ;; esac; }
SITES_CSV=$(IFS=,; echo "${SITES[*]}")
LHS=$(sub $ENV --nice=1000 --job-name=pawn_${PAWN_TAG}_design --partition=intelsr_short --time=00:15:00 --mem=8G \
  --output=logs/%x_%j.out \
  --wrap="$R_PRE; Rscript $PAWN_DIR/build_lhs_design.R $PAWN_TAG $SITES_CSV $PAWN_POINTS $PAWN_REPS $PAWN_RANGES $PAWN_SEED")
log_job "pawn_${PAWN_TAG}_design" "$LHS"
PAWN_ARRS=()
for site in "${SITES[@]}"; do
  a=$(sub --dependency=afterok:$G2:$LHS --kill-on-invalid-dep=yes --nice=1000 --mem=$(pawn_mem $site) \
    --job-name=pawn_${PAWN_TAG}_$site --array=1-${PAWN_TASKS}%6 \
    --export=ALL,CANOPY_OBS_CSV=$OBS,CANOPY_NICHE_CACHE=$NICHE,CANOPY_NICHE_BACKGROUND=$BG,CANOPY_SENS_DESIGN=data/params/sensitivity_design_${PAWN_TAG}_${site}.rds,CANOPY_SENS_OUTDIR=sensitivity_runs_${PAWN_TAG}_${site},CANOPY_ARRAY_N=$PAWN_TASKS \
    $PAWN_DIR/run_design_array.sh)
  log_job "pawn_${PAWN_TAG}_array_$site" "$a"; PAWN_ARRS+=("$a")
done
PAWN=$(sub --dependency="$(dep_any "${PAWN_ARRS[@]}")" --nice=1000 --job-name=pawn_${PAWN_TAG}_analysis \
  $PAWN_DIR/run_pawn_analysis.sh "$PAWN_TAG" "$PAWN_RANGES")
log_job "pawn_${PAWN_TAG}_analysis" "$PAWN"

echo "== A09b held-out validation + crossing height (after the realistic runs) =="
id=$(sub $ENV --dependency="$(dep_any "${BASE_IDS[@]}")" --nice=0 \
  scripts/experiments/A09b_held_out_validation/run_held_out_validation.sh realistic)
log_job A09b_held_out "$id"; MAIN_IDS+=("$id")

echo "== End of chain: manifest, plots, organise, PDF, mail =="
END_DEP="$(dep_any "${MAIN_IDS[@]}")"
M=$(sub $ENV --dependency="$END_DEP" scripts/experiments/A19_results_manifest_build/run_build_results_manifest_data.sh); log_job results_manifest "$M"
P=$(sub $ENV --dependency="$END_DEP" scripts/diagnostics/run_plots.sh); log_job plots "$P"
N=$(sub $ENV --dependency="$END_DEP" scripts/experiments/A17_niche_suitability_plots/run_niche_figures.sh); log_job niche_figures "$N"
O=$(sub --dependency=afterany:$M:$P:$N --job-name=organize_output --partition=intelsr_short --time=00:10:00 --mem=2G \
  --output=logs/%x_%j.out --wrap="cd /home/s38leste_hpc/seres && bash scripts/diagnostics/organize_output.sh --apply"); log_job organize_output "$O"
B=$(sub --dependency=afterany:$O --job-name=build_pdf --partition=intelsr_short --time=00:30:00 --mem=8G \
  --output=logs/%x_%j.out --wrap="sh /home/s38leste_hpc/seres/report/build.sh"); log_job build_pdf "$B"
NONPAWN=(); for i in "${ALL_IDS[@]}"; do case " $LHS ${PAWN_ARRS[*]} $PAWN " in *" $i "*) ;; *) NONPAWN+=("$i");; esac; done
MAIL=$(sub --dependency="$(dep_any "${NONPAWN[@]}")" --job-name=overnight_summary --mail-type=END --mail-user=s38leste@uni-bonn.de \
  --output=logs/%x_%j.out --error=logs/%x_%j.err scripts/diagnostics/batch_summary_mailer.sh overnight_20260929 "${NONPAWN[@]}")
log_job summary_mail "$MAIL"

echo
echo "Submitted $(wc -l < "$MANIFEST") entries -- job map: $MANIFEST"
echo "Watch: squeue -u \$USER -o '%.10i %.30j %.8T %.6y %.20E'"
