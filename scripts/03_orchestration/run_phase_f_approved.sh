#!/bin/bash
# run_phase_f_approved.sh — Phase F, the 5 approved experiments (2026-09-08):
# height resolution, persistence, founder number, competition and isolation,
# horizontal resolution. Reproduction-survival factorial held pending
# separate approval; held-out validation dropped (confirmed post-processing,
# depends on best_combo, which depends on the held factorial).
#
# Every job: CANOPY_CLIM_MODE=voxel (the new default anyway, set explicitly
# for the record), CANOPY_OBS_CSV=combinedv6.csv, height_step=0.4, the 6
# modelled sites (LaElenita excluded, no confirmed species).
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/diagnostics/lib.sh

export CANOPY_OBS_CSV="data/csv/combinedv6.csv"
SITES=("Maquipucuna" "Mashpi" "Yanayacu" "MindoMirador" "MindoTarabita" "Saloya")
# SKIP_SITES (space-separated, default empty): drop sites from this run, e.g. SKIP_SITES="Saloya"
if [ -n "${SKIP_SITES:-}" ]; then
  _keep=(); for _s in "${SITES[@]}"; do case " $SKIP_SITES " in *" $_s "*) echo "SKIPPING $_s (SKIP_SITES)";; *) _keep+=("$_s");; esac; done
  SITES=("${_keep[@]}")
fi
SBATCH_ENV="--export=ALL,CANOPY_CLIM_MODE=voxel,CANOPY_OBS_CSV=data/csv/combinedv6.csv"

mkdir -p logs
MANIFEST="logs/phase_f_approved_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=()
log_job() { echo "$1: $2" | tee -a "$MANIFEST"; ALL_IDS+=("$2"); }

# 2026-09-24: persistence:realistic is now ALSO the competition/isolation
# baseline (30-founder realistic.rds, per the user's explicit decision to
# switch off the old 273-founder baseline) -- submitted first, on its own,
# so it can be sanity-checked (recruited/extinct/finite totals) before
# competition/isolation's dependency is treated as trustworthy, not just
# wall-clock-satisfied. Its job IDs are reused directly below, not
# resubmitted as a separate "competition_baseline" run.
echo "== Persistence: best_case + realistic, every site (realistic doubles as the competition/isolation baseline) =="
REALISTIC_IDS=()
for site in "${SITES[@]}"; do
  id=$(sbatch --parsable $SBATCH_ENV scripts/02_model/run/run_colonization.sh "$site" data/params/best_case.rds best_case 0.4)
  log_job "persistence_best_case_${site}" "$id"
  id=$(sbatch --parsable $SBATCH_ENV scripts/02_model/run/run_colonization.sh "$site" data/params/realistic.rds realistic 0.4)
  log_job "persistence_realistic_${site}" "$id"
  REALISTIC_IDS+=("$id")
done

echo
echo "== Founder number: 20-level sweep, every site =="
for site in "${SITES[@]}"; do
  id=$(sbatch --parsable $SBATCH_ENV scripts/02_model/run/run_colonization.sh "$site" data/params/n_founders.rds founder_number 0.4)
  log_job "founder_number_${site}" "$id"
done

echo
echo "== Competition and isolation: realistic (30-founder) baseline [= persistence_realistic above] + one run per (site, species) =="
ISO_IDS=()
{
  read -r _header
  while IFS=, read -r site species species_file exp_tag; do
    site=$(echo "$site" | tr -d '"')
    species_file=$(echo "$species_file" | tr -d '"')
    exp_tag=$(echo "$exp_tag" | tr -d '"')
    # LaElenita excluded (isolation manifest may still list it if it has any
    # confirmed species -- it doesn't, so this is defensive, not expected to trigger).
    skip=1
    for s in "${SITES[@]}"; do [ "$s" = "$site" ] && skip=0; done
    [ "$skip" = 1 ] && continue
    id=$(sbatch --parsable $SBATCH_ENV scripts/02_model/run/run_colonization.sh "$site" data/params/realistic.rds "$exp_tag" 0.4 "$species_file")
    log_job "isolation_${site}_${exp_tag}" "$id"
    ISO_IDS+=("$id")
  done
} < data/params/isolation_manifest.csv
join_dep_any() { local out="afterany"; for id in "$@"; do out="$out:$id"; done; echo "$out"; }
COMPETITION_ID=$(sbatch --parsable --dependency="$(join_dep_any "${REALISTIC_IDS[@]}" "${ISO_IDS[@]}")" \
  --partition=intelsr_short --account=ag_biob_scabral --time=00:30:00 --ntasks=1 --output=logs/log_%j.out \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; export CANOPY_OBS_CSV=data/csv/combinedv6.csv; Rscript scripts/experiments/A10_competition_isolation/competition_analysis.R realistic")
log_job "competition_analysis" "$COMPETITION_ID"

echo
echo "== Height resolution: every site =="
for site in "${SITES[@]}"; do
  id=$(sbatch --parsable $SBATCH_ENV --time=16:00:00 \
    scripts/experiments/A11_height_resolution/run_height_resolution_experiment.sh \
    "$site" "0.1,0.25,0.5,1.0" data/params/realistic_75founders.rds)
  log_job "height_resolution_${site}" "$id"
done

SUMMARY_ID=$(sbatch --parsable --dependency="$(join_dep_any "${ALL_IDS[@]}")" \
  --job-name=phase_f_approved_summary --mail-type=END --mail-user=s38leste@uni-bonn.de \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  scripts/diagnostics/batch_summary_mailer.sh "phase_f_approved" "${ALL_IDS[@]}")
log_job "summary" "$SUMMARY_ID"

ALL_IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")
echo
echo "Phase F (5 approved experiments) submitted: $(echo "${ALL_IDS[@]}" | wc -w) jobs -- job map: $MANIFEST"
echo "sacct -j $ALL_IDS_CSV --format=JobID,JobName,State,Elapsed"
