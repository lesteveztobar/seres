#!/bin/bash
# run_full_analysis_pipeline.sh — everything needed for the thesis's 4 main
# results (resolution, factorial parameter importance, competition/isolation,
# elevation), chained via SLURM --dependency so the whole thing runs
# unattended: submit it and walk away.
#
# Site list is derived dynamically from combinedv3.csv (Area_or_Site column)
# at submission time -- NOT hardcoded -- so more observations/species at an
# existing site, or an entirely new site, are picked up automatically the
# next time this runs (2026-07-18: confirmed at least one new site is being
# added while this was written). make_sites() (lib.R), which
# the whole microclimate pipeline is built on, was already fully dynamic in
# this same way -- what wasn't dynamic was everything downstream of it in
# this pipeline, which is what this version fixes.
#
# ── What this runs, and why in this order ──────────────────────────────────
#   0. (synchronous, on the login node -- reads a <100-row CSV, instant)
#      make_isolation_species_files.R -- writes isolation_manifest.csv, which
#      this script itself needs to read below to know which isolation jobs
#      to submit. Everything else is a real SLURM job.
#  -1. Base production microenv (0.25m) for any site that doesn't have one
#      yet -- i.e. a genuinely new site. Runs ALONE, before that site's
#      height-step diagnostic variants (step 3): run_microclimate_site.R
#      caches weather/DTM/point-model data per site and reuses it for later
#      height-step variants ONLY if that cache is already populated (see
#      height_res_array.sh's own header note) -- for a brand-new site
#      nothing is cached yet, so submitting the 0.25m base run and the
#      0.1/0.5/1.0 diagnostic variants all at once would race multiple jobs
#      against the same ERA5/DTM cache files. Existing sites skip this step
#      entirely (cache already warm, exactly as before).
#      NOTE: first-time ERA5/DTM acquisition timing is unpredictable
#      (Copernicus CDS queue times especially, outside our control) -- if
#      this times out, resubmit by hand with more time, same as any other
#      timeout in this pipeline.
#   1. characterize_niches.R + make_params.R -- niche cache and params/*.rds
#      regenerated from the CURRENT combinedv3.csv, which now has the
#      "Maxillaria bradei" / "Maxillaria  bradei" (double-space typo) merge
#      fixed (2026-07-18) -- Mashpi and MindoTarabita's species counts
#      changed, so anything species-dependent downstream needs this first.
#      Also waits on step -1 for any new site, since characterize_niches.R
#      reads every site's microenv.
#   2. climate_variation_test.R -- independent of the niche fix (pure
#      climate data). Extends the existing plot_temperature_profile()
#      height-tier test (temperature only, Maquipucuna only) to all 3 niche
#      variables and every site, plus the new cross-site elevation
#      comparison (uses each site's DTM via elevation_helpers.R to fill
#      MindoTarabita/Yanayacu, which have no recorded field elevation at
#      all). Only waits on step -1 (needs each site's microenv to exist),
#      not on the niche fix.
#   3. Missing height-step microenv variants (0.1/0.5/1.0m; 0.25m already
#      exists for sites that had it before this run) -- for a new site,
#      chained after its own step -1 base run; for existing sites, no such
#      dependency (cache already warm). Dominant cost of this whole
#      pipeline: ~2h51m per (site, step) for an already-cached site, based
#      on Maquipucuna's own generation history -- a brand-new site's own
#      base run (step -1) is uncached and could take meaningfully longer.
#   4. realistic_273founders.rds multi-species run, every site (rerun, not
#      reused -- see prof-feedback note below on why 273 founders, and why
#      existing realistic_h0.25 results for Mashpi/MindoTarabita are stale
#      relative to step 1's fixed species list anyway). Depends on step 1.
#   5. Isolation runs, one per (site, species) pair -- each site's species
#      run ALONE (species_subset, see run_colonization.R), same
#      realistic_273founders.rds params as step 4 so species_subset is the
#      only thing that differs. Depends on step 1.
#   6. competition_analysis.R -- compares each species' realized height
#      distribution alone (step 5) vs. with competitors (step 4). Depends on
#      both (afterany, not afterok -- see join_dep_any() below).
#   7. height_resolution_experiment.R, every site, using
#      realistic_273founders.rds -- each site depends on step 1 plus (for
#      any site that needed it) step 3's height-step variants finishing.
#   8. reproduction_factorial_v3, every viable site, run FRESH -- not a
#      resume. The old MindoTarabita checkpoint (500/625 combos) was built
#      against the pre-confidence-filter, pre-FinalID-fix species list and
#      is no longer valid to resume from (moved aside to
#      data/processed/stale_pre_6site_data/, along with Maquipucuna/Mashpi's
#      now-stale completed factorial results). Depends on step 1.
#      NOTE: MindoMirador (renamed MiradorMindo) and Yanayacu are
#      deliberately excluded -- their best_case run still crashes with
#      "missing value where TRUE/FALSE needed" (0/5 replicates, 2026-07-16,
#      still unresolved), and the factorial's upper parameter levels
#      approach best_case's values, so queuing 625 combinations there would
#      likely just hit the same crash 625x for no benefit. Fix that
#      separately before adding these two.
#   9. Final plots (run_plots.sh) + summarize_all_results.R -- depends on
#      everything above (afterany).
#
# Founders: realistic.rds's literature-default 30 founders produces final
# populations of only ~4-8 individuals under realistic vital rates by t=30
# (DECLINING from the 30 founders, not growing) -- too small/noisy for the
# resolution outcome comparison or the isolation-vs-multispecies height-
# distribution comparison. realistic_273founders.rds (make_params.R) keeps
# every vital rate at its literature value and only raises n_founders to
# 273 (the second level from the reproduction_factorial_v3 sweep, so not an
# arbitrary new number) -- more individuals to observe throughout the run,
# without changing the (declining) dynamics being observed. More timesteps
# was considered instead and rejected: on a declining trajectory, running
# longer tends toward FEWER individuals, not more.
#
# CAVEAT this script can't fully check: an existing site's cached microenv
# only covers the bounding box computed from combinedv3.csv WHEN IT WAS
# GENERATED (see make_sites(), lib.R, pad=0.15 degrees). New
# observations added within that same padded extent are fine; if new
# observations extend meaningfully beyond a site's previous extent, that
# site's microenv should be regenerated too (delete
# data/processed/microenv_<site>_h0.25.rds and this script will treat it as
# a new site next run) -- there's no automatic check for this here.
#
# If anything fails partway through, everything chained after it with a
# strict (afterok) dependency is auto-cancelled instead of quietly running
# on bad/incomplete inputs -- check `sacct` (command printed at the end) to
# see what actually completed, fix the failure, and resubmit just that piece
# by hand rather than rerunning this whole script.
#
# Usage: sh scripts/03_orchestration/run_full_analysis_pipeline.sh
# Run from: /home/s38leste_hpc/canopymicroenv/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")/../.."

# Master observations CSV -- mirrors paths.R's OBSERVATIONS_CSV (same
# CANOPY_OBS_CSV override, same default). Keep the two in sync if either
# changes.
CANOPY_OBS_CSV="${CANOPY_OBS_CSV:-data/csv/combined_with_identification.csv}"
export CANOPY_OBS_CSV

mapfile -t SITES < <(awk -F, 'NR>1 && $2!="" {print $2}' "$CANOPY_OBS_CSV" | sort -u)
if [ "${#SITES[@]}" -eq 0 ]; then
  echo "ERROR: no sites found in $CANOPY_OBS_CSV (Area_or_Site column) -- aborting." >&2
  exit 1
fi
echo "Sites found in $CANOPY_OBS_CSV: ${SITES[*]}"

HEIGHT_STEPS_TO_GENERATE=("0.1" "0.5" "1.0")
N_MONTHS=12
RESOLUTION_PARAMS="data/params/realistic_273founders.rds"
COMPETITION_PARAMS="data/params/realistic_273founders.rds"

mkdir -p logs
MANIFEST="logs/full_pipeline_$(date +%Y%m%d_%H%M%S).jobs"
: > "$MANIFEST"
ALL_IDS=()

join_dep() { local out="afterok"; for id in "$@"; do out="$out:$id"; done; echo "$out"; }
# afterany: for aggregation steps specifically written to tolerate missing
# inputs (competition_analysis.R, run_plots.sh's safe_plot(),
# summarize_all_results.R -- all check file.exists()/skip gracefully rather
# than assuming every upstream job succeeded). Using afterok here instead
# would mean a single failed leg out of many permanently blocks the whole
# aggregation with DependencyNeverSatisfied -- the exact stuck-job problem
# the 2026-07-17 plots run hit.
join_dep_any() { local out="afterany"; for id in "$@"; do out="$out:$id"; done; echo "$out"; }
log_job() { echo "$1: $2" | tee -a "$MANIFEST"; ALL_IDS+=("$2"); }

echo "== Step 0: isolation species_subset files (synchronous, instant -- reads \$CANOPY_OBS_CSV) =="
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
Rscript scripts/02_model/setup/make_isolation_species_files.R
ISO_MANIFEST="data/params/isolation_manifest.csv"
if [ ! -f "$ISO_MANIFEST" ]; then
  echo "ERROR: $ISO_MANIFEST was not created -- aborting before submitting anything." >&2
  exit 1
fi

echo
echo "== Step -1: base production microenv (0.25m) for any site that doesn't have one yet =="
declare -A BASE_MICROENV_ID
NEW_SITES=()
for site in "${SITES[@]}"; do
  if [ ! -f "data/processed/microenv_${site}_h0.25.rds" ]; then
    NEW_SITES+=("$site")
    # --partition override required: run_microenv.sh's own #SBATCH default is
    # lm_short, capped at 8h -- requesting 20h there without also overriding
    # partition gets rejected outright by SLURM at submission time (2026-07-22).
    # --cpus-per-task/--mem override added 2026-07-26: without it this falls
    # back to run_microenv.sh's own defaults (4 cpus, 500G) instead of step
    # 3's 32 cpus/900G -- LaElenita's first-time base run got stuck on just 4
    # cores and took ~8x longer than an equivalent 32-core height-step job
    # (Saloya's h0.1, same ~170-500 heights, finished in ~2h on 32 cores).
    # Matching step 3's allocation here so a brand-new site's base run isn't
    # needlessly throttled.
    id=$(sbatch --parsable --partition=lm_medium --time=20:00:00 --cpus-per-task=32 --mem=900G \
      scripts/01_microclimate/run_microenv.sh "$site" "$N_MONTHS" 0.25)
    log_job "microenv_base_${site}" "$id"
    BASE_MICROENV_ID[$site]="$id"
  fi
done
if [ "${#NEW_SITES[@]}" -gt 0 ]; then
  echo "New site(s) detected (no existing 0.25m microenv yet): ${NEW_SITES[*]}"
  echo "First-time ERA5/DTM acquisition timing is unpredictable -- see header note if this times out."
else
  echo "No new sites -- every site already has a 0.25m microenv."
fi

echo
echo "== Step 1: niche re-characterization + params regen =="
NICHE_EXTRA_DEPS=()
for site in "${NEW_SITES[@]}"; do NICHE_EXTRA_DEPS+=("${BASE_MICROENV_ID[$site]}"); done
if [ "${#NICHE_EXTRA_DEPS[@]}" -gt 0 ]; then
  NICHE_ID=$(sbatch --parsable --dependency="$(join_dep "${NICHE_EXTRA_DEPS[@]}")" \
    scripts/02_model/setup/characterize_niches.sh 0.25)
else
  NICHE_ID=$(sbatch --parsable scripts/02_model/setup/characterize_niches.sh 0.25)
fi
log_job "niche" "$NICHE_ID"
PARAMS_ID=$(sbatch --parsable --partition=intelsr_short --account=ag_biob_scabral --time=00:10:00 \
  --ntasks=1 --output=logs/log_%j.out \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; Rscript scripts/02_model/setup/make_params.R")
log_job "params" "$PARAMS_ID"
# NICHE_ID already transitively depends on every new site's base microenv
# (see above), so anything gated on NICHE_PARAMS_DEP doesn't need to also
# list BASE_MICROENV_ID directly -- SLURM dependencies are transitive in
# effect (if C waits on B and B waited on A, C running already implies A
# succeeded).
NICHE_PARAMS_DEP="$(join_dep "$NICHE_ID" "$PARAMS_ID")"

echo
echo "== Step 2: climate variation test (temp/relhum/swdown x height x elevation, every site) =="
CLIMVAR_DEPS=("${BASE_MICROENV_ID[@]}")
if [ "${#CLIMVAR_DEPS[@]}" -gt 0 ]; then
  CLIMVAR_DEP_ARG="--dependency=$(join_dep "${CLIMVAR_DEPS[@]}")"
else
  CLIMVAR_DEP_ARG=""
fi
  # --cpus-per-task=32 added 2026-07-27: lookup_climate_by_height() read each height tier
  # sequentially/single-threaded -- even after bumping --time 02:00:00 ->
  # 06:00:00, the job still hit its own time limit a second time (died on
  # Saloya, 5 sites in). Root cause wasn't the time budget, it was zero
  # parallelism on 100-400+ per-height file reads; climate_variation_test.R's
  # read loop (and build_clim_cache(), get_colonization.R) now parallelize
  # via mclapply, so this needs real cores to use.
  # --mem=900G (was 64G): same 32-worker x ~6GB/height-file OOM as
  # characterize_niches.sh -- see that script's own comment for the math.
CLIMVAR_ID=$(sbatch --parsable $CLIMVAR_DEP_ARG --partition=lm_medium --account=ag_biob_scabral --time=06:00:00 \
  --ntasks=1 --cpus-per-task=32 --mem=900G --output=logs/log_%j.out \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; export LD_LIBRARY_PATH=\"/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:\$LD_LIBRARY_PATH\"; Rscript scripts/02_model/analysis/climate_variation_test.R")
log_job "climate_variation" "$CLIMVAR_ID"

echo
echo "== Step 3: missing height-step microenv variants =="
# HEIGHT_JOB_IDS[site]: space-separated raw job IDs (not pre-joined into an
# "afterok:..." string) so step 7 can freely combine them with NICHE_ID/
# PARAMS_ID via join_dep() at the point of use.
declare -A HEIGHT_JOB_IDS
for site in "${SITES[@]}"; do
  ids=()
  # Chained sequentially per site (each step also depends on the previous
  # step's job, not just the base 0.25m job) -- these all share one ERA5-
  # merged-file cache per site (get_weather(), lib.R) with no
  # locking around its extend/rename step, so running two height-step
  # variants for the same site concurrently can race on that shared file
  # and silently corrupt it: 2026-07-24, three Saloya height-step jobs hit
  # the ERA5 extend path within the same second, and both non-production
  # manifests came back with 0 valid height tiers ("weather[[k]] :
  # subscript out of bounds") despite the SLURM jobs exiting 0. Different
  # SITES' chains still run in parallel with each other -- only steps
  # within the same site's chain are serialized.
  prev_step_id=""
  for step in "${HEIGHT_STEPS_TO_GENERATE[@]}"; do
    suffix="_h${step}"; [ "$step" = "0.1" ] && suffix=""
    if [ -f "data/processed/microenv_${site}${suffix}.rds" ]; then
      continue  # already generated (e.g. Maquipucuna, from earlier work)
    fi
    dep_ids=()
    if [ -n "${BASE_MICROENV_ID[$site]+x}" ]; then dep_ids+=("${BASE_MICROENV_ID[$site]}"); fi
    if [ -n "$prev_step_id" ]; then dep_ids+=("$prev_step_id"); fi
    if [ "${#dep_ids[@]}" -gt 0 ]; then
      dep_arg="--dependency=$(join_dep "${dep_ids[@]}")"
    else
      dep_arg=""
    fi
    # --mem=900G (was 500G): Saloya's h0.1 run (500 height tiers, by far the
    # most of any step here) OOM-killed at 500G after 2h13m on 2026-07-25 --
    # lm_medium nodes have 2TB available (sinfo), so this leaves real margin
    # without meaningfully affecting scheduling (well under node capacity).
    id=$(sbatch --parsable $dep_arg --partition=lm_medium --time=24:00:00 --cpus-per-task=32 --mem=900G \
      --account=ag_biob_scabral --output=logs/log_%j.out \
      scripts/01_microclimate/run_microenv.sh "$site" "$N_MONTHS" "$step")
    log_job "microenv_${site}_h${step}" "$id"
    ids+=("$id")
    prev_step_id="$id"
  done
  if [ "${#ids[@]}" -gt 0 ]; then HEIGHT_JOB_IDS[$site]="${ids[*]}"; fi
done

echo
echo "== Step 4: realistic_273founders.rds multi-species baseline, every site (rerun) =="
declare -A REALISTIC_ID
for site in "${SITES[@]}"; do
  id=$(sbatch --parsable --dependency="$NICHE_PARAMS_DEP" \
    scripts/02_model/run/run_colonization.sh "$site" "$COMPETITION_PARAMS" realistic_273founders 0.25)
  log_job "realistic_${site}" "$id"
  REALISTIC_ID[$site]="$id"
done

echo
echo "== Step 5: isolation runs, one per (site, species) pair =="
ISO_IDS=()
{
  read -r _header
  while IFS=, read -r site species species_file exp_tag; do
    site=$(echo "$site" | tr -d '"')
    species_file=$(echo "$species_file" | tr -d '"')
    exp_tag=$(echo "$exp_tag" | tr -d '"')
    id=$(sbatch --parsable --dependency="$NICHE_PARAMS_DEP" \
      scripts/02_model/run/run_colonization.sh "$site" "$COMPETITION_PARAMS" "$exp_tag" 0.25 "$species_file")
    log_job "isolation_${site}_${exp_tag}" "$id"
    ISO_IDS+=("$id")
  done
} < "$ISO_MANIFEST"

echo
echo "== Step 6: competition analysis (depends on all isolation + realistic runs) =="
COMPETITION_DEPS=("${ISO_IDS[@]}")
for site in "${SITES[@]}"; do COMPETITION_DEPS+=("${REALISTIC_ID[$site]}"); done
COMPETITION_ID=$(sbatch --parsable --dependency="$(join_dep_any "${COMPETITION_DEPS[@]}")" \
  --partition=intelsr_short --account=ag_biob_scabral --time=00:30:00 --ntasks=1 --output=logs/log_%j.out \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; Rscript scripts/02_model/analysis/competition_analysis.R realistic_273founders")
log_job "competition_analysis" "$COMPETITION_ID"

echo
echo "== Step 7: height-resolution outcome test, every site =="
# params_tag mirrors height_resolution_experiment.R's own derivation
# (tools::file_path_sans_ext(basename(params_file))) so this skip check
# matches exactly what that script would name its own output file.
RESOLUTION_PARAMS_TAG="$(basename "$RESOLUTION_PARAMS" .rds)"
RESOLUTION_IDS=()
for site in "${SITES[@]}"; do
  if [ -f "data/processed/height_resolution_${site}_${RESOLUTION_PARAMS_TAG}.rds" ]; then
    echo "  $site: already have height_resolution_${site}_${RESOLUTION_PARAMS_TAG}.rds -- skipping"
    continue
  fi
  if [ -n "${HEIGHT_JOB_IDS[$site]+x}" ]; then
    dep="$(join_dep "$NICHE_ID" "$PARAMS_ID" ${HEIGHT_JOB_IDS[$site]})"
  else
    dep="$NICHE_PARAMS_DEP"
  fi
  # --time override: run_height_resolution_experiment.sh's own 6h default was
  # calibrated only from Maquipucuna's actual run (3h13m total for all 4
  # resolutions); other sites have noticeably larger landscapes/climate-cache
  # costs (e.g. Mashpi's climate-cache-dominated single-config runs already
  # ran longer), so reusing that budget blindly risks the same class of
  # timeout MindoTarabita's factorial hit on lm_short's 8h cap. lm_medium
  # allows up to 24h; 16h leaves real margin without over-claiming it.
  id=$(sbatch --parsable --dependency="$dep" --time=16:00:00 \
    scripts/02_model/resolution/run_height_resolution_experiment.sh \
    "$site" "0.1,0.25,0.5,1.0" "$RESOLUTION_PARAMS")
  log_job "resolution_${site}" "$id"
  RESOLUTION_IDS+=("$id")
done

echo
echo "== Step 8: reproduction factorial, every viable site (fresh -- old MindoTarabita checkpoint was built against the pre-6-site species list, no longer valid to resume from; moved aside to data/processed/stale_pre_6site_data/) =="
# MindoMirador (renamed MiradorMindo) and Yanayacu excluded: their best_case
# run still crashes with "missing value where TRUE/FALSE needed" (0/5
# replicates, 2026-07-16, still unresolved), and the factorial's upper
# parameter levels approach best_case's values, so queuing 625 combinations
# there would likely just hit the same crash 625x for no benefit. Fix that
# separately before adding these two.
FACTORIAL_IDS=()
for site in "${SITES[@]}"; do
  if [ "$site" = "MindoMirador" ] || [ "$site" = "Yanayacu" ]; then
    echo "  $site: excluded (see note above)"
    continue
  fi
  id=$(sbatch --parsable --dependency="$NICHE_PARAMS_DEP" --partition=lm_medium --time=10:00:00 \
    scripts/02_model/run/run_colonization.sh "$site" data/params/reproduction_factorial_v3.rds \
    reproduction_factorial_v3 0.25)
  log_job "factorial_${site}" "$id"
  FACTORIAL_IDS+=("$id")
done

echo
echo "== Step 9: final plots + summary (depends on everything above) =="
FINAL_DEPS=("$CLIMVAR_ID" "$COMPETITION_ID" "${FACTORIAL_IDS[@]}" "${RESOLUTION_IDS[@]}")
FINAL_DEP="$(join_dep_any "${FINAL_DEPS[@]}")"
PLOT_ID=$(sbatch --parsable --dependency="$FINAL_DEP" scripts/02_model/plots/run_plots.sh)
log_job "plots" "$PLOT_ID"
SUMMARY_ID=$(sbatch --parsable --dependency="$FINAL_DEP" --partition=intelsr_short \
  --account=ag_biob_scabral --time=00:10:00 --ntasks=1 --output=logs/log_%j.out \
  --wrap="module purge; module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a; Rscript scripts/02_model/analysis/summarize_all_results.R > logs/summary_\$SLURM_JOB_ID.txt")
log_job "summary" "$SUMMARY_ID"

ALL_IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")
echo
echo "All steps submitted ($(echo "${ALL_IDS[@]}" | wc -w) jobs) -- job map saved to $MANIFEST"
echo "You can close this terminal now; SLURM will run the chain unattended."
echo "Check progress any time with: squeue -u \$USER"
echo "Or after the fact with:       sacct -j $ALL_IDS_CSV --format=JobID,JobName,State,Elapsed"
