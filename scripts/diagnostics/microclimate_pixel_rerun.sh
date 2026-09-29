#!/bin/bash
# microclimate_pixel_rerun.sh — v7 3D re-run, Phase A: per-pixel-retained
# microclimate output over each site's own landscape footprint (see
# docs/methods_update_report.md, "Four pricing questions" + "Approved:
# three-dimensional microenvironment").
#
# Full ~33km domain kept for horizon/terrain physics; only the write-time
# reduction changes (run_microclimate_site.R / .compute_pixel_means(),
# get_colonization.R) -- means, not quantiles, restricted to each site's own
# footprint.
#
# 2026-09-07 (MERGE RETRACTED): this originally ran MindoMirador+
# MindoTarabita+Saloya+LaElenita as ONE combined-domain job (A3) to save
# compute and give the 4 sites identical forcing. Retracted once that job
# sat queued on vlm (only 5 cluster-wide nodes, occupied by other groups'
# week-long jobs) while its own measured-then-extrapolated memory need
# (~2.0-2.1TB, from the merged domain's real cell count -- see below) turned
# out to exceed lm's ~2048GB/node ceiling. Compute was never the actual
# constraint here -- vlm availability was -- and 4 separate lm jobs start in
# seconds where 1 vlm job waits days. Back to 6 solo-site runs (Maquipucuna,
# Mashpi, Yanayacu, MindoMirador, MindoTarabita, Saloya, LaElenita --
# 7 total; LaElenita excluded from colonization but still gets the
# microclimate re-run for E2-style microclimate analyses). Cost: the 4
# Mindo-cluster sites no longer share one forcing run -- each gets its own
# independent ERA5 pull/point model. Their existing rasters already overlap
# 55-95% and are 7-8km apart, so the climate difference this reintroduces
# is expected to be small, not zero -- flagged, not measured.
#
# Merged-domain memory arithmetic (for the record, since it drove the
# retraction): merged padded bbox (union of the 4 sites' raw obs bboxes,
# +0.15deg each side, same convention as every solo site) = 44.67km x
# 40.52km at 90m cells = 496 x 450 = 223,460 cells, vs each solo site's own
# ~137,640-140,624 cells (~138,758 average) -- a 1.61x ratio, not the
# "modestly larger" the 55-95% raster overlap alone would suggest, because
# overlap describes shared AREA, not the padded bounding box's shape (the 4
# sites are strung out linearly, so their union's bounding box grows faster
# than their shared area). Scaling the measured single-site peak (1169-1305
# GiB at ~138k cells, this project's own sacct/seff data, jobs
# 27063885-27063891) linearly by that 1.61x ratio gives ~1990-2119 GiB --
# over the ~1800GB threshold -- hence the retraction to 4 solo runs.
#
# 2026-09-07 (concurrency): every run below is independent, no dependency
# between them -- the shared-scratch-capacity caution from
# microenv_array.sh's original 7-site chain still applies in principle (see
# earlier git history) but is accepted here given real historical use at
# similar concurrency without incident since that one 2026-08-04 night.
#
# ONLY_LABELS (optional env var, comma-separated): restrict this invocation
# to a subset of the 7 labels below -- e.g. `ONLY_LABELS=MindoMirador,
# MindoTarabita,Saloya,LaElenita sbatch ...` to (re)launch just the 4 that
# were on the retracted merged domain, without resubmitting sites already
# running under an earlier invocation of this script.
#
# Run from: /home/s38leste_hpc/seres/
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --output=/home/s38leste_hpc/seres/logs/log_array_%j.out

# Sized from this project's own measured sacct/seff MaxRSS+Elapsed (see
# docs/methods_update_report.md for the full table): every solo site's peak
# is 1169-1305 GiB (+20% headroom -> 1500-1600G, comfortably inside lm's
# ~2048GB/node ceiling). Time = ~1.5x each site's own historical elapsed,
# rounded up; any site whose 1.5x figure lands within ~1h of a partition's
# time cap is bumped to the next tier rather than requesting right at the
# wall (this project's own history shows the SAME site taking 2-3x longer
# on a bad attempt).
# 2026-09-08 (grid-snap re-run): uniform lm_long/1500G/concurrent across all
# 7 sites, per explicit instruction -- domain WIDTH is unchanged by the
# grid-snap (only shifted), so this project's own measured peak (1169-1305
# GiB, +20% headroom) still applies uniformly; not resizing per-site this
# time for simplicity/consistency now that this is the intended standard
# domain-construction method going forward, not a one-off.
ALL_LABELS=(     "Maquipucuna" "Mashpi"    "Yanayacu"  "MindoMirador" "MindoTarabita" "Saloya"    "LaElenita")
ALL_RUNS=(       "Maquipucuna" "Mashpi"    "Yanayacu"  "MindoMirador" "MindoTarabita" "Saloya"    "LaElenita")
ALL_PARTITIONS=( "lm_long"     "lm_long"   "lm_long"   "lm_long"      "lm_long"       "lm_long"   "lm_long")
ALL_MEM=(        "1500G"       "1500G"     "1500G"     "1500G"        "1500G"         "1500G"     "1500G")
ALL_TIME=(       "30:00:00"    "30:00:00"  "30:00:00"  "30:00:00"     "30:00:00"      "45:00:00"  "36:00:00")

N_MONTHS=${1:-12}
HEIGHT_STEP=${2:-0.1}
VERSION="v7pix"

EXCLUDE_NODES=${EXCLUDE_NODES-vlmnode219}
excl_arg=${EXCLUDE_NODES:+--exclude=${EXCLUDE_NODES}}

IFS=',' read -ra ONLY_ARR <<< "${ONLY_LABELS:-}"

JOB_IDS=()
for i in "${!ALL_LABELS[@]}"; do
  LABEL="${ALL_LABELS[$i]}"
  if [ -n "${ONLY_LABELS:-}" ]; then
    match=0
    for want in "${ONLY_ARR[@]}"; do [ "$want" = "$LABEL" ] && match=1; done
    [ "$match" = 0 ] && continue
  fi
  RUN="${ALL_RUNS[$i]}"
  PART="${ALL_PARTITIONS[$i]}"
  MEM="${ALL_MEM[$i]}"
  TIME="${ALL_TIME[$i]}"

  JOB_SCRIPT=$(mktemp)
  cat > "$JOB_SCRIPT" <<EOS
#!/bin/bash
#SBATCH --job-name=microclimate_pixel_rerun_${LABEL}_v7pix
#SBATCH --partition=${PART}
#SBATCH --account=ag_biob_scabral
#SBATCH --time=${TIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=${MEM}
#SBATCH --requeue
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=FAIL,TIME_LIMIT,TIME_LIMIT_80
#SBATCH --mail-user=s38leste@uni-bonn.de
module purge
module load GCCcore/13.3.0
module load R/4.4.2-gfbf-2024a
module load Miniforge3/24.1.2-0
module load UDUNITS/2.2.28-GCCcore-13.2.0
export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:\$LD_LIBRARY_PATH"
unset PYTHONPATH
export LD_PRELOAD="/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libcrypto.so.3:/home/s38leste_hpc/.conda/envs/canopy_rgee/lib/libssl.so.3"
export CANOPY_PYTHON="/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python"
export CANOPY_OBS_CSV="data/csv/combinedv6.csv"
RUN_SITES=\$1
N_MONTHS=\${2:-12}
HEIGHT_STEP=\${3:-0.1}
export CANOPY_SCRATCH=\$(ws_allocate seres 90)
echo "Scratch workspace: \$CANOPY_SCRATCH"
mkdir -p "\$CANOPY_SCRATCH/tmp"
export TMPDIR="\$CANOPY_SCRATCH/tmp"
cd /home/\$USER/seres
\$CANOPY_PYTHON -c "import ee; print('ee import OK')" 2>&1
Rscript scripts/01_microclimate/run_microclimate_site.R "\$RUN_SITES" "\$N_MONTHS" "\$HEIGHT_STEP"
EOS

  id=$(sbatch --parsable $excl_arg \
    --job-name="microclimate_pixel_rerun_${LABEL}_${VERSION}" \
    "$JOB_SCRIPT" "$RUN" "$N_MONTHS" "$HEIGHT_STEP")
  rm -f "$JOB_SCRIPT"
  echo "Submitted ${LABEL} (sites: ${RUN}) on ${PART} mem=${MEM} time=${TIME}: job ${id} (concurrent, no dependency)"
  JOB_IDS+=("$id")
done

echo "ALL_JOB_IDS=${JOB_IDS[*]}"
# NOTE: this invocation does NOT submit the batch summary job itself when
# ONLY_LABELS restricts the run -- the caller (e.g. a partial re-launch
# alongside already-running jobs from a prior invocation) is responsible for
# submitting one afterany summary covering the FULL real job-ID set, not
# just this invocation's subset. See docs/methods_update_report.md for the
# worked example from the MindoCluster retraction.
if [ -z "${ONLY_LABELS:-}" ]; then
  summary_id=$(sbatch --parsable --dependency=afterany:$(IFS=:; echo "${JOB_IDS[*]}") \
    --job-name="microclimate_pixel_rerun_${VERSION}_summary" \
    --mail-type=END --mail-user=s38leste@uni-bonn.de \
    --output=logs/%x_%j.out --error=logs/%x_%j.err \
    scripts/diagnostics/batch_summary_mailer.sh "microclimate_pixel_rerun_${VERSION}" "${JOB_IDS[@]}")
  echo "Submitted summary job: ${summary_id} (depends on afterany: ${JOB_IDS[*]})"
  echo "SUMMARY_JOB_ID=${summary_id}"
fi
