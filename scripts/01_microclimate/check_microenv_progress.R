# check_microenv_progress.R
# Reports, per site (and optionally per height-tier spacing), how many of the
# expected per-height climate files have been computed on scratch, and
# whether the final microenv_<site>[_h<step>].rds manifest has been saved
# yet. Read-only — safe to run any time.
#
# Usage:
#   export CANOPY_SCRATCH=$(ws_find seres)
#   module load GCCcore/13.3.0 R/4.4.2-gfbf-2024a
#   # Needed so terra can actually load vhgt.tif for sites without measured
#   # CanopyHeight_m — without this, height_ceiling() (lib.R) silently falls
#   # back to hObs_max and under-reports the expected height count for those sites
#   # specifically (the ones WITH measured CanopyHeight_m are unaffected,
#   # since they never touch terra in this script at all).
#   export LD_LIBRARY_PATH="/opt/software/easybuild-INTEL/software/PROJ/9.3.1-GCCcore-13.2.0/lib:/opt/software/easybuild-INTEL/software/GDAL/3.9.0-foss-2023b/lib:/opt/software/easybuild-INTEL/software/GEOS/3.12.1-GCC-13.2.0/lib:$LD_LIBRARY_PATH"
#
#   Rscript scripts/01_microclimate/check_microenv_progress.R                                  # all sites, 0.1m
#   Rscript scripts/01_microclimate/check_microenv_progress.R Maquipucuna 0.1,0.25,0.5,1.0     # one site, several height steps
#   Rscript scripts/01_microclimate/check_microenv_progress.R Maquipucuna,Mashpi,Yanayacu 0.25 # several sites (comma-separated), one step
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")  # .filter_maxillariinae() -- 2026-08-29
source("scripts/01_microclimate/lib.R")

args <- commandArgs(trailingOnly = TRUE)
HEIGHT_STEPS <- if (length(args) >= 2) as.numeric(strsplit(args[2], ",")[[1]]) else 0.1

scratch <- Sys.getenv("CANOPY_SCRATCH", unset = "")
if (nchar(scratch) == 0 || !dir.exists(scratch)) {
  stop("CANOPY_SCRATCH not set to an existing directory.\n",
       "Run: export CANOPY_SCRATCH=$(ws_find seres)")
}

niches <- load_observations()
niches <- niches[!is.na(niches$Height_m), ]

# Default site list derived from OBSERVATIONS_CSV itself, not hardcoded --
# a newly added site is picked up automatically. Explicit sites arg still
# overrides this.
SITES <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else
  sort(unique(niches$Area_or_Site[!is.na(niches$Area_or_Site) & nzchar(niches$Area_or_Site)]))

# Height ceiling now comes from measured CanopyHeight_m first, then the
# canopy-height raster (vhgt.tif), then hObs_max -- see height_ceiling()
# (lib.R), shared with run_microclimate_site.R so "expected" counts match
# reality.

cat(sprintf("%-15s %-6s %-14s %-10s %s\n", "Site", "Step", "Heights", "Manifest", "Height dir"))
cat(strrep("-", 80), "\n")

for (site_name in SITES) {
  hmax <- suppressWarnings(max(niches$Height_m[niches$Area_or_Site == site_name], na.rm = TRUE))
  canopy_max <- suppressWarnings(max(niches$CanopyHeight_m[niches$Area_or_Site == site_name], na.rm = TRUE))
  for (step in HEIGHT_STEPS) {
    suffix <- if (step != 0.1) sprintf("_h%.2f", step) else ""
    height_dir <- file.path(scratch, sprintf("microenv_%s%s_heights", site_name, suffix))
    manifest_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, suffix))
    manifest_status <- if (file.exists(manifest_path)) "done" else "pending"

    if (!is.finite(hmax)) {
      cat(sprintf("%-15s %-6.2f %-14s %-10s %s\n", site_name, step, "no obs in CSV", manifest_status, height_dir))
      next
    }
    ceiling <- height_ceiling(canopy_max, file.path(RAW_DIR, site_name, "vhgt.tif"), hmax)
    n_expected <- length(seq(0.1, ceiling, by = step))

    if (!dir.exists(height_dir)) {
      cat(sprintf("%-15s %-6.2f %-14s %-10s %s\n", site_name, step, sprintf("0/%d", n_expected), manifest_status,
                  paste0(height_dir, " (not created yet)")))
      next
    }
    n_done <- length(list.files(height_dir, pattern = "\\.rds$"))
    cat(sprintf("%-15s %-6.2f %-14s %-10s %s\n", site_name, step, sprintf("%d/%d", n_done, n_expected),
                manifest_status, height_dir))
  }
}
