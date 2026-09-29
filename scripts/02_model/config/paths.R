# paths.R
# Project directory constants for seres
# All paths derived from BASE_DIR — change only BASE_DIR if project moves
# Per-site subdirectories (era5, dtm, soil, etc.) are built dynamically
# inside the site loop in getmicroenv.R using BASE_DIR as the root
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

BASE_DIR <- "/home/s38leste_hpc/seres"

# data
RAW_DIR <- file.path(BASE_DIR, "data", "raw")
CSV_DIR <- file.path(BASE_DIR, "data", "csv")
PROCESSED_DIR <- file.path(BASE_DIR, "data", "processed")
PARAMS_DIR <- file.path(BASE_DIR, "data", "params")

# Master field-observation CSV, read by every script that needs
# species/site/height data (niche characterization, colonization runs,
# microenv site bounding boxes, diagnostics, ...). Overridable via
# CANOPY_OBS_CSV so a SLURM job or shell script can point a whole run at a
# different snapshot (e.g. combinedv3.csv, or a future post-migration file)
# without editing R source. Default switched from combinedv3.csv to
# combined_with_identification.csv on 2026-07-24, once the latter had all
# seven sites (post rebuild_combined_csv.py LaElenita/MindoMirador/Saloya
# split -- see scripts/data_prep/rebuild_combined_csv.py) and combinedv3.csv did not.
#
# 2026-08-28: briefly pointed straight at combinedv6.csv, then reverted --
# Saloya's still-pending OAT-chain jobs (27169761-27169766) each freshly
# source() this file when they actually start (paths.R is re-read from disk
# at job START, not at sbatch submission), so changing the *default* here
# would have silently switched them mid-sweep to a different dataset than
# whatever the already-running/already-completed steps of that same chain
# used -- an inconsistent, uninterpretable mix. combinedv6.csv work uses the
# CANOPY_OBS_CSV override explicitly instead (see run_colonization_v6.sh /
# any v6-tagged submission), leaving this default -- and every currently
# in-flight job -- untouched.
# 2026-09-06 (v7 rebuild, CORRECTED): a deprecation COMMENT does not stop a
# script that forgets to set CANOPY_OBS_CSV from silently resolving to the
# old, pre-v7-schema file (combined_with_identification.csv -- fewer
# columns, an explicit FinalID column rather than one derived from
# Identification, missing every v7 fix applied to combinedv6.csv). A
# comment is not a control; this is not a silent behaviour change, it is
# the removal of one that was already silently possible. No unset default
# at all now -- any script that reaches this line without CANOPY_OBS_CSV
# set stops immediately with a message naming the file to use.
OBSERVATIONS_CSV <- Sys.getenv("CANOPY_OBS_CSV", unset = NA_character_)
if (is.na(OBSERVATIONS_CSV)) {
  stop("CANOPY_OBS_CSV is not set. This pipeline no longer has a silent ",
       "default observations CSV (removed 2026-09-06 -- see paths.R). ",
       "Set CANOPY_OBS_CSV=data/csv/combinedv6.csv (the current v7 dataset) ",
       "explicitly before sourcing paths.R.")
}

# Loads OBSERVATIONS_CSV and populates FinalID from Identification
# (2026-07-24: combined_with_identification.csv leaves FinalID entirely
# unpopulated on every row -- species IDs live in Identification instead,
# pending manual review/promotion into FinalID, same as every prior
# field-data drop; without this, characterize_niches.R and everything
# downstream finds zero usable observations at any site, not just the
# newly added ones -- see the 2026-07-24 characterize_niches.sh failure).
# Rows where Identification is itself blank/NA (e.g. LaElenita: no photo,
# or photo pending review) get a hardcoded default of "Maxillaria
# acutifolia" rather than being left blank, which otherwise silently
# produces a species with an empty-string ID downstream. This default is
# applied only in memory here, never written back to the CSV. NOTE: this
# feeds Identification's guesses -- including ones Confidence marks
# "unverified"/AI-only -- directly into species-level niche and
# colonization modeling. If combinedv3.csv-quality curation matters more
# than the extra sites, override with CANOPY_OBS_CSV=data/csv/combinedv3.csv
# instead of relying on this fallback.
load_observations <- function(path = OBSERVATIONS_CSV) {
  niches <- read.csv(path)
  # 2026-08-29/30: Maxillariinae scope filter (.filter_maxillariinae(),
  # shared_helpers.R) -- renames synonym genera to their accepted name, then
  # drops any row whose genus isn't Maxillariinae, AND (2026-08-30) drops
  # any row with no Identification at all. Must run first, directly on the
  # raw Identification column, before FinalID is even derived from it --
  # every downstream consumer (characterize_niches.R, every site/species
  # helper added this session) needs to see already-scoped, already-renamed
  # data with no fabricated placeholder species in it.
  #
  # The blank-ID-defaults-to-"Maxillaria acutifolia" fallback that used to
  # live here is REMOVED (2026-08-30, explicit decision) -- it fabricated a
  # specific species identity for an individual nobody actually identified,
  # which is exactly what let unidentified-heavy sites (LaElenita/
  # MindoMirador/Saloya) silently borrow OTHER sites' real acutifolia niche
  # evidence in earlier suitability figures. Every row reaching this point
  # now has a genuine, in-scope Identification -- FinalID is just a direct
  # copy of it, no defaulting left to do.
  if ("Identification" %in% names(niches)) {
    niches <- .filter_maxillariinae(niches)
    niches$FinalID <- niches$Identification
  }
  niches
}

# Cross-site species-niche cache (characterize_niches.R's output, loaded by
# init_colonization() (get_colonization.R), .niche_plot_context() (plot_
# functions.R), check_niche_suitability.R, and summarize_all_results.R).
# Same override pattern as OBSERVATIONS_CSV, added 2026-08-28 for the same
# reason: characterize_niches.R must be re-run against combinedv6.csv, but
# overwriting species_niches.rds/niche_background_density.rds in place would
# silently swap the niche model out from under every currently-running or
# still-pending SLURM job that reads it fresh at its own start (same hazard
# as OBSERVATIONS_CSV -- see that constant's comment above). Left unset,
# both fall back to today's exact literal paths, so nothing in flight is
# affected; v6 work sets CANOPY_NICHE_CACHE/CANOPY_NICHE_BACKGROUND
# explicitly (see characterize_niches.R's invocation for v6).
NICHE_CACHE_PATH <- Sys.getenv("CANOPY_NICHE_CACHE",
  unset = file.path(PROCESSED_DIR, "species_niches.rds")
)
NICHE_BACKGROUND_PATH <- Sys.getenv("CANOPY_NICHE_BACKGROUND",
  unset = file.path(PROCESSED_DIR, "niche_background_density.rds")
)

# project
SCRIPTS_DIR <- file.path(BASE_DIR, "scripts")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
LOGS_DIR <- file.path(BASE_DIR, "logs")

# 02_model subdirectories — mirrors the current scripts/02_model/ layout
# (analysis, config, diagnostics, engine, plots, resolution, run, setup,
# tests). Update here if that layout changes rather than hardcoding
# sub-paths at each call site.
MODEL_DIR <- file.path(SCRIPTS_DIR, "02_model")
MODEL_ANALYSIS_DIR <- file.path(MODEL_DIR, "analysis")
MODEL_CONFIG_DIR <- file.path(MODEL_DIR, "config")
MODEL_DIAGNOSTICS_DIR <- file.path(MODEL_DIR, "diagnostics")
MODEL_ENGINE_DIR <- file.path(MODEL_DIR, "engine")
MODEL_PLOTS_DIR <- file.path(MODEL_DIR, "plots")
MODEL_RESOLUTION_DIR <- file.path(MODEL_DIR, "resolution")
MODEL_RUN_DIR <- file.path(MODEL_DIR, "run")
MODEL_SETUP_DIR <- file.path(MODEL_DIR, "setup")
MODEL_TESTS_DIR <- file.path(MODEL_DIR, "tests")

# create all directories if they don't exist
# showWarnings = FALSE silently skips dirs that already exist
# NOTE: only data/output/log dirs are created here -- the scripts/02_model/*
# subdirectories above are source-code directories checked into the repo,
# not runtime output, so they are intentionally excluded from auto-creation.
for (d in c(RAW_DIR, CSV_DIR, PROCESSED_DIR, PARAMS_DIR, SCRIPTS_DIR, OUTPUT_DIR, LOGS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Analysis-ID-aware output path helper (added 2026-09-29, repo reorg pass).
# Opt-in only: with CANOPY_OUTPUT_ROOT unset, this is a strict no-op that
# returns `legacy_path` unchanged -- every caller's current hardcoded path
# is untouched unless a run explicitly opts in. Setting CANOPY_OUTPUT_ROOT
# routes output to <root>/by_analysis/<analysis_id>/<state>/<filename>
# instead, creating the directory on demand and never silently overwriting
# an existing file there (a numeric suffix is appended and a warning
# printed, rather than clobbering it -- same hazard this file's
# NICHE_CACHE_PATH/OBSERVATIONS_CSV comments already describe for in-flight
# jobs reading a path fresh at their own start).
resolve_output_path <- function(analysis_id, filename, legacy_path, state = "current") {
  output_root <- Sys.getenv("CANOPY_OUTPUT_ROOT", unset = NA_character_)
  if (is.na(output_root)) {
    return(legacy_path)
  }
  dest_dir <- file.path(output_root, "by_analysis", analysis_id, state)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest_path <- file.path(dest_dir, filename)
  if (!file.exists(dest_path)) {
    return(dest_path)
  }
  ext <- tools::file_ext(filename)
  base <- tools::file_path_sans_ext(filename)
  suffix <- 1
  repeat {
    candidate <- if (nzchar(ext)) {
      file.path(dest_dir, sprintf("%s_%d.%s", base, suffix, ext))
    } else {
      file.path(dest_dir, sprintf("%s_%d", base, suffix))
    }
    if (!file.exists(candidate)) {
      warning(sprintf(
        "resolve_output_path(): %s already exists, writing to %s instead",
        dest_path, candidate
      ))
      return(candidate)
    }
    suffix <- suffix + 1
  }
}
