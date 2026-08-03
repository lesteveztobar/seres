# paths.R
# Project directory constants for canopymicroenv
# All paths derived from BASE_DIR — change only BASE_DIR if project moves
# Per-site subdirectories (era5, dtm, soil, etc.) are built dynamically
# inside the site loop in getmicroenv.R using BASE_DIR as the root
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

BASE_DIR <- "/home/s38leste_hpc/canopymicroenv"

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
# split -- see scripts/00_data_conversion/rebuild_combined_csv.py) and combinedv3.csv did not.
OBSERVATIONS_CSV <- Sys.getenv("CANOPY_OBS_CSV",
  unset = file.path(CSV_DIR, "combined_with_identification.csv")
)

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
  if ("Identification" %in% names(niches)) {
    niches$FinalID <- niches$Identification
  }
  needs_default <- is.na(niches$FinalID) | !nzchar(trimws(niches$FinalID))
  niches$FinalID[needs_default] <- "Maxillaria acutifolia"
  niches
}

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
