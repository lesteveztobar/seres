# export_niche_tables.R — CSV exports of the niche inputs and per-timestep
# run summaries behind the Results figures.
#
#   output/tables/niche/niche_curves.csv      species, axis, x, score, n_train
#   output/tables/niche/niche_background.csv  axis, x, density
#   output/tables/niche/niche_by_height.csv   site, species, height,
#                                             relative_height, before, after
#   output/tables/run_summaries/<run file>_summary.csv
#                                             rep, t, totalS, totalJ, totalA,
#                                             total, extinct, recruited, regime
#
# Niche version: everything comes from the cache the Results runs consumed,
# i.e. NICHE_CACHE_PATH as the run launchers set it (species_niches.rds, the
# v8 rebuild of 2026-09-25). Do not point this at species_niches_v6/_v7 --
# the curve and by-height figures must describe the niches the model ran on.
#
# niche_curves: the cache stores one raw 0-100 score per axis; the ceiling
# rescale is applied per SITE to the combined (geometric-mean) score, never
# per axis, so there is no per-axis rescaled curve. before/after in
# niche_by_height is the rescale comparison.
#
# n_train: training records per species after the fixed held-out split
# (get_held_out_split(), same call characterize_niches.R makes). Checked
# against the counts printed by the cache's own build log, if given.
#
# niche_by_height: the cached .niche_score_rows_all_sites_cached() table
# (same one the niche figures use), refused if it is older than any input.
# relative_height = height / site_canopy_ceiling() -- the ceiling the runs
# used for canopy_z (run_colonization.R).
#
# Usage: sbatch scripts/experiments/A17_niche_suitability_plots/run_export_niche_tables.sh

source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")

NICHE_OUT <- file.path(OUTPUT_DIR, "tables", "niche")
SUMMARY_OUT <- file.path(OUTPUT_DIR, "tables", "run_summaries")
dir.create(NICHE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_OUT, recursive = TRUE, showWarnings = FALSE)
BUILD_LOG <- Sys.getenv("CANOPY_NICHE_BUILD_LOG", unset = NA_character_)
HEIGHT_STEP <- 0.4

cat("Observations:", OBSERVATIONS_CSV, "\n")
cat("Niche cache: ", NICHE_CACHE_PATH, "\n")
cat("Background:  ", NICHE_BACKGROUND_PATH, "\n\n")

# ── n_train ──────────────────────────────────────────────────────────────────
obs <- load_observations()
held_out <- get_held_out_split(obs)
n_train <- table(obs$FinalID[!held_out])
cat(sprintf("Training records: %d of %d\n", sum(!held_out), length(held_out)))

if (!is.na(BUILD_LOG)) {
  lines <- grep("^\\S.* +[0-9]+ observations \\(", readLines(BUILD_LOG), value = TRUE)
  logged <- setNames(as.integer(sub("^.*? +([0-9]+) observations \\(.*$", "\\1", lines, perl = TRUE)),
                     trimws(sub(" +[0-9]+ observations \\(.*$", "", lines)))
  mismatch <- names(logged)[logged != as.integer(n_train[names(logged)])]
  if (length(mismatch) > 0 || !setequal(names(logged), names(n_train)))
    stop("n_train does not match the cache build log for: ",
         paste(c(mismatch, setdiff(union(names(logged), names(n_train)),
                                   intersect(names(logged), names(n_train)))), collapse = ", "))
  cat("n_train matches", BUILD_LOG, "for all", length(logged), "species\n")
}

# ── niche_curves.csv ─────────────────────────────────────────────────────────
niche_cache <- readRDS(NICHE_CACHE_PATH)
curves <- do.call(rbind, lapply(names(niche_cache), function(sp) {
  nc <- niche_cache[[sp]]
  if (!identical(names(nc), "axes")) stop(sp, ": unexpected cache fields: ", paste(names(nc), collapse = ", "))
  do.call(rbind, lapply(names(nc$axes), function(ax) {
    a <- nc$axes[[ax]]
    if (!identical(sort(names(a)), c("score", "x"))) stop(sp, "/", ax, ": unexpected axis fields")
    data.frame(species = sp, axis = ax, x = a$x, score = a$score)
  }))
}))
missing_n <- setdiff(unique(curves$species), names(n_train))
if (length(missing_n) > 0) stop("No training records for cached species: ", paste(missing_n, collapse = ", "))
curves$n_train <- as.integer(n_train[curves$species])
write.csv(curves, file.path(NICHE_OUT, "niche_curves.csv"), row.names = FALSE)
cat(sprintf("niche_curves.csv: %d rows, %d species\n", nrow(curves), length(unique(curves$species))))

# ── niche_background.csv ─────────────────────────────────────────────────────
bg <- readRDS(NICHE_BACKGROUND_PATH)
background <- do.call(rbind, lapply(names(bg), function(ax)
  data.frame(axis = ax, x = bg[[ax]]$x, density = bg[[ax]]$y)))
write.csv(background, file.path(NICHE_OUT, "niche_background.csv"), row.names = FALSE)
cat(sprintf("niche_background.csv: %d rows\n", nrow(background)))

# ── niche_by_height.csv ──────────────────────────────────────────────────────
# Rebuilds (slow, ~1 h) only if the cached table is stale -- the cached
# wrapper's own mtime check.
built <- .niche_score_rows_all_sites_cached(height_step = HEIGHT_STEP)
rows <- built$rows
if (length(built$excluded) > 0) cat("Sites with no confirmed species (no rows):", paste(built$excluded, collapse = ", "), "\n")
ceilings <- vapply(unique(rows$site), function(s) {
  microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h%.2f.rds", s, HEIGHT_STEP)))
  site_canopy_ceiling(s, obs, max(microenv_heights(microenv)))
}, numeric(1))
cat("Canopy ceilings (m):", paste(sprintf("%s=%.1f", names(ceilings), ceilings), collapse = ", "), "\n")
rows$relative_height <- rows$height / ceilings[rows$site]
if (any(rows$relative_height > 1 + 1e-9)) warning("Some heights exceed the site canopy ceiling")
by_height <- rows[, c("site", "species", "height", "relative_height", "before", "after")]
write.csv(by_height, file.path(NICHE_OUT, "niche_by_height.csv"), row.names = FALSE)
cat(sprintf("niche_by_height.csv: %d rows\n", nrow(by_height)))

# ── run summaries ────────────────────────────────────────────────────────────
summary_cols <- c("rep", "t", "totalS", "totalJ", "totalA", "total", "extinct", "recruited", "regime")
run_files <- list.files(PROCESSED_DIR, full.names = TRUE, pattern =
  "^colonization_.*_(realistic|best_case|isolation_.*|horizontal_resolution_(pooled|voxel))_h0\\.40\\.rds$")
cat(sprintf("\nExporting $summary from %d run files\n", length(run_files)))
for (f in run_files) {
  r <- readRDS(f)
  s <- r$summary
  rm(r); invisible(gc())
  if (is.null(s) || !all(summary_cols %in% names(s))) {
    message("SKIPPED -- ", basename(f), ": no $summary with the expected columns")
    next
  }
  out <- file.path(SUMMARY_OUT, sub("\\.rds$", "_summary.csv", basename(f)))
  write.csv(s[, summary_cols], out, row.names = FALSE)
  cat(sprintf("  %s: %d rows\n", basename(out), nrow(s)))
}
cat("\nDone.\n")
