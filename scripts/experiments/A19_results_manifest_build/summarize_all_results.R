# summarize_all_results.R
# Consolidated text summary across every site and experiment that's actually
# been run so far: niche characterization coverage, persistence validation
# (best_case/realistic), and reproduction factorial (v3) outcomes. Same
# per-result logic as check_colonization_run.R (kept in sync deliberately --
# same aggregate formula, same wording), just looped over every site instead
# of one site/tag pair at a time, so you get one consolidated report instead
# of running it 10+ times by hand.
#
# Read-only, no cluster job needed -- only needs paths.R (doesn't even source
# get_colonization.R, since it only reads already-saved RDS/CSV results, not
# raw climate), so no module-loading wrapper is required.
#
# Usage: Rscript scripts/experiments/A19_results_manifest_build/summarize_all_results.R
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
# Sites derived from OBSERVATIONS_CSV itself, not hardcoded -- a newly added
# site is picked up automatically the next time this runs.
SITES <- sort(unique(niches$Area_or_Site))

# ── Niche characterization coverage ─────────────────────────────────────────
cat("========================================\n")
cat("Niche characterization\n")
cat("========================================\n")
niche_cache_path <- NICHE_CACHE_PATH
if (file.exists(niche_cache_path)) {
  niche_cache <- readRDS(niche_cache_path)
  n_cached <- sum(!vapply(niche_cache, is.null, logical(1)))
  cat(sprintf("species_niches.rds: %d/%d species have a usable niche (cross-site cache)\n",
              n_cached, length(niche_cache)))

  for (site in SITES) {
    site_species <- sort(unique(niches$FinalID[niches$Area_or_Site == site]))
    if (length(site_species) == 0) next
    has_niche <- sum(!vapply(niche_cache[site_species], is.null, logical(1)))
    cat(sprintf("  %-15s %d species observed, %d with a cached niche\n",
                site, length(site_species), has_niche))
  }
} else {
  cat("No species_niches.rds yet -- run characterize_niches.R first.\n")
}

# ── One result's worth of summary (sweep/factorial data.frame, or a
# run_replicated() list) -- same logic as check_colonization_run.R, factored
# into a function so it can be looped instead of run once per invocation. ──
summarize_one <- function(in_path, label) {
  if (!file.exists(in_path)) return(invisible(FALSE))
  result <- readRDS(in_path)
  cat(sprintf("\n-- %s --\n", label))

  shape <- .classify_result_shape(result)  # shared_helpers.R

  if (shape$shape == "sweep") {
    if (shape$degenerate) {
      cat("  (empty/degenerate result)\n")
      return(invisible(TRUE))
    }
    combo <- shape$combo
    cat(sprintf("  Swept: %s | %d combination(s) | final year t=%d\n",
                paste(shape$swept, collapse = " x "), nrow(combo), shape$t_max))
    cat(sprintf("  Persisted: %d/%d (%.1f%%)\n",
                sum(combo$persisted), nrow(combo), 100 * mean(combo$persisted)))

    for (v in shape$swept) {
      rate <- tapply(combo$persisted, combo[[v]], mean)
      cat(sprintf("    Persistence by %s: %s\n", v,
                  paste(sprintf("%s=%.1f%%", names(rate), 100 * rate), collapse = ", ")))
    }

    persisting <- shape$persisting
    if (nrow(persisting) > 0) {
      cat(sprintf("  Among persisting: final total abundance %.0f-%.0f (median %.0f)\n",
                  min(persisting$total), max(persisting$total), median(persisting$total)))
      for (v in shape$swept) {
        cat(sprintf("    %s bounding box among persisting: %s to %s\n",
                    v, format(min(persisting[[v]])), format(max(persisting[[v]]))))
      }
    } else {
      cat("  No combinations persisted.\n")
    }
    return(invisible(TRUE))
  }

  # run_replicated() list shape: list(runs = <one per replicate>, summary = <tidy df>)
  cat(sprintf("  %d/%d replicate(s) succeeded\n", shape$n_ok, shape$n_runs))
  if (shape$n_ok == 0 || is.null(shape$summary)) return(invisible(TRUE))

  final <- shape$final
  cat(sprintf("  Persisted: %d/%d replicates (adults present in 2nd half of run)\n",
              sum(!final$extinct), nrow(final)))
  # n_recruiting (2026-08-28, .classify_result_shape()): distinguishes
  # genuine self-sustaining persistence from a founder cohort just decaying
  # with zero seedling/juvenile replacement -- see get_colonization.R's
  # `recruited` column.
  if (!is.null(shape$n_recruiting) && !is.na(shape$n_recruiting)) {
    cat(sprintf("  Recruiting: %d/%d replicates (any S/J in 2nd half of run)\n",
                shape$n_recruiting, nrow(final)))
  }
  cat(sprintf("  Final-year (t=%d) total abundance by replicate: %s\n",
              shape$final_t, paste(round(final$total), collapse = ", ")))
  invisible(TRUE)
}

cat("\n\n========================================\n")
cat("Persistence validation (best_case / realistic / realistic_273founders)\n")
cat("========================================\n")
any_persistence <- FALSE
for (site in SITES) {
  ok1 <- summarize_one(file.path(PROCESSED_DIR, sprintf("colonization_%s_best_case_h0.40.rds", site)),
                       sprintf("%s / best_case", site))
  ok2 <- summarize_one(file.path(PROCESSED_DIR, sprintf("colonization_%s_realistic_h0.40.rds", site)),
                       sprintf("%s / realistic", site))
  ok3 <- summarize_one(file.path(PROCESSED_DIR, sprintf("colonization_%s_realistic_273founders_h0.40.rds", site)),
                       sprintf("%s / realistic_273founders", site))
  any_persistence <- any_persistence || isTRUE(ok1) || isTRUE(ok2) || isTRUE(ok3)
}
if (!any_persistence) cat("(no persistence validation results yet)\n")

cat("\n\n========================================\n")
cat("Reproduction factorial (v3)\n")
cat("========================================\n")
any_factorial <- FALSE
for (site in SITES) {
  final_path <- file.path(PROCESSED_DIR,
    sprintf("colonization_%s_reproduction_factorial_v3_h0.40.rds", site))
  ckpt_path <- file.path(PROCESSED_DIR,
    sprintf("colonization_%s_reproduction_factorial_v3_h0.40_checkpoint.rds", site))

  if (file.exists(final_path)) {
    ok <- summarize_one(final_path, sprintf("%s / reproduction_factorial_v3", site))
    any_factorial <- any_factorial || isTRUE(ok)
  } else if (file.exists(ckpt_path)) {
    summarize_one(ckpt_path, sprintf("%s / reproduction_factorial_v3 (INCOMPLETE -- checkpoint only, rerun to resume)", site))
    any_factorial <- TRUE
  }
}
if (!any_factorial) cat("(no factorial results yet)\n")

cat("\nDone.\n")
