# held_out_validation.R -- held-out validation (post-processing only; no
# colonization runs of its own). Reads an already-completed `best_combo`
# result per site and compares each site's held-out observations against
# it. Depends on `colonization_<site>_best_combo_h0.40.rds` existing --
# which in turn depends on the reproduction-survival factorial (to know
# which combination is "best") -- so this cannot run until that chain
# completes. Rewritten 2026-09-08 per the following (Phase E):
#   - modelled heights are ABUNDANCE-weighted, not occupancy-weighted (a
#     voxel with 8 individuals must count 8x, not 1x, in the modelled
#     height distribution);
#   - the primary comparison is computed PER REPLICATE, never pooled
#     across replicates first -- pooling before comparing understates
#     between-replicate variance;
#   - primary statistic: mean shift (model - observed), with a bootstrap
#     CI (>=10,000 resamples); KS is secondary, reported alongside;
#   - n_observed and n_modelled_individuals are counted as INDIVIDUALS
#     (sum of abundance), not voxels/rows.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

# Input run tag (2026-09-29): the best_combo chain this script used to read
# is void (pre-K-recalibration factorial); default is the literature-parameter
# `realistic` run, which carries the same held-out split (obs_val).
args <- commandArgs(trailingOnly = TRUE)
RUN_TAG <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "realistic"
N_BOOT <- 10000L
SITES <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")

all_rows <- list()
height_compare_rows <- list()   # one row per (site, rep, individual) -- abundance-expanded

for (site in SITES) {
  # Naming per this session's "drop E-numbering, no v-suffix" convention --
  # matches whatever run_colonization.R's own exp_tag="best_combo" produces.
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s_h0.40.rds", site, RUN_TAG))
  if (!file.exists(f)) { message("Skipping ", site, " -- no ", f); next }
  result <- readRDS(f)
  microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site)))
  cc_voxel <- build_clim_cache_voxel(microenv)
  niche_cache <- readRDS(NICHE_CACHE_PATH)

  for (rep_i in seq_along(result$runs)) {
    run <- result$runs[[rep_i]]
    if (is.null(run) || is.null(run$obs_val) || nrow(run$obs_val) == 0) next
    val <- run$obs_val

    # (a)/(b): suitability scores under 3 characterizations, unchanged from
    # before -- these are per-observation, not a place abundance-weighting
    # applies (a held-out INDIVIDUAL's own suitability doesn't get more
    # weight just because the model happens to be crowded elsewhere).
    clim_own <- voxel_climate_table(cc_voxel, microenv,
      data.frame(lon = val$lon, lat = val$lat, height = val$Height_m), NICHE_VARS, months = "annual")
    clim_height_only <- voxel_climate_table(cc_voxel, microenv,
      data.frame(lon = NA_real_, lat = NA_real_, height = val$Height_m), NICHE_VARS, months = "annual")
    site_mean <- colMeans(clim_own, na.rm = TRUE)

    for (i in seq_len(nrow(val))) {
      sp <- val$FinalID[i]
      niche_sp <- niche_cache[[sp]]
      if (is.null(niche_sp)) next
      score_own    <- niche_raw_score(as.list(clim_own[i, ]), niche_sp)
      score_height <- niche_raw_score(as.list(clim_height_only[i, ]), niche_sp)
      score_site   <- niche_raw_score(as.list(site_mean), niche_sp)
      all_rows[[length(all_rows) + 1]] <- data.frame(
        site = site, rep = rep_i, species = sp, observed_height = val$Height_m[i],
        score_own_position = score_own, score_height_only = score_height, score_site_mean = score_site
      )
    }

    # (c) REWRITTEN: abundance-weighted, per-replicate. Every occupied
    # voxel's height is repeated once per individual actually there (not
    # once per voxel), across S+J+A stages, at the final timestep of THIS
    # replicate specifically -- height_compare_rows keeps `rep` as its own
    # column precisely so nothing downstream can accidentally pool reps
    # before computing the primary statistic.
    heights <- microenv_heights(microenv)
    final_t <- dim(run$abundanceA)[4]
    rep_modelled_heights <- numeric(0)
    for (stage_arr in list(A = run$abundanceA, J = run$abundanceJ, S = run$abundanceS)) {
      slice <- stage_arr[, , , final_t, , drop = FALSE]
      occ <- which(slice > 0, arr.ind = TRUE)
      if (nrow(occ) == 0) next
      abund <- slice[occ]                      # count of individuals at each occupied voxel
      h_for_voxel <- heights[occ[, 3]]
      # abundance-weighted expansion: this voxel's height appears `abund`
      # times, once per individual, not once per voxel.
      rep_modelled_heights <- c(rep_modelled_heights, rep(h_for_voxel, times = abund))
    }
    if (length(rep_modelled_heights) > 0) {
      height_compare_rows[[length(height_compare_rows) + 1]] <- data.frame(
        site = site, rep = rep_i, modelled_height = rep_modelled_heights)
    }
  }
}

val_df <- do.call(rbind, all_rows)
cat("=== Part (a)/(b): held-out suitability scores, 3 predictor characterizations ===\n")
cat("Pooled sample size (held-out individual x replicate rows):", nrow(val_df), "\n")
print(table(val_df$site, val_df$species))
write.csv(val_df, file.path(OUTPUT_DIR, sprintf("held_out_validation_scores_%s.csv", RUN_TAG)), row.names = FALSE)

cat("\n=== Part (b): paired comparison of suitability under each characterization ===\n")
cat("(Wilcoxon signed-rank, paired by held-out individual x rep -- non-parametric given small/unequal n)\n")
per_sp <- split(val_df, val_df$species)
for (sp in names(per_sp)) {
  d <- per_sp[[sp]]
  cat(sprintf("\n-- %s (n=%d) --\n", sp, nrow(d)))
  if (nrow(d) < 6) {
    cat("  n<6 -- too small for a species-level test, reporting raw scores only:\n")
    print(d[, c("site", "observed_height", "score_own_position", "score_height_only", "score_site_mean")], row.names = FALSE)
    next
  }
  w1 <- tryCatch(wilcox.test(d$score_own_position, d$score_height_only, paired = TRUE), error = function(e) NULL)
  w2 <- tryCatch(wilcox.test(d$score_own_position, d$score_site_mean, paired = TRUE), error = function(e) NULL)
  cat("  own-position vs height-only:", if (!is.null(w1)) sprintf("V=%.1f p=%.4f", w1$statistic, w1$p.value) else "NA", "\n")
  cat("  own-position vs site-mean:  ", if (!is.null(w2)) sprintf("V=%.1f p=%.4f", w2$statistic, w2$p.value) else "NA", "\n")
}

cat("\n\n=== Part (c) REWRITTEN: abundance-weighted, per-replicate mean shift + bootstrap CI ===\n")
realized_df <- if (length(height_compare_rows) > 0) do.call(rbind, height_compare_rows) else data.frame()

boot_mean_shift <- function(obs_h, mod_h, n_boot = N_BOOT) {
  # Bootstrap the mean shift (model - observed) by resampling each side
  # independently with replacement -- standard percentile bootstrap for a
  # difference of means from two independent samples.
  obs_boot <- vapply(seq_len(n_boot), function(i) mean(sample(obs_h, length(obs_h), replace = TRUE)), numeric(1))
  mod_boot <- vapply(seq_len(n_boot), function(i) mean(sample(mod_h, length(mod_h), replace = TRUE)), numeric(1))
  diffs <- mod_boot - obs_boot
  list(mean_shift = mean(mod_h) - mean(obs_h),
       ci_lo = quantile(diffs, 0.025, names = FALSE),
       ci_hi = quantile(diffs, 0.975, names = FALSE))
}

per_rep_rows <- list()
for (site in unique(val_df$site)) {
  # val_df has one row per held-out individual PER REPLICATE; take one
  # replicate only so each observed individual is counted once (2026-09-30:
  # pooling all replicates inflated n_observed by n_reps and narrowed the CI).
  first_rep <- min(val_df$rep[val_df$site == site])
  obs_h <- val_df$observed_height[val_df$site == site & val_df$rep == first_rep]
  n_observed <- length(obs_h)   # already one row per individual (no voxel pooling on the observed side)
  if (n_observed == 0) next
  site_reps <- if (nrow(realized_df) > 0) unique(realized_df$rep[realized_df$site == site]) else integer(0)
  cat(sprintf("\n-- %s: n_observed=%d individuals, %d replicate(s) with modelled data --\n",
              site, n_observed, length(site_reps)))
  for (rep_i in site_reps) {
    mod_h <- realized_df$modelled_height[realized_df$site == site & realized_df$rep == rep_i]
    n_modelled_individuals <- length(mod_h)  # abundance-expanded above -- this IS individuals, not voxels
    if (n_observed < 4 || n_modelled_individuals < 4) {
      cat(sprintf("  rep %d: n_modelled_individuals=%d -- too small for KS/bootstrap, skipping\n", rep_i, n_modelled_individuals))
      next
    }
    bs <- boot_mean_shift(obs_h, mod_h)
    ks <- suppressWarnings(ks.test(obs_h, mod_h))
    cat(sprintf("  rep %d: n_modelled_individuals=%d | mean_shift(model-observed)=%.2fm [%.2f, %.2f] (95%% bootstrap CI, %d resamples) | KS D=%.3f p=%.4f (secondary)\n",
                rep_i, n_modelled_individuals, bs$mean_shift, bs$ci_lo, bs$ci_hi, N_BOOT, ks$statistic, ks$p.value))
    per_rep_rows[[length(per_rep_rows) + 1]] <- data.frame(
      site = site, rep = rep_i, n_observed = n_observed, n_modelled_individuals = n_modelled_individuals,
      mean_shift = bs$mean_shift, ci_lo = bs$ci_lo, ci_hi = bs$ci_hi,
      ks_D = unname(ks$statistic), ks_p = ks$p.value)
  }
}
per_rep_df <- if (length(per_rep_rows) > 0) do.call(rbind, per_rep_rows) else data.frame()
write.csv(per_rep_df, file.path(OUTPUT_DIR, sprintf("held_out_realized_heights_per_rep_%s.csv", RUN_TAG)), row.names = FALSE)
write.csv(realized_df, file.path(OUTPUT_DIR, sprintf("held_out_realized_heights_%s.csv", RUN_TAG)), row.names = FALSE)
cat("\nDone.\n")
