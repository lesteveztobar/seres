# mashpi_langcanopy_taskc.R
# Task 3(c) recomputed for Mashpi under the Lang et al. site-specific canopy
# height rerun (Task 6, CANOPY_SITE_TREE_HEIGHT=1), using the EXACT same
# method as held_out_validation.R's Part (c) (held-out obs_val heights vs.
# realized S/J/A voxel heights at the final timestep) -- run against both
# the original best_combo_v6 result and the langcanopy rerun, for direct
# comparison.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

microenv <- readRDS(file.path(PROCESSED_DIR, "microenv_Mashpi_h0.40.rds"))
heights <- microenv_heights(microenv)

for (tag in c("best_combo_v6", "best_combo_v6_langcanopy")) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_Mashpi_%s_h0.40.rds", tag))
  if (!file.exists(f)) { message(tag, ": no result"); next }
  result <- readRDS(f)

  obs_rows <- list()
  realized_rows <- list()
  for (rep_i in seq_along(result$runs)) {
    run <- result$runs[[rep_i]]
    if (is.null(run)) next
    if (!is.null(run$obs_val) && nrow(run$obs_val) > 0) {
      obs_rows[[length(obs_rows) + 1]] <- data.frame(rep = rep_i, observed_height = run$obs_val$Height_m)
    }
    final_t <- dim(run$abundanceA)[4]
    for (stage_arr in list(A = run$abundanceA, J = run$abundanceJ, S = run$abundanceS)) {
      occ <- which(stage_arr[, , , final_t, , drop = FALSE] > 0, arr.ind = TRUE)
      if (nrow(occ) == 0) next
      realized_rows[[length(realized_rows) + 1]] <- data.frame(rep = rep_i, realized_height = heights[occ[, 3]])
    }
  }
  obs_h <- if (length(obs_rows) > 0) do.call(rbind, obs_rows)$observed_height else numeric(0)
  mod_h <- if (length(realized_rows) > 0) do.call(rbind, realized_rows)$realized_height else numeric(0)

  cat(sprintf("\n=== %s ===\nn_observed=%d n_realized=%d\n", tag, length(obs_h), length(mod_h)))
  if (length(obs_h) < 4 || length(mod_h) < 4) { cat("n too small for KS test\n"); next }
  ks <- suppressWarnings(ks.test(obs_h, mod_h))
  cat(sprintf("KS D=%.3f p=%.4f | mean_shift(model-observed)=%.2fm | obs_mean=%.2f mod_mean=%.2f\n",
              ks$statistic, ks$p.value, mean(mod_h) - mean(obs_h), mean(obs_h), mean(mod_h)))
}
