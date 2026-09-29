# b4_equivalence_test.R -- Phase B4, the equivalence-test GATE. TWO arms,
# both must pass (2026-09-07: added arm 2 after the author flagged that arm
# 1 alone, in `pooled` mode, cannot detect a change confined to the
# per-pixel path -- which is exactly where Phase A's new microclimate
# output lands, and exactly the path every production consumer actually
# uses (stochastic=FALSE, the runcolonization() default never overridden
# by run_colonization.R/run_replicated()/run_experiment())).
#
# Arm 1 ("pooled"): old vs new microenv input, CANOPY_CLIM_MODE=pooled.
# .clim_voxel_slice() returns NA unconditionally in this mode, so both runs
# collapse to the identical flat/pooled climate math regardless of which
# microenv they read from -- confirms the new write-time reduction didn't
# corrupt the pooled *_mean fields every consumer's fallback path reads.
#
# Arm 2 ("voxel"): old vs new microenv input, CANOPY_CLIM_MODE=voxel. This
# is real production behaviour today (per-pixel quantile system, median
# under stochastic=FALSE) -- get_colonization.R's .clim_voxel_slice() does
# NOT read the new mean-based cache in this pass (see its 2026-09-07
# correction comment), so this arm is testing whether Phase A's re-run
# produced a statistically similar per-pixel quantile summary to the
# archived pre-rerun one for the same site/domain -- not testing the new
# mean-based pathway, which isn't wired into any consumer yet.
#
# Same seeds both arms/both inputs: run_replicated() seeds each replicate
# as seed=rep (1,2,3, from n_reps=3 in the shared params file), identical
# regardless of which microenv or clim_mode is active.
source("scripts/02_model/config/patches.R")
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

old_path <- "data/processed/archive_pre_v7pix/microenv_MindoMirador_h0.40.rds"
new_path <- "data/processed/microenv_MindoMirador_h0.40.rds"
if (!file.exists(old_path)) stop("Missing archived old-input manifest: ", old_path)
if (!file.exists(new_path)) stop("Missing new (Phase A) manifest: ", new_path, " -- has the re-run written MindoMirador yet?")

run_arm <- function(microenv_path, tag, clim_mode) {
  Sys.setenv(CANOPY_MICROENV_OVERRIDE = microenv_path, CANOPY_CLIM_MODE = clim_mode)
  status <- system2("Rscript",
    c("scripts/02_model/run/run_colonization.R", "MindoMirador",
      "data/params/b4_equivalence_params.rds", tag, "0.4"))
  if (status != 0) stop("run_colonization.R failed for arm: ", tag)
  readRDS(sprintf("data/processed/colonization_MindoMirador_%s_h0.40.rds", tag))
}

final_stats <- function(result) {
  s <- result$summary
  final_t <- max(s$t)
  tapply(s$totalA[s$t == final_t], s$rep[s$t == final_t], sum)
}

run_gate <- function(clim_mode, label) {
  message(sprintf("=== Arm '%s' (CANOPY_CLIM_MODE=%s) ===", label, clim_mode))
  message("Running OLD-input arm...")
  old_result <- run_arm(old_path, sprintf("b4old_%s", label), clim_mode)
  message("Running NEW-input arm...")
  new_result <- run_arm(new_path, sprintf("b4new_%s", label), clim_mode)

  old_final <- final_stats(old_result)
  new_final <- final_stats(new_result)

  cat(sprintf("\n--- %s: mean/SD final-timestep total adult abundance across 3 reps ---\n", label))
  cat(sprintf("  old: mean=%.3f sd=%.3f | new: mean=%.3f sd=%.3f\n",
              mean(old_final), sd(old_final), mean(new_final), sd(new_final)))

  tol <- max(sd(old_final), 1e-6)
  passed <- abs(mean(old_final) - mean(new_final)) <= 2 * tol
  cat(sprintf("GATE (%s): %s (|mean diff|=%.3f, 2x old-run SD=%.3f)\n",
              label, if (passed) "PASSED" else "FAILED",
              abs(mean(old_final) - mean(new_final)), 2 * tol))
  passed
}

pooled_ok <- run_gate("pooled", "pooled")
voxel_ok  <- run_gate("voxel",  "voxel")

cat(sprintf("\n=== B4 OVERALL: %s ===\n",
            if (pooled_ok && voxel_ok) "PASSED (both arms)" else "FAILED"))
if (!(pooled_ok && voxel_ok)) {
  quit(status = 1)
}
