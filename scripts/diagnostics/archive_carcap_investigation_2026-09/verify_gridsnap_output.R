# verify_gridsnap_output.R -- confirms the "pre-snap results are still
# valid" reasoning (0% pooled fallback + home cell always valid) actually
# holds once real output exists, rather than resting on it untested.
# MindoMirador: 1 site, 3 replicates, pooled mode (deterministic given
# stochastic=FALSE, so "same seeds" is automatic), old (pre-gridsnap,
# already-approved-for-Phase-F) manifest vs new (grid-snapped +
# nearest-valid-cell-substitution) manifest.
source("scripts/02_model/config/patches.R")
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

SITE <- "MindoMirador"
old_path <- "data/processed/microenv_MindoMirador_h0.40.rds"       # what Phase F actually used
new_path <- "data/processed/archive_pre_gridsnap_v2/microenv_MindoMirador_h0.40.rds"  # copied in before overwrite, see run_verify_gridsnap_output.sh

run_arm <- function(microenv_path, tag) {
  Sys.setenv(CANOPY_MICROENV_OVERRIDE = microenv_path, CANOPY_CLIM_MODE = "voxel")
  status <- system2("Rscript",
    c("scripts/02_model/run/run_colonization.R", SITE,
      "data/params/b4_equivalence_params.rds", tag, "0.4"))
  if (status != 0) stop("run_colonization.R failed for arm: ", tag)
  readRDS(sprintf("data/processed/colonization_%s_%s_h0.40.rds", SITE, tag))
}
final_stats <- function(result) {
  s <- result$summary
  final_t <- max(s$t)
  tapply(s$totalA[s$t == final_t], s$rep[s$t == final_t], sum)
}

message("Running OLD (pre-gridsnap) arm...")
old_result <- run_arm(old_path, "gridsnap_verify_old")
old_final <- final_stats(old_result)

message("Running NEW (grid-snapped) arm...")
new_result <- run_arm(new_path, "gridsnap_verify_new")
new_final <- final_stats(new_result)

cat(sprintf("\n=== GRID-SNAP VERIFICATION: %s, pooled/voxel mode, old vs new manifest ===\n", SITE))
cat(sprintf("OLD (pre-gridsnap): mean=%.3f sd=%.3f (per-rep: %s)\n",
            mean(old_final), sd(old_final), paste(round(old_final, 1), collapse = ", ")))
cat(sprintf("NEW (grid-snapped): mean=%.3f sd=%.3f (per-rep: %s)\n",
            mean(new_final), sd(new_final), paste(round(new_final, 1), collapse = ", ")))
match <- abs(mean(old_final) - mean(new_final)) < 1e-6
cat(sprintf("\nRESULT: %s (mean diff = %.6f)\n",
            if (match) "MATCH -- Phase F's use of pre-gridsnap manifests is confirmed safe."
            else "MISMATCH -- Phase F's batch needs to be re-run against grid-snapped data.",
            mean(old_final) - mean(new_final)))
