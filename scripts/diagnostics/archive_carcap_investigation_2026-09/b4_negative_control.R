# b4_negative_control.R -- B4 negative control (per the author's instruction,
# 2026-09-08). MindoMirador retracted as the test site: its live
# voxel_quantiles is only 0.03% populated (see the valid-grid-cells finding
# below), so a weak/absent response there would be indistinguishable from
# "nothing to perturb". Uses Maquipucuna instead -- 4/4 valid grid cells,
# matching its own historical baseline exactly, genuinely healthy data.
#
# Perturbs the NEW manifest's voxel_quantiles (+1C to every finite temp
# value, in every day/night/both/month/annual key -- this IS what "voxel"
# mode's .clim_voxel_slice() actually reads, not pixel_means, which no
# consumer reads yet), then re-runs the voxel arm against the perturbed
# copy and compares to an UNPERTURBED voxel-mode baseline run against the
# same (real, new) manifest, same seeds (run_replicated() seeds by rep,
# unaffected by which manifest is read) -- computed fresh here rather than
# reused from B4 (B4 never ran a Maquipucuna arm).
source("scripts/02_model/config/patches.R")
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

SITE      <- "Maquipucuna"
live_dir  <- sprintf("/lustre/scratch/data/s38leste_hpc-seres/microenv_%s_h0.40_heights", SITE)
pert_dir  <- sprintf("/lustre/scratch/data/s38leste_hpc-seres/microenv_%s_h0.40_heights_PERTURBED", SITE)
dir.create(pert_dir, showWarnings = FALSE)

new_manifest <- readRDS(sprintf("data/processed/microenv_%s_h0.40.rds", SITE))
heights <- new_manifest$.heights

already_perturbed <- all(file.exists(file.path(pert_dir, sprintf("h%.2f.rds", heights))))
if (already_perturbed) {
  message("Perturbed height files already exist from a previous attempt -- skipping re-write.")
} else {
  message(sprintf("Perturbing %d height files (+1C to every finite *_temp quantile)...", length(heights)))
  for (h in heights) {
    h_key <- sprintf("h%.2f.rds", h)
    src <- file.path(live_dir, h_key)
    dst <- file.path(pert_dir, h_key)
    hd <- readRDS(src)
    temp_keys <- grep("_temp$", names(hd$voxel_quantiles$quantiles), value = TRUE)
    for (k in temp_keys) {
      m <- hd$voxel_quantiles$quantiles[[k]]
      m[is.finite(m)] <- m[is.finite(m)] + 1.0
      hd$voxel_quantiles$quantiles[[k]] <- m
    }
    # Also perturb temp_mean (the pooled fallback series) for completeness --
    # not what this control targets, but keeps the perturbed file internally
    # consistent rather than a mismatched pooled/voxel pair.
    hd$temp_mean <- hd$temp_mean + 1.0
    saveRDS(hd, dst)
  }
  message("Perturbation done.")
}

pert_manifest <- new_manifest
pert_manifest$.height_dir <- pert_dir
pert_path <- sprintf("data/processed/microenv_%s_h0.40_PERTURBED.rds", SITE)
saveRDS(pert_manifest, pert_path)
message(sprintf("Perturbed manifest saved: %s", pert_path))

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

# 2026-09-08: the first attempt at this job TIMED OUT at 4h with the
# baseline arm already complete and saved -- reuse it rather than
# re-spending ~3h recomputing an identical result (same manifest, same
# seeds, deterministic given stochastic=FALSE).
base_out_path <- sprintf("data/processed/colonization_%s_b4_negctrl_base_h0.40.rds", SITE)
if (file.exists(base_out_path)) {
  message("Reusing already-completed UNPERTURBED baseline result (previous attempt's timeout)...")
  base_result <- readRDS(base_out_path)
} else {
  message("Running UNPERTURBED baseline (real new manifest, voxel mode)...")
  base_result <- run_arm(sprintf("data/processed/microenv_%s_h0.40.rds", SITE), "b4_negctrl_base")
}
base_final  <- final_stats(base_result)

message("Running PERTURBED arm (+1C temp, voxel mode)...")
pert_result <- run_arm(pert_path, "b4_negctrl_voxel")
pert_final  <- final_stats(pert_result)

cat(sprintf("\n=== B4 NEGATIVE CONTROL: %s, voxel mode, +1C temp perturbation ===\n", SITE))
cat(sprintf("Unperturbed (real new manifest): mean=%.3f sd=%.3f  (per-rep: %s)\n",
            mean(base_final), sd(base_final), paste(round(base_final, 1), collapse = ", ")))
cat(sprintf("Perturbed (+1C temp):            mean=%.3f sd=%.3f  (per-rep: %s)\n",
            mean(pert_final), sd(pert_final), paste(round(pert_final, 1), collapse = ", ")))

moved <- abs(mean(pert_final) - mean(base_final)) > 1e-6
cat(sprintf("\nRESULT: output %s after perturbation (mean diff = %.3f)\n",
            if (moved) "MOVED" else "DID NOT MOVE -- model is not reading the perturbed file",
            mean(pert_final) - mean(base_final)))
