# grid_origin_sensitivity.R -- the "noise floor" measurement asked for in
# the option-3 decision. One site (MindoMirador -- footprint is only a
# handful of 90 m climate pixels, so it is the most aliasing-prone), four
# sub-pixel grid origins (the deterministic ERA5 snap at 0 m, plus 22, 45
# and 68 m single-axis shifts already built as separate microclimate
# manifests), 3 replicates each, voxel mode, identical seeds across all
# arms. Reports the spread in (a) final total abundance and (b) the
# abundance-weighted realised-height distribution, BETWEEN grid origins
# versus BETWEEN replicates of one origin -- the former is the grid-origin
# noise floor, the latter the ordinary stochastic noise it must be judged
# against.
source("scripts/02_model/config/patches.R")
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

SITE     <- "MindoMirador"
PARAMS_F <- "data/params/b4_equivalence_params.rds"   # scalar params, n_reps=3, literature values
ARMS <- c(
  "0"  = file.path(PROCESSED_DIR, "microenv_MindoMirador_h0.40.rds"),
  "22" = file.path(PROCESSED_DIR, "microenv_MindoMirador_h0.40_gridoffset22m.rds"),
  "45" = file.path(PROCESSED_DIR, "microenv_MindoMirador_h0.40_gridoffset45m.rds"),
  "68" = file.path(PROCESSED_DIR, "microenv_MindoMirador_h0.40_gridoffset68m.rds"))
stopifnot(all(file.exists(ARMS)))

run_arm <- function(microenv_path, tag) {
  Sys.setenv(CANOPY_MICROENV_OVERRIDE = microenv_path, CANOPY_CLIM_MODE = "voxel")
  st <- system2("Rscript",
    c("scripts/02_model/run/run_colonization.R", SITE, PARAMS_F, tag, "0.4"))
  if (st != 0) stop("run_colonization.R failed for arm ", tag)
  readRDS(file.path(PROCESSED_DIR, sprintf("colonization_%s_%s_h0.40.rds", SITE, tag)))
}

# abundance-weighted realised heights for one run object (final timestep,
# S+J+A pooled) -- identical construction to held_out_validation.R part (c)
realised_heights <- function(run, heights) {
  if (is.null(run)) return(numeric(0))
  final_t <- dim(run$abundanceA)[4]
  out <- numeric(0)
  for (arr in list(run$abundanceA, run$abundanceJ, run$abundanceS)) {
    slice <- arr[, , , final_t, , drop = FALSE]
    occ <- which(slice > 0, arr.ind = TRUE)
    if (nrow(occ) == 0) next
    out <- c(out, rep(heights[occ[, 3]], times = slice[occ]))
  }
  out
}
final_total <- function(run) {
  if (is.null(run)) return(NA_real_)
  ft <- dim(run$abundanceA)[4]
  sum(run$abundanceA[, , , ft, ], run$abundanceJ[, , , ft, ], run$abundanceS[, , , ft, ])
}

microenv0 <- readRDS(ARMS[["0"]])
heights   <- microenv_heights(microenv0)

rows <- list(); hrows <- list()
for (tag in names(ARMS)) {
  message("=== arm offset ", tag, " m ===")
  res <- run_arm(ARMS[[tag]], sprintf("gridorigin_%sm", tag))
  for (i in seq_along(res$runs)) {
    run <- res$runs[[i]]
    tot <- final_total(run)
    h   <- realised_heights(run, heights)
    rows[[length(rows) + 1]] <- data.frame(
      offset_m = as.integer(tag), rep = i, final_total = tot,
      n_ind = length(h),
      mean_h = if (length(h)) mean(h) else NA_real_,
      median_h = if (length(h)) median(h) else NA_real_,
      p10_h = if (length(h)) quantile(h, .10, names = FALSE) else NA_real_,
      p90_h = if (length(h)) quantile(h, .90, names = FALSE) else NA_real_)
    if (length(h)) hrows[[length(hrows) + 1]] <- data.frame(
      offset_m = as.integer(tag), rep = i, height = h)
  }
}
df  <- do.call(rbind, rows)
hdf <- if (length(hrows)) do.call(rbind, hrows) else data.frame()
write.csv(df,  file.path(OUTPUT_DIR, "grid_origin_sensitivity_summary.csv"), row.names = FALSE)
write.csv(hdf, file.path(OUTPUT_DIR, "grid_origin_sensitivity_heights.csv"), row.names = FALSE)

cat("\n\n================  GRID-ORIGIN SENSITIVITY  ================\n")
cat(sprintf("Site: %s | params: %s | voxel mode | 4 grid origins x %d reps, identical seeds\n\n",
            SITE, PARAMS_F, max(df$rep)))
print(df, row.names = FALSE)

# per-arm means
arm_mean_tot <- tapply(df$final_total, df$offset_m, mean)
arm_mean_h   <- tapply(df$mean_h,      df$offset_m, function(x) mean(x, na.rm = TRUE))
# within-arm (between-replicate) sd, averaged over arms
within_tot <- mean(tapply(df$final_total, df$offset_m, sd), na.rm = TRUE)
within_h   <- mean(tapply(df$mean_h,      df$offset_m, function(x) sd(x, na.rm = TRUE)), na.rm = TRUE)

cat("\n-- Final total abundance --\n")
cat(sprintf("  per-origin means      : %s\n", paste(sprintf("%dm=%.1f", as.integer(names(arm_mean_tot)), arm_mean_tot), collapse = "  ")))
cat(sprintf("  BETWEEN-origin sd      : %.2f  (range %.1f)\n", sd(arm_mean_tot), diff(range(arm_mean_tot))))
cat(sprintf("  WITHIN-origin sd (reps): %.2f  (mean over origins)\n", within_tot))
cat(sprintf("  ratio between/within   : %.2f\n", sd(arm_mean_tot) / within_tot))

cat("\n-- Realised height (abundance-weighted mean, m) --\n")
cat(sprintf("  per-origin means      : %s\n", paste(sprintf("%dm=%.2f", as.integer(names(arm_mean_h)), arm_mean_h), collapse = "  ")))
cat(sprintf("  BETWEEN-origin sd      : %.3f  (range %.2f)\n", sd(arm_mean_h), diff(range(arm_mean_h))))
cat(sprintf("  WITHIN-origin sd (reps): %.3f  (mean over origins)\n", within_h))
if (nrow(hdf)) {
  cat("\n-- Pooled realised-height distribution per origin (quantiles) --\n")
  for (o in sort(unique(hdf$offset_m))) {
    hh <- hdf$height[hdf$offset_m == o]
    cat(sprintf("  %2dm: n=%5d  min=%.1f  p25=%.1f  median=%.1f  p75=%.1f  max=%.1f  mean=%.2f\n",
                o, length(hh), min(hh), quantile(hh,.25), median(hh), quantile(hh,.75), max(hh), mean(hh)))
  }
  ks <- suppressWarnings(ks.test(hdf$height[hdf$offset_m == 0], hdf$height[hdf$offset_m == 68]))
  cat(sprintf("\n  KS 0m vs 68m realised-height distributions: D=%.3f p=%.4f\n", ks$statistic, ks$p.value))
}
cat("\nDone.\n")
