# Task 4b -- elevation x canopy-height exchange rate. Task 4a (fill missing
# CanopyHeight_m from the Lang et al. vhgt.tif raster) is BLOCKED -- vhgt.tif
# is not persisted on disk anywhere in this project (confirmed: find over
# data/raw turns up nothing); it's fetched on demand via a Google Earth
# Engine export task per site (rgee::ee$batch$Export$image$toDrive(), see
# lib.R get_vegetation()) then downloaded from Drive -- the same EE
# authentication mechanism that had expired and needed a fresh interactive
# browser re-auth earlier this session, plus a real per-site GEE export/
# Drive-download wait (each one is a genuine remote job, not instant). Given
# time constraints this was not attempted. This script uses only the
# ALREADY-MEASURED CanopyHeight_m field -- no raster fill -- and reports
# sample size honestly at every step per the task's own instruction.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

niches <- load_observations()
cat("Total rows (post Maxillariinae filter):", nrow(niches), "\n")
cat("Rows with CanopyHeight_m:", sum(!is.na(niches$CanopyHeight_m)), "\n")
cat("Rows with Height_m:", sum(!is.na(niches$Height_m)), "\n")
niches <- augment_elevation(niches)
cat("Rows with Height_m + CanopyHeight_m + Elevation_final_m all present:",
    sum(!is.na(niches$Height_m) & !is.na(niches$CanopyHeight_m) & !is.na(niches$Elevation_final_m)), "\n\n")

SITES <- c("Maquipucuna","Mashpi","MindoMirador","MindoTarabita","Saloya","Yanayacu")
niches_f <- niches[!is.na(niches$lat) & !is.na(niches$lon) & !is.na(niches$Height_m), ]
sc <- .build_site_climate_series(SITES, niches = niches_f)
site_pixels <- sc$site_pixels
site_elev <- sc$site_elev

# Crossing height: for a target VPD-like value (proxy: site's own median
# relhum, since VPD isn't directly computed anywhere in this pipeline --
# flagged explicitly, not silently substituted) find the height at which
# each site's smoothed relhum profile crosses it.
crossing_height <- function(px, target, ceiling) {
  agg <- aggregate(relhum ~ height, data = px, FUN = median)
  agg <- agg[order(agg$height), ]
  # find first sign change of (relhum - target) -- linear interpolation
  d <- agg$relhum - target
  sign_change <- which(diff(sign(d)) != 0)
  if (length(sign_change) == 0) return(c(abs = NA_real_, rel = NA_real_))
  i <- sign_change[1]
  h1 <- agg$height[i]; h2 <- agg$height[i+1]
  d1 <- d[i]; d2 <- d[i+1]
  h_cross <- h1 + (0 - d1) * (h2 - h1) / (d2 - d1)
  c(abs = h_cross, rel = h_cross / ceiling)
}

# v7 (Phase 1.6): shared ceiling formula, site_canopy_ceiling() (shared_helpers.R) --
# was an independently-duplicated mean(...)-based formula here.
ceiling_for <- function(site, heights) site_canopy_ceiling(site, niches, max(heights))

# Target: overall median relhum across ALL sites pooled (one fixed,
# site-independent target value, per task wording "a target value of each
# climate variable, e.g. site median VPD" -- using relhum as the tractable
# proxy this pipeline actually computes).
all_relhum <- unlist(lapply(site_pixels, function(px) px$relhum))
target_relhum <- median(all_relhum, na.rm = TRUE)
cat("Target relhum (pooled median across all sites' pixel-height-hours):", round(target_relhum, 2), "\n\n")

# Primary threshold (2026-09-29): each site's OWN median relhum, so the
# crossing height measures where that site's profile passes its own midpoint
# rather than a pooled value some sites never reach. The pooled-median
# crossing is kept alongside (suffix _pooled) for comparison.
rows <- lapply(names(site_pixels), function(s) {
  ceiling <- ceiling_for(s, site_pixels[[s]]$height)
  site_target <- median(site_pixels[[s]]$relhum, na.rm = TRUE)
  cr <- crossing_height(site_pixels[[s]], site_target, ceiling)
  cp <- crossing_height(site_pixels[[s]], target_relhum, ceiling)
  data.frame(site = s, elevation = site_elev[s], ceiling = ceiling,
             target_relhum_site = site_target,
             crossing_height_abs = cr["abs"], crossing_height_rel = cr["rel"],
             target_relhum_pooled = target_relhum,
             crossing_height_abs_pooled = cp["abs"], crossing_height_rel_pooled = cp["rel"])
})
result_df <- do.call(rbind, rows)
cat("n sites with a usable crossing height:", sum(!is.na(result_df$crossing_height_abs)), "of", nrow(result_df), "\n")
print(result_df, row.names = FALSE)
write.csv(result_df, file.path(OUTPUT_DIR, "elevation_canopy_crossing_height.csv"), row.names = FALSE)

n_usable <- sum(!is.na(result_df$crossing_height_abs))
cat(sprintf("\nSample size for the elevation~crossing-height regression: n=%d sites.\n", n_usable))
if (n_usable < 4) {
  cat("n<4 -- NOT reporting a regression here; too small to defend a slope+CI. Reporting only the raw table above.\n")
} else {
  d <- result_df[!is.na(result_df$crossing_height_abs), ]
  fit_abs <- lm(crossing_height_abs ~ elevation, data = d)
  fit_rel <- lm(crossing_height_rel ~ elevation, data = d)
  cat("\n=== Absolute crossing height ~ elevation ===\n")
  print(summary(fit_abs)$coefficients)
  cat(sprintf("R^2=%.3f, n=%d\n", summary(fit_abs)$r.squared, nrow(d)))
  cat("\n=== Relative crossing height ~ elevation ===\n")
  print(summary(fit_rel)$coefficients)
  cat(sprintf("R^2=%.3f, n=%d\n", summary(fit_rel)$r.squared, nrow(d)))
  cat("\nSpearman (rank-based, more defensible at this n):\n")
  print(cor.test(d$elevation, d$crossing_height_abs, method="spearman"))
  print(cor.test(d$elevation, d$crossing_height_rel, method="spearman"))
}
