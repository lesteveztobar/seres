# carcap_audit.R -- LOOKUP ONLY, no colonization runs. Reports the
# RECALIBRATED carrying capacity (2026-09-10) per site: dimensional fix in
# build_forest() (crown branch area once per tree, distributed across the
# tree's crown voxels), plus occupiable_bark_fraction and
# maxillariinae_community_share. Rebuilds each site's forest exactly as
# init_colonization() does (same xDim/yDim/zDim, same forestparams, same
# build_forest()), fixed seed, so the audit is reproducible.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("carcap_audit_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))))

SITES <- c("Maquipucuna", "Mashpi", "Yanayacu", "MindoMirador", "MindoTarabita", "Saloya")
RESOLUTION <- 10

niches_all <- load_observations()
niches_all <- niches_all[!is.na(niches_all$lat) & !is.na(niches_all$lon) &
                         !is.na(niches_all$Height_m) & !is.na(niches_all$FinalID), ]
field_total <- nrow(niches_all)

fp0 <- default_forestparams()
cat(sprintf("Field observations total (all sites, Maxillariinae, identified): %d\n", field_total))
cat(sprintf("Forest params: branch_density=%.1f  trunk_r=%.3f  epiphyte_footprint_m2=%.3f\n", fp0$branch_density, fp0$trunk_r, fp0$epiphyte_footprint_m2))
cat(sprintf("  occupiable_bark_fraction=%.3f   maxillariinae_community_share=%.4f\n", fp0$occupiable_bark_fraction, fp0$maxillariinae_community_share))
cat(sprintf("  => per-tree Maxillariinae capacity from a mean crown (r=%.1f m): pi*r^2*bd*obf*msh/fp = %.2f\n\n",
            fp0$mean_crown_r, pi * fp0$mean_crown_r^2 * fp0$branch_density * fp0$occupiable_bark_fraction * fp0$maxillariinae_community_share / fp0$epiphyte_footprint_m2))

rows <- list()
for (site in SITES) {
  f <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site))
  if (!file.exists(f)) { message("skip ", site, " -- no manifest"); next }
  microenv <- readRDS(f)
  heights  <- microenv_heights(microenv)
  zDim     <- length(heights)
  so <- niches_all[niches_all$Area_or_Site == site, ]
  if (nrow(so) == 0) { message("skip ", site, " -- no obs"); next }
  # landscape extent from ALL raw observations (matches init_colonization())
  bb <- site_landscape_bbox(site)
  if (is.null(bb)) bb <- list(lat = range(so$lat), lon = range(so$lon), n = nrow(so))

  lat_range_m <- diff(bb$lat) * 111000
  lon_range_m <- diff(bb$lon) * 111000 * cos(mean(bb$lat) * pi / 180)
  xDim <- max(round(lon_range_m / RESOLUTION), 10) + 4
  yDim <- max(round(lat_range_m / RESOLUTION), 10) + 4
  area_ha <- (lat_range_m * lon_range_m) / 10000
  ground_m2 <- xDim * yDim * RESOLUTION^2

  ceiling <- site_canopy_ceiling(site, niches_all, max(heights))
  fp <- site_forestparams(site, ceiling)

  set.seed(1)
  landscape <- array(FALSE, dim = c(xDim, yDim, zDim))
  forest <- build_forest(landscape, heights, fp, so, RESOLUTION, land_bbox = bb)
  cc <- forest$carCap_voxel
  lz <- forest$landscape

  cap_occ <- cc[cc > 0L]
  total_K <- sum(cc)
  qs <- if (length(cap_occ)) quantile(cap_occ, c(0, .5, .9, .99, 1)) else rep(0, 5)

  rows[[site]] <- data.frame(
    site = site, xyz = sprintf("%dx%dx%d", xDim, yDim, zDim),
    area_ha = round(area_ha, 2), n_trees = forest$n_trees,
    n_valid_voxels = sum(lz), n_K_voxels = sum(cc > 0L),
    K_vox_median = qs[2], K_vox_p90 = qs[3], K_vox_max = qs[5],
    total_K = total_K,
    K_per_tree = round(total_K / forest$n_trees, 2),
    K_per_m2_ground = round(total_K / ground_m2, 3),
    K_per_ha = round(total_K / (area_ha), 1))
  cat(sprintf("%-14s trees=%d  area=%.2f ha  valid voxels=%d  K-carrying voxels=%d (%.1f%%)\n",
              site, forest$n_trees, area_ha, sum(lz), sum(cc > 0L), 100 * sum(cc > 0L) / max(1, sum(lz))))
  cat(sprintf("   TOTAL K = %s | %.2f per tree | %.3f per m2 ground | %.1f per ha | K/voxel median=%.0f p90=%.0f max=%.0f\n\n",
              format(total_K, big.mark = ","), total_K / forest$n_trees,
              total_K / ground_m2, total_K / area_ha, qs[2], qs[3], qs[5]))
}

df <- do.call(rbind, rows)
write.csv(df, file.path(OUTPUT_DIR, "carcap_audit.csv"), row.names = FALSE)

# ── Crown-scaling verification: K/tree should track ceiling^2 ─────────────
cat("\n================  CROWN-SCALING CHAIN (verification)  ================\n")
cat(sprintf("%-14s %8s %9s %9s %10s %9s %11s\n",
            "site", "ceiling", "mean_hgt", "crown_r", "crownArea", "K/tree", "pred K/tree"))
ref <- NULL
scl <- list()
for (site in names(rows)) {
  microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site)))
  heights <- microenv_heights(microenv)
  ceil <- site_canopy_ceiling(site, niches_all, max(heights))
  fp <- site_forestparams(site, ceil)
  ca <- pi * fp$mean_crown_r^2
  kpt <- df$K_per_tree[df$site == site]
  scl[[site]] <- list(ceil = ceil, mh = fp$mean_hgt, cr = fp$mean_crown_r, ca = ca, kpt = kpt)
  if (site == "Maquipucuna") ref <- scl[[site]]
}
for (site in names(scl)) {
  x <- scl[[site]]
  pred <- ref$kpt * (x$ceil / ref$ceil)^2
  cat(sprintf("%-14s %8.1f %9.2f %9.2f %10.1f %9.2f %11.2f\n",
              site, x$ceil, x$mh, x$cr, x$ca, x$kpt, pred))
}
cat("(pred K/tree = Maquipucuna's K/tree x (ceiling / 19.0)^2; deviations come from\n",
    " per-tree integer rounding in build_forest(), which lifts small-K sites more.)\n")

cat("\n================  RECALIBRATED CARRYING CAPACITY -- CROSS-SITE  ================\n")
print(df[, c("site", "n_trees", "n_valid_voxels", "n_K_voxels", "total_K", "K_per_tree", "K_per_m2_ground")], row.names = FALSE)
cat(sprintf("\nField total across ALL 7 sites: %d individuals.\n", field_total))
cat(sprintf("Model total K, summed over the 6 modelled sites: %s.\n", format(sum(df$total_K), big.mark = ",")))
cat(sprintf("(Pre-recalibration total over the same 6 sites was 169,012,168 -- a %.0fx reduction.)\n",
            169012168 / sum(df$total_K)))
cat("\nDone.\n")
