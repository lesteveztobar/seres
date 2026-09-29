source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg()

set.seed(42)
sites <- c("Maquipucuna","Mashpi","MindoMirador","MindoTarabita","Saloya","Yanayacu")
niches_raw <- read.csv(OBSERVATIONS_CSV)  # raw, pre-filter, for lat/lon extent (site_obs isn't species-filtered for landscape purposes)

for (s in sites) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", s))
  if (!file.exists(microenv_path)) { cat(s, ": no microenv\n"); next }
  microenv <- readRDS(microenv_path)
  heights <- microenv_heights(microenv)
  site_obs <- niches_raw[niches_raw$Area_or_Site == s & !is.na(niches_raw$lat) & !is.na(niches_raw$lon), ]
  if (nrow(site_obs) == 0) { cat(s, ": no raw observations at all for landscape extent\n"); next }

  # v7 (Phase 2.2 sanity gate, revised): one MAX-based ceiling drives
  # zDim (via the microenv, already built), canopy_z, AND mean_hgt/sd_hgt
  # (site_forestparams()'s ceiling-scaled rule) -- computed once here and
  # threaded through, not three independently-sourced quantities.
  ceiling <- site_canopy_ceiling(s, niches_raw, max(heights))
  forestparams <- site_forestparams(s, ceiling)
  cat(sprintf("%-15s ceiling(max)=%.2fm -> forestparams: mean_hgt=%.2fm sd_hgt=%.2fm\n",
              s, ceiling, forestparams$mean_hgt, forestparams$sd_hgt))

  resolution <- 10
  lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
  lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 * cos(mean(site_obs$lat) * pi / 180)
  xDim <- max(round(lon_range_m / resolution), 10) + 4
  yDim <- max(round(lat_range_m / resolution), 10) + 4
  zDim <- length(heights)
  landscape <- array(FALSE, dim = c(xDim, yDim, zDim))

  forest <- build_forest(landscape, heights, forestparams, site_obs, resolution)
  max_tree_height <- max(forest$trees$height)
  occupied_z <- which(apply(forest$landscape, 3, any))
  n_zero_canopy_tiers <- zDim - length(occupied_z)
  area_ha <- (lon_range_m * lat_range_m) / 10000

  # 2026-09-06 (item 3, revised metric): "zero-canopy tiers" is arithmetic
  # under the ceiling-scaled forestparams rule, not validation -- with
  # mean_hgt=ratio*ceiling and sd_hgt=cv*mean_hgt (ratio=0.442, cv=0.417),
  # the ceiling sits (1/ratio - 1)/cv = (1/0.442-1)/0.417 = 3.02 SD above
  # mean tree height at EVERY site by construction, so with thousands of
  # trees, a >=3sigma draw at the very top tier is routine and "exactly 1
  # empty tier everywhere" was guaranteed by the formula, not a property of
  # any one site's data. Replaced with a metric that CAN differ between
  # sites: the height (as a fraction of the site's own ceiling) at which
  # generated canopy area falls below 10% of that site's own maximum
  # per-tier canopy area.
  canopy_area_by_z <- apply(forest$landscape, 3, sum)  # occupied-voxel count per height tier
  max_area <- max(canopy_area_by_z)
  below_10pct <- which(canopy_area_by_z < 0.10 * max_area)
  # Height at which the profile FIRST drops below 10% and stays of interest
  # -- report the lowest height (closest to ground) where this threshold is
  # crossed going up, i.e. the top of the "substantial canopy" zone.
  above_10pct <- which(canopy_area_by_z >= 0.10 * max_area)
  z_10pct <- if (length(above_10pct) > 0) max(above_10pct) + 1 else 1
  z_10pct <- min(z_10pct, zDim)
  height_10pct <- heights[z_10pct]
  frac_of_ceiling <- height_10pct / ceiling

  cat(sprintf("%-15s ceiling(max)=%.1fm | max_simulated_tree=%.1fm | zDim=%d | occupied_z_tiers=%d | zero-canopy tiers=%d | area=%.1fha | n_trees=%d | height@10%%-of-max-canopy-area=%.2fm (%.1f%% of ceiling)\n",
    s, ceiling, max_tree_height, zDim, length(occupied_z), n_zero_canopy_tiers,
    area_ha, nrow(forest$trees), height_10pct, 100 * frac_of_ceiling))
}
