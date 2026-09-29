source("scripts/02_model/config/patches.R")
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")
niches <- load_observations(OBSERVATIONS_CSV)
sites <- unique(niches$Area_or_Site[!is.na(niches$Area_or_Site)])
sites <- sites[sites != "LaElenita"]
VARS <- c("temp", "relhum", "swdown")
for (s in sites) {
  mp <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", s))
  if (!file.exists(mp)) next
  microenv <- readRDS(mp)
  obs_site <- niches[niches$Area_or_Site == s, ]
  footprint <- if (.spatial_extent_usable(microenv)) {
    px <- .lonlat_to_pixel(obs_site$lon, obs_site$lat, microenv$.spatial)
    unique(data.frame(row = px$row, col = px$col))
  } else NULL
  cc_voxel <- build_clim_cache_voxel(microenv, footprint = footprint)
  bg <- as.data.frame(voxel_background_table(cc_voxel, microenv, footprint, VARS))
  cat(sprintf("%-15s n_rows=%5d  ", s, nrow(bg)))
  for (v in VARS) {
    n_fin <- sum(is.finite(bg[[v]]))
    cat(sprintf("%s: %d/%d finite  ", v, n_fin, nrow(bg)))
  }
  cat("\n")
}
