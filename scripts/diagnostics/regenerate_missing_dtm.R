# regenerate_missing_dtm.R
# Re-fetches data/raw/<site>/dtm/dtm.tif for any site missing it, WITHOUT
# rerunning the rest of the microclimate pipeline (weather/albedo/LAI/point
# model). Needed because a 2026-07-25 raw-data cleanup (intended to reclaim
# space from sites whose microenv generation was already fully complete)
# deleted dtm.tif along with the rest of those sites' raw non-.rds files --
# missing climate_variation_test.R's elevation-helper dependency on it, which is separate from and
# runs after microenv generation. get_dtm() (scripts/01_microclimate/lib.R)
# itself already skips any site whose dtm.tif is still present, so this is
# safe to rerun for all sites at any time.
#
# Usage: Rscript scripts/diagnostics/regenerate_missing_dtm.R
# Run from: /home/s38leste_hpc/seres/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(readr)
library(terra)
library(microclimdata)
source("scripts/02_model/config/paths.R")
source("scripts/01_microclimate/lib.R")

sites <- make_sites(OBSERVATIONS_CSV, pad = 0.15)

for (i in seq_len(nrow(sites))) {
  site <- sites[i, ]
  dtm_dir  <- file.path(RAW_DIR, site$Site, "dtm")
  dtm_path <- file.path(dtm_dir, "dtm.tif")
  if (file.exists(dtm_path)) {
    message(site$Site, ": dtm.tif already present -- skipping.")
    next
  }
  message(site$Site, ": dtm.tif missing -- regenerating...")
  dir.create(dtm_dir, recursive = TRUE, showWarnings = FALSE)
  raster <- terra::rast(
    nrows = 2, ncols = 2,
    xmin  = site$lon_min, xmax = site$lon_max,
    ymin  = site$lat_min, ymax = site$lat_max,
    crs   = "EPSG:4326"
  )
  terra::values(raster) <- 1
  get_dtm(r = raster, dir = dtm_dir, mask = FALSE)
  message(site$Site, ": done.")
}
