# raster_extent_overlap.R
# Item 0 (blocking, per instruction): per-site microclimate raster extent,
# elevation range from the DTM, confirmation of lookup_climate_by_height()'s
# pooling behavior, and a 7x7 pairwise overlap matrix -- to decide whether
# E2's 6-7 sites are independent samples or heavily-overlapping regional
# averages.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

sites <- c("Maquipucuna","Mashpi","MindoMirador","MindoTarabita","Saloya","Yanayacu","LaElenita")

info <- list()
for (s in sites) {
  m <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", s)))
  sp <- m$.spatial
  ext <- sp$ext
  lat0 <- mean(c(ext["ymin"], ext["ymax"]))
  lon0 <- mean(c(ext["xmin"], ext["xmax"]))
  ext_lat_deg <- ext["ymax"] - ext["ymin"]
  ext_lon_deg <- ext["xmax"] - ext["xmin"]
  ext_lat_km <- ext_lat_deg * 111
  ext_lon_km <- ext_lon_deg * 111 * cos(lat0 * pi / 180)
  cell_lat_m <- ext_lat_deg * 111000 / sp$nrow
  cell_lon_m <- ext_lon_deg * 111000 * cos(lat0 * pi / 180) / sp$ncol
  info[[s]] <- list(
    lat0 = lat0, lon0 = lon0, ext = ext,
    ext_lat_deg = ext_lat_deg, ext_lon_deg = ext_lon_deg,
    ext_lat_km = ext_lat_km, ext_lon_km = ext_lon_km,
    nrow = sp$nrow, ncol = sp$ncol,
    cell_lat_m = cell_lat_m, cell_lon_m = cell_lon_m
  )
}

cat("=== Per-site raster extent ===\n")
summary_df <- do.call(rbind, lapply(names(info), function(s) {
  x <- info[[s]]
  data.frame(site = s, lat0 = round(x$lat0, 4), lon0 = round(x$lon0, 4),
             ext_lat_deg = round(x$ext_lat_deg, 4), ext_lon_deg = round(x$ext_lon_deg, 4),
             ext_lat_km = round(x$ext_lat_km, 2), ext_lon_km = round(x$ext_lon_km, 2),
             nrow = x$nrow, ncol = x$ncol,
             cell_lat_m = round(x$cell_lat_m, 2), cell_lon_m = round(x$cell_lon_m, 2))
}))
print(summary_df, row.names = FALSE)
write.csv(summary_df, "output/raster_extent_summary.csv", row.names = FALSE)

# ── Elevation range from the DTM, per site ──────────────────────────────
cat("\n=== Elevation range (from DTM) per site ===\n")
elev_rows <- list()
for (s in sites) {
  dtm_path <- file.path(RAW_DIR, s, "dtm", "dtm.tif")
  if (!file.exists(dtm_path)) {
    # Fall back to any dtm-like file in that site's raw dir
    cand <- list.files(file.path(RAW_DIR, s), pattern = "dtm.*\\.tif$", full.names = TRUE, ignore.case = TRUE)
    if (length(cand) > 0) dtm_path <- cand[1]
  }
  if (!file.exists(dtm_path)) {
    message(s, ": no DTM file found at ", dtm_path, " -- skipping elevation range")
    elev_rows[[s]] <- data.frame(site = s, elev_min = NA, elev_median = NA, elev_max = NA)
    next
  }
  vals <- tryCatch(terra::values(terra::rast(dtm_path), na.rm = TRUE), error = function(e) numeric(0))
  if (length(vals) == 0) {
    elev_rows[[s]] <- data.frame(site = s, elev_min = NA, elev_median = NA, elev_max = NA)
  } else {
    elev_rows[[s]] <- data.frame(site = s, elev_min = min(vals), elev_median = median(vals), elev_max = max(vals))
  }
}
elev_df <- do.call(rbind, elev_rows)
print(elev_df, row.names = FALSE)
write.csv(elev_df, "output/raster_elevation_range.csv", row.names = FALSE)

# ── Pooling behavior confirmation ───────────────────────────────────────
cat("\n=== lookup_climate_by_height() pooling behavior ===\n")
cat("Confirmed by direct code trace (get_colonization.R, lookup_climate_by_height(),\n")
cat("header comment at the function definition):\n")
cat("  \"Spatially averaged across the raster at that height -> one value per\n")
cat("   hourly timestep. Reads the pre-averaged hourly-mean fields ... written\n")
cat("   by run_microclimate_site.R's write-time reduction.\"\n")
cat("This is the FULL raster, not footprint pixels -- there is no footprint-\n")
cat("restricted read path for this per-height climate series anywhere in\n")
cat("get_colonization.R. (voxel_background_table()/build_clim_cache_voxel()\n")
cat("are a SEPARATE, footprint-aware pathway used only for niche background\n")
cat("sampling -- not what E2/climate_variation_relative_height.R reads.)\n")

# ── Pairwise overlap matrix ──────────────────────────────────────────────
cat("\n=== Pairwise raster overlap (% of each site's OWN area) ===\n")
rect_overlap_area_km2 <- function(a, b) {
  # Rectangle overlap in lon/lat degrees, converted to km^2 using each
  # rectangle's own latitude for the lon->km conversion (approximate, fine
  # at this scale).
  ox1 <- max(a$ext["xmin"], b$ext["xmin"]); ox2 <- min(a$ext["xmax"], b$ext["xmax"])
  oy1 <- max(a$ext["ymin"], b$ext["ymin"]); oy2 <- min(a$ext["ymax"], b$ext["ymax"])
  if (ox2 <= ox1 || oy2 <= oy1) return(0)
  lat_mid <- (oy1 + oy2) / 2
  ((ox2 - ox1) * 111 * cos(lat_mid * pi / 180)) * ((oy2 - oy1) * 111)
}
area_km2 <- function(x) x$ext_lat_km * x$ext_lon_km

overlap_pct <- matrix(NA_real_, nrow = length(sites), ncol = length(sites), dimnames = list(sites, sites))
for (i in sites) for (j in sites) {
  ov <- rect_overlap_area_km2(info[[i]], info[[j]])
  overlap_pct[i, j] <- 100 * ov / area_km2(info[[i]])  # % of ROW site's own area
}
print(round(overlap_pct, 1))
write.csv(overlap_pct, "output/raster_overlap_pct_matrix.csv")

cat("\nDone.\n")
