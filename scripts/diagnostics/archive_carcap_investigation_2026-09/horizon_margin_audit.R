# horizon_margin_audit.R -- quantifies the horizon-shading cost of the
# ERA5 grid-snap (2026-09-08), per the author's request: rather than
# growing the domain to avoid a thin margin, compute horizon angles at
# each site's own coordinates and see whether the margin loss actually
# matters.
#
# Uses the EXISTING cached DTM (data/raw/<site>/dtm/dtm.tif -- the OLD,
# pre-snap domain, ~16.7km symmetric margin every side) since the NEW
# (post-snap) DTM hasn't been re-fetched yet (that's what the currently-
# running Phase A jobs are doing). The relevant comparison for the 4
# shifted sites is: how much does the computed horizon angle change if
# terrain is truncated at the NEW domain's own (asymmetric) margin instead
# of the OLD symmetric one? For directions where the new margin is LARGER
# than 16.7km (the "thick" side of the shift), more terrain becomes
# available, not less -- horizon estimates there can only improve or stay
# the same, so this only computes the THIN-side truncation, the actual
# risk.
suppressMessages(library(terra))

sites <- list(
  MindoMirador  = list(lat = -0.0218, lon = -78.7622, old_margin_km = 16.7, new_thin_km = 5.21, thin_az = 0),    # north
  LaElenita     = list(lat = -0.0202, lon = -78.7760, old_margin_km = 16.7, new_thin_km = 5.66, thin_az = 90),   # east
  MindoTarabita = list(lat = -0.0820, lon = -78.7609, old_margin_km = 16.7, new_thin_km = 11.92, thin_az = 0),   # north
  Saloya        = list(lat = 0.0092,  lon = -78.8248, old_margin_km = 16.7, new_thin_km = 11.03, thin_az = 90)   # east
)

km_per_deg_lat <- 111.32
azimuths <- seq(0, 330, by = 30)  # 12 directions, degrees from north, clockwise

horizon_angle <- function(dtm, site_lat, site_lon, azimuth_deg, max_dist_km, step_km = 0.2) {
  km_per_deg_lon <- 111.32 * cos(site_lat * pi / 180)
  site_elev <- terra::extract(dtm, cbind(site_lon, site_lat))[1, 1]
  if (is.na(site_elev)) return(NA_real_)
  dists <- seq(step_km, max_dist_km, by = step_km)
  az_rad <- azimuth_deg * pi / 180
  lats <- site_lat + (dists / km_per_deg_lat) * cos(az_rad)
  lons <- site_lon + (dists / km_per_deg_lon) * sin(az_rad)
  elevs <- terra::extract(dtm, cbind(lons, lats))[, 1]
  ok <- !is.na(elevs)
  if (!any(ok)) return(NA_real_)
  angles <- atan2(elevs[ok] - site_elev, dists[ok] * 1000) * 180 / pi
  max(angles, na.rm = TRUE)
}

rows <- list()
for (s in names(sites)) {
  d <- sites[[s]]
  dtm_path <- sprintf("data/raw/%s/dtm/dtm.tif", s)
  if (!file.exists(dtm_path)) { message("Skipping ", s, " -- no cached DTM"); next }
  dtm <- terra::rast(dtm_path)
  for (az in azimuths) {
    h_old  <- horizon_angle(dtm, d$lat, d$lon, az, d$old_margin_km)
    h_thin <- horizon_angle(dtm, d$lat, d$lon, az, min(d$new_thin_km, d$old_margin_km))
    rows[[length(rows) + 1]] <- data.frame(
      site = s, azimuth = az, is_thin_direction = (az == d$thin_az),
      horizon_deg_old_margin = h_old, horizon_deg_thin_margin = h_thin,
      diff_deg = h_old - h_thin)
  }
}
out <- do.call(rbind, rows)
write.csv(out, "output/horizon_margin_audit.csv", row.names = FALSE)
cat("\n=== HORIZON ANGLE, OLD (16.7km) MARGIN vs NEW THIN-SIDE MARGIN, PER AZIMUTH ===\n")
print(out, row.names = FALSE)

cat("\n=== Summary: max |diff| per site, and specifically at the thin direction ===\n")
for (s in unique(out$site)) {
  sub <- out[out$site == s, ]
  thin_row <- sub[sub$is_thin_direction, ]
  cat(sprintf("%-15s max|diff| across all azimuths = %.2f deg | at the thin direction itself = %.2f deg\n",
              s, max(abs(sub$diff_deg), na.rm = TRUE),
              if (nrow(thin_row) > 0) thin_row$diff_deg[1] else NA))
}
