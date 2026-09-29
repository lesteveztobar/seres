# home_cell_validity_lookup.R -- for each of 7 sites: which ERA5 cell
# contains the site's own (tiny, ~500m) landscape footprint, and is that
# SPECIFIC cell's hour-1 temperature finite (the first of runpointmodela()'s
# two validity conditions -- the other, vegp$hgt not NA, isn't checked here,
# see the report). Uses each site's CURRENT (grid-snapped) domain bounds,
# reads cached ERA5 data (no new download -- get_weather() hits its own
# on-disk cache).
suppressMessages({
  library(readr); library(rgee); library(mcera5); library(microclimf)
  library(microclimdata); library(terra); library(luna); library(parallel)
})
PYTHON_PATH <- Sys.getenv("CANOPY_PYTHON", unset = "/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python")
reticulate::use_python(PYTHON_PATH, required = TRUE)
source("scripts/02_model/config/patches.R")
source("scripts/01_microclimate/lib.R")
source("scripts/02_model/config/paths.R")
source("scripts/data_prep/helper_functions.R")

# Site coordinate (its own observation centroid -- the footprint is built
# around this) and current (post-snap) domain bounds, exactly as computed/
# applied by the currently-running Phase A jobs.
sites <- list(
  Maquipucuna   = list(lat = 0.1208,  lon = -78.6318,
                       lat_min = 0.1208-0.15245, lat_max = 0.1208+0.15245,
                       lon_min = -78.6318-0.15165, lon_max = -78.6318+0.15165),
  Mashpi        = list(lat = 0.1640,  lon = -78.8794,
                       lat_min = 0.1640-0.15285, lat_max = 0.1640+0.15285,
                       lon_min = -78.8794-0.15205, lon_max = -78.8794+0.15205),
  Yanayacu      = list(lat = -0.5947, lon = -77.8924,
                       lat_min = -0.5947-0.15315, lat_max = -0.5947+0.15315,
                       lon_min = -77.8924-0.15235, lon_max = -77.8924+0.15235),
  MindoMirador  = list(lat = -0.0218, lon = -78.7622,
                       lat_min = -0.2767, lat_max = 0.0267,
                       lon_min = -78.7622-0.15085, lon_max = -78.7622+0.15085),
  MindoTarabita = list(lat = -0.0820, lon = -78.7609,
                       lat_min = -0.2774, lat_max = 0.0275,
                       lon_min = -78.7609-0.15165, lon_max = -78.7609+0.15165),
  Saloya        = list(lat = 0.0092,  lon = -78.8248,
                       lat_min = 0.0092-0.15085, lat_max = 0.0092+0.15085,
                       lon_min = -79.0250, lon_max = -78.7250),
  LaElenita     = list(lat = -0.0202, lon = -78.7760,
                       lat_min = -0.0202-0.15085, lat_max = -0.0202+0.15085,
                       lon_min = -79.0250, lon_max = -78.7250)
)

mycredentials <- readRDS(file.path(BASE_DIR, "credentials.rds"))
sites_df <- make_sites(OBSERVATIONS_CSV, pad = 0.15)

rows <- list()
for (s in names(sites)) {
  d <- sites[[s]]
  site_row <- sites_df[sites_df$Site == s, ]
  if (nrow(site_row) == 0) next
  site_row <- split(site_row, seq_len(nrow(site_row)))[[1]]
  site_row$lat_min <- d$lat_min; site_row$lat_max <- d$lat_max
  site_row$lon_min <- d$lon_min; site_row$lon_max <- d$lon_max

  raster <- terra::rast(nrows = 2, ncols = 2, xmin = site_row$lon_min, xmax = site_row$lon_max,
                        ymin = site_row$lat_min, ymax = site_row$lat_max, crs = "EPSG:4326")
  terra::values(raster) <- 1
  tme_obs_end <- as.POSIXlt(site_row$tme_end, tz = "UTC")
  tme_to <- as.POSIXlt(format(seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
    by = "month", length.out = 2L)[2L] - 3600L, "%Y-%m-%d %H:00:00"), tz = "UTC")
  tme_from <- as.POSIXlt(format(seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
    by = "-1 month", length.out = 12L)[12L], "%Y-%m-01 00:00:00"), tz = "UTC")
  tme <- as.POSIXlt(seq(from = tme_from, to = tme_to, by = "hour"), tz = "UTC")
  site_era5_dir <- file.path(RAW_DIR, s, "era5")

  weatherdata <- tryCatch(
    get_weather(site = site_row, credentials = mycredentials, r = raster, tme = tme,
               dir = site_era5_dir, output = "grid"),
    error = function(e) { message(s, ": get_weather() failed -- ", conditionMessage(e)); NULL })
  if (is.null(weatherdata)) { rows[[length(rows)+1]] <- data.frame(site=s, note="get_weather failed"); next }
  wtemp <- if (inherits(weatherdata$temp, "PackedSpatRaster")) terra::unwrap(weatherdata$temp) else weatherdata$temp

  # The site's own home ERA5 grid point (nearest multiple of 0.25deg to its
  # own coordinate) -- its footprint sits inside this one cell.
  home_lat <- round(d$lat / 0.25) * 0.25
  home_lon <- round(d$lon / 0.25) * 0.25
  # Row/col of that grid point's CENTER within this run's raster.
  rc <- terra::rowColFromCell(wtemp, terra::cellFromXY(wtemp, cbind(home_lon, home_lat)))
  hour1_temp <- terra::extract(wtemp[[1]], cbind(home_lon, home_lat))[1, 1]
  rows[[length(rows) + 1]] <- data.frame(
    site = s, home_cell_lat = home_lat, home_cell_lon = home_lon,
    row = rc[1, 1], col = rc[1, 2],
    hour1_temp_value = hour1_temp, hour1_valid = is.finite(hour1_temp))
}
out <- do.call(rbind, rows)
write.csv(out, "output/home_cell_validity_lookup.csv", row.names = FALSE)
cat("\n=== Home-cell hour-1 temperature validity, all 7 sites ===\n")
print(out, row.names = FALSE)
