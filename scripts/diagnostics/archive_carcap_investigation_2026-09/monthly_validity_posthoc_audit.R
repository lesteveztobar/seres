# monthly_validity_posthoc_audit.R -- post-hoc: for each of the 7
# COMPLETED (archived, pre-gridsnap) sites, re-derive the monthly-sampled
# valid-cell count from the already-cached ERA5 data on disk (get_weather()
# reads the cached .nc, no new download/compute), and compare against the
# hour-1-only count already recorded in that site's own archived manifest.
# No microclimf model run involved.
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

sites_df <- make_sites(OBSERVATIONS_CSV, pad = 0.15)
mycredentials <- readRDS(file.path(BASE_DIR, "credentials.rds"))

rows <- list()
for (s in c("Maquipucuna","Mashpi","Yanayacu","MindoMirador","MindoTarabita","Saloya","LaElenita")) {
  site <- sites_df[sites_df$Site == s, ]
  if (nrow(site) == 0) next
  site <- split(site, seq_len(nrow(site)))[[1]]
  manifest_path <- file.path(PROCESSED_DIR, "archive_pre_gridsnap", sprintf("microenv_%s_h0.40.rds", s))
  hour1_count <- NA_integer_; hour1_total <- NA_integer_
  if (file.exists(manifest_path)) {
    m <- readRDS(manifest_path)
    hour1_count <- if (is.null(m$.n_valid_era5_cells)) NA_integer_ else m$.n_valid_era5_cells
    hour1_total <- if (is.null(m$.n_total_era5_cells)) NA_integer_ else m$.n_total_era5_cells
  }
  raster <- terra::rast(nrows = 2, ncols = 2, xmin = site$lon_min, xmax = site$lon_max,
                        ymin = site$lat_min, ymax = site$lat_max, crs = "EPSG:4326")
  terra::values(raster) <- 1
  tme_obs_end <- as.POSIXlt(site$tme_end, tz = "UTC")
  tme_to <- as.POSIXlt(format(seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
    by = "month", length.out = 2L)[2L] - 3600L, "%Y-%m-%d %H:00:00"), tz = "UTC")
  tme_from <- as.POSIXlt(format(seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
    by = "-1 month", length.out = 12L)[12L], "%Y-%m-01 00:00:00"), tz = "UTC")
  tme <- as.POSIXlt(seq(from = tme_from, to = tme_to, by = "hour"), tz = "UTC")
  site_era5_dir <- file.path(RAW_DIR, s, "era5")
  weatherdata <- tryCatch(
    get_weather(site = site, credentials = mycredentials, r = raster, tme = tme,
               dir = site_era5_dir, output = "grid"),
    error = function(e) { message(s, ": get_weather() failed -- ", conditionMessage(e)); NULL })
  if (is.null(weatherdata)) { rows[[length(rows)+1]] <- data.frame(site=s, hour1_valid=hour1_count, hour1_total=hour1_total, monthly_valid=NA, monthly_total=NA); next }
  wtemp <- if (inherits(weatherdata$temp, "PackedSpatRaster")) terra::unwrap(weatherdata$temp) else weatherdata$temp
  n_time <- terra::nlyr(wtemp)
  month_idx <- unique(round(seq(1, n_time, length.out = 12)))
  temp_samples <- terra::values(wtemp[[month_idx]])
  valid_monthly <- apply(temp_samples, 1, function(x) all(is.finite(x)))
  rows[[length(rows)+1]] <- data.frame(
    site = s, hour1_valid = hour1_count, hour1_total = hour1_total,
    monthly_valid = sum(valid_monthly), monthly_total = length(valid_monthly))
}
out <- do.call(rbind, rows)
write.csv(out, "output/monthly_validity_posthoc_audit.csv", row.names = FALSE)
cat("\n=== POST-HOC: hour-1-only vs 12-timestep monthly-sampled valid-cell counts ===\n")
print(out, row.names = FALSE)
