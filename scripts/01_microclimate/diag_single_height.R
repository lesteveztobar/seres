# diag_single_height.R — isolated single-height, single-worker timing/memory diagnostic.
# Not part of the pipeline. Reuses cached weather/DTM/vegetation/soil/point-model
# data for one site and calls runmicro() exactly once, so wall-clock and peak RSS
# (measured externally via `/usr/bin/time -v`) reflect ONE height, ONE core, no
# mclapply contention -- the baseline needed before sizing cpus-per-task/--mem/--time.
args <- commandArgs(trailingOnly = TRUE)
TARGET_SITE <- if (length(args) >= 1) args[1] else "LaElenita"
HEIGHT_STEP <- if (length(args) >= 2) as.numeric(args[2]) else 0.25
TEST_HEIGHT <- if (length(args) >= 3) as.numeric(args[3]) else 5.10

source("scripts/02_model/config/patches.R")
library(rgee)
library(readr)
library(mcera5)
library(microclimf)
library(microclimdata)
library(terra)
library(luna)

PYTHON_PATH <- Sys.getenv("CANOPY_PYTHON",
  unset = "/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python")
reticulate::use_python(PYTHON_PATH, required = TRUE)

source("scripts/01_microclimate/lib.R")
source("scripts/02_model/config/paths.R")
source("scripts/00_data_conversion/helper_functions.R")

mycredentials <- readRDS(file.path(BASE_DIR, "credentials.rds"))
cds_row <- mycredentials[mycredentials$Site == "CDS", ]
ecmwfr::wf_set_key(key = cds_row$password, user = cds_row$username)
sites <- make_sites(OBSERVATIONS_CSV, pad = 0.15)
site  <- sites[sites$Site == TARGET_SITE, ]
site  <- split(site, seq_len(nrow(site)))[[1]]

ee$Initialize(project = "ee-lizethestevezt")

cat(sprintf("[%s] Setup starting for %s\n", format(Sys.time(), "%H:%M:%S"), TARGET_SITE))
t_setup0 <- Sys.time()

site_dir        <- file.path(RAW_DIR,  site$Site)
site_dtm_dir    <- file.path(site_dir, "dtm")
site_soil_dir   <- file.path(site_dir, "soil")
site_era5_dir   <- file.path(site_dir, "era5")
site_alb_dir    <- file.path(site_dir, "albedo")
site_lai_dir    <- file.path(site_dir, "lai")
site_lcover_dir <- file.path(site_dir, "landcover")
site_model_path <- file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", site$Site))

raster <- terra::rast(nrows = 2, ncols = 2,
  xmin = site$lon_min, xmax = site$lon_max,
  ymin = site$lat_min, ymax = site$lat_max, crs = "EPSG:4326")
terra::values(raster) <- 1

tme_obs_end <- as.POSIXlt(site$tme_end, tz = "UTC")
tme_to <- as.POSIXlt(format(
  seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
      by = "month", length.out = 2L)[2L] - 3600L,
  "%Y-%m-%d %H:00:00"), tz = "UTC")
tme_from <- as.POSIXlt(format(
  seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
      by = "-1 month", length.out = 12L)[12L],
  "%Y-%m-01 00:00:00"), tz = "UTC")
tme <- as.POSIXlt(seq(from = tme_from, to = tme_to, by = "hour"), tz = "UTC")
site$tme_start <- tme_from
site$tme_end   <- tme_to

weatherdata    <- get_weather(site = site, credentials = mycredentials,
                              r = raster, tme = tme, dir = site_era5_dir, output = "grid")
dtmdata        <- get_dtm(r = raster, dir = site_dtm_dir, mask = FALSE)
landcoverdata  <- get_landcover(site = site, r = raster, out_dir = site_lcover_dir)
laidata        <- get_lai(r = raster, tme = tme, pathout = site_lai_dir, credentials = mycredentials)
albedodata     <- get_albedo(r = raster, tme = tme, pathout = site_alb_dir, credentials = mycredentials)
landcover_dtm  <- terra::resample(landcoverdata, dtmdata, method = "near")
refldata       <- get_reflectance(lai = laidata, alb = albedodata, landcover = landcover_dtm,
                                  cachefile = file.path(site_dir, "reflectance.rds"))
vegetationdata <- get_vegetation(r = raster, lcover = landcover_dtm, lai = laidata,
                                 refldata = refldata, dir = site_dir, site_name = site$Site)
soildata       <- get_soil(r = raster, dir = site_soil_dir,
                           landcover = landcover_dtm, refldata = refldata)

stopifnot(file.exists(site_model_path))
model <- readRDS(site_model_path)
n_before <- length(model)
valid_mp <- Filter(function(x) inherits(x, "micropoint"), model)
model <- lapply(model, function(x) if (inherits(x, "micropoint")) x else valid_mp[[1]])
cat(sprintf("Valid grid cells: %d / %d\n", length(valid_mp), n_before))

era5_template <- terra::rast(weatherdata[[1]])[[1]]
dtmc          <- terra::resample(dtmdata, era5_template, method = "bilinear")

t_setup1 <- Sys.time()
cat(sprintf("[%s] Setup done in %.1f s\n", format(Sys.time(), "%H:%M:%S"),
            as.numeric(difftime(t_setup1, t_setup0, units = "secs"))))

OUT_VARS <- c(
  Tz = TRUE, tleaf = FALSE, relhum = TRUE, soilm = FALSE, windspeed = TRUE,
  Rdirdown = TRUE, Rdifdown = TRUE, Rlwdown = FALSE, Rswup = FALSE, Rlwup = FALSE)

cat(sprintf("[%s] Starting SINGLE runmicro() call at height=%.2f m (full year, 1 core, no mclapply)\n",
            format(Sys.time(), "%H:%M:%S"), TEST_HEIGHT))
t0 <- Sys.time()
mout <- microclimf::runmicro(micropoint = model, reqhgt = TEST_HEIGHT,
  vegp = vegetationdata, soilc = soildata, dtm = dtmdata, dtmc = dtmc,
  altcorrect = 1, method = "R", out = unname(OUT_VARS))
t1 <- Sys.time()
elapsed_s <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("[%s] runmicro() done. Elapsed: %.1f s (%.2f min). n timesteps: %d\n",
            format(Sys.time(), "%H:%M:%S"), elapsed_s, elapsed_s / 60, length(mout$tme)))

out_list <- list(tme = mout$tme, Tz = mout$Tz, relhum = mout$relhum,
                  windspeed = mout$windspeed, Rdirdown = mout$Rdirdown, Rdifdown = mout$Rdifdown)
cat(sprintf("In-memory object size (5 vars, this height): %.2f GB\n",
            as.numeric(object.size(out_list)) / 1024^3))

diag_path <- file.path(PROCESSED_DIR, sprintf("diag_single_height_%s.rds", TARGET_SITE))
saveRDS(list(site = TARGET_SITE, height = TEST_HEIGHT, elapsed_s = elapsed_s,
             obj_size_gb = as.numeric(object.size(out_list)) / 1024^3), diag_path)
cat(sprintf("Diagnostic summary saved to %s\n", diag_path))
