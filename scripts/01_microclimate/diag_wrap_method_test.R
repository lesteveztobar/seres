# diag_wrap_method_test.R — copy of diag_wrap_collision.R, parameterized by
# runmicro()'s `method` argument ("R" | "Cpp") so the same wrap()/unwrap()+
# mclapply fork path can be timed and correctness-checked under both.
# Not part of the pipeline; submitted manually. Does NOT touch
# run_microclimate_site.R or diag_wrap_collision.R.
#
# Unlike diag_wrap_collision.R (which only saved tme+n per height, since it
# was built purely to check for try-errors), this script saves the FULL
# per-pixel/per-hour output arrays (Tz, relhum, windspeed, Rdirdown, Rdifdown)
# for every height tested, so method="R" vs method="Cpp" outputs for the same
# site/height can be compared value-for-value afterward.
#
# MODE "serial": mc.cores=1 -- see diag_wrap_collision.R's header for why.
# MODE "parallel": mc.cores=N_CORES_TEST, N_HEIGHTS_TEST heights.
args <- commandArgs(trailingOnly = TRUE)
TARGET_SITE     <- if (length(args) >= 1) args[1] else "LaElenita"
MODE            <- if (length(args) >= 2) args[2] else "serial"      # "serial" | "parallel"
N_CORES_TEST    <- if (length(args) >= 3) as.integer(args[3]) else 1L
N_HEIGHTS_TEST  <- if (length(args) >= 4) as.integer(args[4]) else 1L
HEIGHT_START    <- if (length(args) >= 5) as.numeric(args[5]) else 5.10
HEIGHT_STEP_TST <- if (length(args) >= 6) as.numeric(args[6]) else 0.10
METHOD_TEST     <- if (length(args) >= 7) args[7] else "Cpp"          # "R" | "Cpp"

source("scripts/02_model/config/patches.R")
library(rgee)
library(readr)
library(mcera5)
library(microclimf)
library(microclimdata)
library(terra)
library(luna)
library(parallel)

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

cat(sprintf("[%s] diag_wrap_method_test.R — site=%s mode=%s method=%s cores=%d n_heights=%d\n",
            format(Sys.time(), "%H:%M:%S"), TARGET_SITE, MODE, METHOD_TEST, N_CORES_TEST, N_HEIGHTS_TEST))

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

OUT_VARS <- c(
  Tz = TRUE, tleaf = FALSE, relhum = TRUE, soilm = FALSE, windspeed = TRUE,
  Rdirdown = TRUE, Rdifdown = TRUE, Rlwdown = FALSE, Rswup = FALSE, Rlwup = FALSE)

# ── EXACTLY the real loop's wrap step (run_microclimate_site.R L242-245) ──────
era5_template <- terra::rast(weatherdata[[1]])[[1]]
dtmc          <- terra::resample(dtmdata, era5_template, method = "bilinear")
dtmdata_w     <- terra::wrap(dtmdata)
dtmc_w        <- terra::wrap(dtmc)

heights_test <- seq(HEIGHT_START, by = HEIGHT_STEP_TST, length.out = N_HEIGHTS_TEST)
tag <- sprintf("%s_%s", tolower(METHOD_TEST), MODE)
out_dir <- file.path(PROCESSED_DIR, "diag_wrap_heights")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(out_dir, pattern = sprintf("^%s_", tag), full.names = TRUE))

.log_file <- file.path(LOGS_DIR, sprintf("diag_wrap_%s_%s.log", tag, Sys.getenv("SLURM_JOB_ID", "local")))
.wlog <- wlog(.log_file, with_pid = TRUE)  # from lib.R

t0 <- Sys.time()
results <- parallel::mclapply(seq_along(heights_test), function(i) {
  h <- heights_test[i]
  h_out <- file.path(out_dir, sprintf("%s_h%.2f.rds", tag, h))
  .wlog(sprintf("[%d/%d] %.2f m — starting unwrap + runmicro(method=%s)...", i, N_HEIGHTS_TEST, h, METHOD_TEST))
  dtm_  <- tryCatch(terra::unwrap(dtmdata_w), error = function(e) e)
  if (inherits(dtm_, "error")) { .wlog(sprintf("[%d/%d] %.2f m — UNWRAP(dtmdata_w) FAILED: %s", i, N_HEIGHTS_TEST, h, conditionMessage(dtm_))); return(dtm_) }
  dtmc_ <- tryCatch(terra::unwrap(dtmc_w), error = function(e) e)
  if (inherits(dtmc_, "error")) { .wlog(sprintf("[%d/%d] %.2f m — UNWRAP(dtmc_w) FAILED: %s", i, N_HEIGHTS_TEST, h, conditionMessage(dtmc_))); return(dtmc_) }
  t_h0 <- Sys.time()
  mout <- tryCatch(
    microclimf::runmicro(micropoint = model, reqhgt = h,
      vegp = vegetationdata, soilc = soildata, dtm = dtm_, dtmc = dtmc_,
      altcorrect = 1, method = METHOD_TEST, out = unname(OUT_VARS)),
    error = function(e) e)
  if (inherits(mout, "error")) { .wlog(sprintf("[%d/%d] %.2f m — RUNMICRO FAILED: %s", i, N_HEIGHTS_TEST, h, conditionMessage(mout))); return(mout) }
  t_h1 <- Sys.time()
  h_elapsed_s <- as.numeric(difftime(t_h1, t_h0, units = "secs"))
  # Full field values saved (not just tme+n) -- needed for the R-vs-Cpp
  # correctness comparison this script exists for.
  saveRDS(list(tme = mout$tme, n = length(mout$tme), method = METHOD_TEST,
               elapsed_s = h_elapsed_s,
               Tz = mout$Tz, relhum = mout$relhum, windspeed = mout$windspeed,
               Rdirdown = mout$Rdirdown, Rdifdown = mout$Rdifdown),
          h_out)
  .wlog(sprintf("[%d/%d] %.2f m — done. %d timesteps, %.1f s. Saved: %s",
                i, N_HEIGHTS_TEST, h, length(mout$tme), h_elapsed_s, h_out))
  invisible(NULL)
}, mc.cores = N_CORES_TEST)
t1 <- Sys.time()

n_ok  <- sum(vapply(results, function(x) is.null(x), logical(1)))
n_err <- N_HEIGHTS_TEST - n_ok
cat(sprintf("[%s] Done in %.1f s. %d/%d heights OK, %d errored. method=%s\n",
            format(Sys.time(), "%H:%M:%S"),
            as.numeric(difftime(t1, t0, units = "secs")), n_ok, N_HEIGHTS_TEST, n_err, METHOD_TEST))
cat("=== mclapply return value (auto-printed for transcript) ===\n")
print(results)

summary_path <- file.path(PROCESSED_DIR,
  sprintf("diag_wrap_method_%s_%s_%s.rds", tag, TARGET_SITE, Sys.getenv("SLURM_JOB_ID", "local")))
saveRDS(list(mode = MODE, method = METHOD_TEST, n_cores = N_CORES_TEST, n_heights = N_HEIGHTS_TEST,
             n_ok = n_ok, n_err = n_err, elapsed_s = as.numeric(difftime(t1, t0, units = "secs")),
             results = results), summary_path)
cat(sprintf("Summary saved to %s\n", summary_path))
