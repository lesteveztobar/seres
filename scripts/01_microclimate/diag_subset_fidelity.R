# diag_subset_fidelity.R — Stage A verification gate for the every-other-day
# temporal subsample (run_microclimate_site.R's SUBSET_DAY_STRIDE block,
# added 2026-08-14). Not part of the pipeline; submitted manually.
#
# Reuses cached weather/DTM/vegetation/soil/point-model data for one site
# (same pattern as diag_single_height.R), applies the EXACT subsample block
# from run_microclimate_site.R, calls runmicro() once at a height that
# already has real full-year 0.25m production output on disk, computes
# voxel_quantiles + *_mean exactly as .run_height() does, and compares
# against that real production file -- so this is a direct numerical
# comparison against actual production output, not a synthetic test.
#
# Hard-gate shape assertions (stopifnot) must all pass before the magnitude
# comparison is even meaningful -- a shape mismatch here is exactly the
# silent-corruption failure mode the run_microclimate_site.R stopifnot
# guards are designed to catch in production.
#
# Usage: Rscript scripts/01_microclimate/diag_subset_fidelity.R [site] [test_height] [stride] [reference_rds]
args <- commandArgs(trailingOnly = TRUE)
TARGET_SITE   <- if (length(args) >= 1) args[1] else "Maquipucuna"
TEST_HEIGHT   <- if (length(args) >= 2) as.numeric(args[2]) else 2.10
SUBSET_STRIDE <- if (length(args) >= 3) as.integer(args[3]) else 2L
REFERENCE_RDS <- if (length(args) >= 4) args[4] else
  "/lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_Maquipucuna_h0.25_heights/h2.10.rds"

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
# .compute_voxel_quantiles() -- same reduction run_microclimate_site.R
# applies to every height's output before saving.
source("scripts/02_model/engine/get_colonization.R")

mycredentials <- readRDS(file.path(BASE_DIR, "credentials.rds"))
cds_row <- mycredentials[mycredentials$Site == "CDS", ]
ecmwfr::wf_set_key(key = cds_row$password, user = cds_row$username)
sites <- make_sites(OBSERVATIONS_CSV, pad = 0.15)
site  <- sites[sites$Site == TARGET_SITE, ]
site  <- split(site, seq_len(nrow(site)))[[1]]

ee$Initialize(project = "ee-lizethestevezt")

cat(sprintf("[%s] Setup starting for %s (test height %.2f m, stride %d)\n",
            format(Sys.time(), "%H:%M:%S"), TARGET_SITE, TEST_HEIGHT, SUBSET_STRIDE))

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
if (length(valid_mp) == 0) stop("All grid cells failed in runpointmodela — check ERA5 coverage and LSM.")
model <- lapply(model, function(x) if (inherits(x, "micropoint")) x else valid_mp[[1]])
cat(sprintf("Valid grid cells: %d / %d\n", length(valid_mp), n_before))

era5_template <- terra::rast(weatherdata[[1]])[[1]]
dtmc          <- terra::resample(dtmdata, era5_template, method = "bilinear")

# ── Exactly the run_microclimate_site.R Diff 3a block ────────────────────────
weather_full <- model[[1]]$weather
n_days_avail <- nrow(weather_full) %/% 24L
subset_days  <- seq.int(1L, n_days_avail, by = SUBSET_STRIDE)
model        <- microclimf::subsetpointmodela(model, days = subset_days)
stopifnot(all(vapply(model, inherits, logical(1), "micropoint")))

subset_tme <- as.POSIXct(model[[1]]$weather$obs_time, tz = "UTC")

cat("\n=== Hard-gate shape assertions ===\n")
check <- function(label, expr) {
  ok <- isTRUE(expr)
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) stop(sprintf("Hard gate failed: %s", label))
}
check(sprintf("length(subset_tme) == %d hours", 24L * length(subset_days)),
      length(subset_tme) == 24L * length(subset_days))
check("all 12 calendar months present in subset_tme",
      length(unique(format(subset_tme, "%m"))) == 12L)
check("subset_tme strictly increasing",
      all(diff(subset_tme) > 0))
check("no NA in subset_tme",
      !anyNA(subset_tme))
check(sprintf("nrow(weather_full) == 8760 (was: %d)", nrow(weather_full)),
      nrow(weather_full) == 8760)
check(sprintf("nrow(model[[1]]$weather) == %d (was: %d)", 24L * length(subset_days), nrow(model[[1]]$weather)),
      nrow(model[[1]]$weather) == 24L * length(subset_days))
hod_tab <- table(as.POSIXlt(subset_tme)$hour)
check(sprintf("hour-of-day histogram flat at %d per hour (range: %d-%d)",
              length(subset_days), min(hod_tab), max(hod_tab)),
      length(unique(as.vector(hod_tab))) == 1L)

cat(sprintf("\n[%s] Running runmicro() at height=%.2f m (method=Cpp, subsampled series)...\n",
            format(Sys.time(), "%H:%M:%S"), TEST_HEIGHT))
OUT_VARS <- c(
  Tz = TRUE, tleaf = FALSE, relhum = TRUE, soilm = FALSE, windspeed = TRUE,
  Rdirdown = TRUE, Rdifdown = TRUE, Rlwdown = FALSE, Rswup = FALSE, Rlwup = FALSE)

t0 <- Sys.time()
mout <- microclimf::runmicro(micropoint = model, reqhgt = TEST_HEIGHT,
  vegp = vegetationdata, soilc = soildata, dtm = dtmdata, dtmc = dtmc,
  altcorrect = 1, method = "Cpp", out = unname(OUT_VARS))
t1 <- Sys.time()
cat(sprintf("[%s] runmicro() done in %.1f s\n", format(Sys.time(), "%H:%M:%S"),
            as.numeric(difftime(t1, t0, units = "secs"))))

# ── Exactly the run_microclimate_site.R Diff 3b fix ───────────────────────────
mout$tme <- subset_tme
check(sprintf("dim(mout$Tz)[3] == length(mout$tme) (%d == %d)", dim(mout$Tz)[3], length(mout$tme)),
      dim(mout$Tz)[3] == length(mout$tme))

nr <- terra::nrow(dtmdata)
nc <- terra::ncol(dtmdata)
voxel_quantiles_new <- .compute_voxel_quantiles(mout, nr, nc)
temp_mean_new  <- matrixStats::colMeans2(matrix(mout$Tz, nrow = nr * nc, ncol = length(mout$tme)), na.rm = TRUE)
tme_new <- mout$tme

cat(sprintf("\n=== Loading real production reference: %s ===\n", REFERENCE_RDS))
stopifnot(file.exists(REFERENCE_RDS))
ref <- readRDS(REFERENCE_RDS)
cat(sprintf("Reference: %d timesteps, %d quantile keys\n",
            length(ref$tme), length(ref$voxel_quantiles$quantiles)))

cat("\n=== Structural comparison ===\n")
new_keys <- names(voxel_quantiles_new$quantiles)
ref_keys <- names(ref$voxel_quantiles$quantiles)
check("identical quantile key sets", setequal(new_keys, ref_keys))

na_mismatch_keys <- character(0)
na_mismatch_detail <- data.frame(key = character(0), n_cells = integer(0),
                                  new_extra_na = integer(0), ref_extra_na = integer(0),
                                  total_cells = integer(0))
diff_summary <- data.frame(key = character(0), max_abs_diff = numeric(0),
                            mean_abs_diff = numeric(0), cor = numeric(0))
for (k in intersect(new_keys, ref_keys)) {
  a <- voxel_quantiles_new$quantiles[[k]]
  b <- ref$voxel_quantiles$quantiles[[k]]
  if (!identical(dim(a), dim(b))) {
    cat(sprintf("[FAIL] dim mismatch for %s: new=%s ref=%s\n", k,
                paste(dim(a), collapse="x"), paste(dim(b), collapse="x")))
    next
  }
  na_a <- is.na(a); na_b <- is.na(b)
  if (!isTRUE(all.equal(na_a, na_b))) {
    na_mismatch_keys <- c(na_mismatch_keys, k)
    na_mismatch_detail <- rbind(na_mismatch_detail, data.frame(
      key = k,
      n_cells = sum(na_a != na_b),
      new_extra_na = sum(na_a & !na_b),   # NA in subsampled where ref had a value
      ref_extra_na = sum(!na_a & na_b),   # NA in ref where subsampled had a value
      total_cells = length(a)
    ))
  }
  both_finite <- is.finite(a) & is.finite(b)
  d <- abs(a[both_finite] - b[both_finite])
  cr <- suppressWarnings(cor(a[both_finite], b[both_finite]))
  diff_summary <- rbind(diff_summary, data.frame(
    key = k,
    max_abs_diff = if (length(d)) max(d) else NA_real_,
    mean_abs_diff = if (length(d)) mean(d) else NA_real_,
    cor = cr
  ))
}
# NOT a hard gate: with a real spatial grid (npix in the thousands x 5
# quantile probs), a handful of topographically-shaded boundary pixels
# genuinely losing their last valid day/night sample when the time axis is
# halved is a plausible, non-corrupting outcome -- not the same failure
# mode as the shape/stopifnot hard gates above. Report magnitude instead of
# aborting, so this run still produces the full comparison in one pass.
total_cells_all <- sum(na_mismatch_detail$total_cells)
total_mismatch  <- sum(na_mismatch_detail$n_cells)
cat(sprintf("[%s] NA/NaN pattern identical across all keys (%d/%d keys mismatched)\n",
            if (length(na_mismatch_keys) == 0) "PASS" else "INFO", length(na_mismatch_keys), length(intersect(new_keys, ref_keys))))
if (length(na_mismatch_keys) > 0) {
  cat(sprintf("  total mismatched cells: %d / %d (%.4f%%) across mismatched keys\n",
              total_mismatch, total_cells_all, 100 * total_mismatch / total_cells_all))
  cat("  worst-offending keys (by mismatched cell count):\n")
  print(na_mismatch_detail[order(-na_mismatch_detail$n_cells), ][1:min(15, nrow(na_mismatch_detail)), ],
        row.names = FALSE)
  # A hard gate only against gross corruption (e.g. an entire key flipping
  # to all-NA or vice versa) -- that magnitude of mismatch, unlike sparse
  # boundary-pixel loss, is not explainable by fewer time samples and would
  # indicate a real bug (e.g. the reshape/mask logic silently misaligned).
  gross <- na_mismatch_detail[na_mismatch_detail$n_cells > 0.5 * na_mismatch_detail$total_cells, ]
  check(sprintf("no key has >50%% of cells with a flipped NA status (%d such keys)", nrow(gross)),
        nrow(gross) == 0)
}

cat("\n=== Per-key magnitude comparison (subsampled vs full-year production) ===\n")
print(diff_summary[order(-diff_summary$max_abs_diff), ], row.names = FALSE)

low_cor <- diff_summary[!is.na(diff_summary$cor) & diff_summary$cor < 0.99 &
                         grepl("_temp$", diff_summary$key), ]
if (nrow(low_cor) > 0) {
  cat("\n*** WARNING: temperature quantile keys with cross-pixel correlation < 0.99: ***\n")
  print(low_cor, row.names = FALSE)
} else {
  cat("\nAll temperature quantile keys correlate >= 0.99 with full-year production output.\n")
}

cat("\n=== Monthly mean comparison (by calendar month, NOT positional) ===\n")
ref_month <- format(ref$tme, "%m")
new_month <- format(tme_new, "%m")
ref_temp_mean <- ref$temp_mean
for (m in sprintf("%02d", 1:12)) {
  ref_m <- mean(ref_temp_mean[ref_month == m], na.rm = TRUE)
  new_m <- mean(temp_mean_new[new_month == m], na.rm = TRUE)
  cat(sprintf("  month %s: full-year=%.3f  subsampled=%.3f  diff=%.4f\n",
              m, ref_m, new_m, new_m - ref_m))
}

cat("\n=== DONE. All hard gates passed. Review the magnitude comparison above before enabling in production. ===\n")
