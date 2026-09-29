# run_microclimate_site.R
# Runs the full microclimate pipeline for ONE site, specified as a command-line arg.
# Called by hpc_job.sh: Rscript run_microclimate_site.R Maquipucuna
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript run_microclimate_site.R <SiteName>[,<SiteName>,...] [n_months] [height_step]")
# 2026-09-06 (v7 3D re-run, Phase A): a comma-separated TARGET_SITE list
# means "these sites share one raster domain" -- the Mindo-cluster case
# (MindoMirador, MindoTarabita, Saloya, LaElenita: ~7-8km apart, existing
# rasters already overlap 55-95%). One microclimf run, one set of forcing
# data, one DTM; each listed site still gets its own manifest (below) and
# its own per-pixel footprint within the shared output, so downstream code
# addresses them exactly as it would separate single-site runs. This also
# gives the 4 cluster sites identical forcing/model run by construction --
# a deliberate secondary benefit (see the pricing report, Phase A3).
CLUSTER_SITES <- strsplit(args[1], ",", fixed = TRUE)[[1]]
IS_CLUSTER     <- length(CLUSTER_SITES) > 1L
TARGET_SITE   <- if (IS_CLUSTER) paste0("MindoCluster_", paste(CLUSTER_SITES, collapse = "-")) else CLUSTER_SITES[1]
N_MONTHS    <- if (length(args) >= 2) as.integer(args[2]) else 12L
# Vertical spacing (m) between height tiers. Default 0.1 matches production
# sites. A coarser step (e.g. 0.5) is for the height-resolution efficiency
# diagnostic — it reuses the cached weather/DTM/point-model data (those don't
# depend on height spacing) and only redoes the per-height loop below, so
# testing a coarser resolution is proportionally cheaper, not free.
HEIGHT_STEP <- if (length(args) >= 3) as.numeric(args[3]) else 0.1

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

# 2026-08-10: keep terra's scratch rasters off node-local /tmp -- it filled
# up during Maquipucuna's first production run (76 heights x 4 concurrent
# mclapply workers, each doing several wrap/unwrap/resample/mask ops, none
# of it ever cleaned up), causing every worker to fail simultaneously with
# "[resample] cannot write file" / "[mask] cannot read from /tmp/..." once
# node-local /tmp filled (logs/run_microclimate_26922710.err, 2026-08-09/10)
# -- silently truncating that site's manifest to 16/76 heights. microenv_
# array.sh already exports TMPDIR under the much larger Lustre scratch
# workspace before this script runs; this is a belt-and-suspenders guard
# in case some code path doesn't reliably inherit TMPDIR through the
# mclapply() fork below.
terraOptions(tempdir = Sys.getenv("TMPDIR", unset = tempdir()))

source("scripts/01_microclimate/lib.R")
source("scripts/02_model/config/paths.R")
source("scripts/data_prep/helper_functions.R")
# .compute_voxel_quantiles() -- the per-pixel/month/daypart quantile
# reduction applied to each height's output below, before it's written to
# scratch. get_colonization.R has no top-level library()/source() calls of
# its own, so this is a lightweight addition.
source("scripts/02_model/engine/get_colonization.R")

mycredentials <- readRDS(file.path(BASE_DIR, "credentials.rds"))
cds_row <- mycredentials[mycredentials$Site == "CDS", ]
ecmwfr::wf_set_key(key = cds_row$password, user = cds_row$username)
if (IS_CLUSTER) {
  # Same summary make_sites() computes per site, but grouped across the
  # whole cluster instead of by individual Area_or_Site -- i.e. one padded
  # bounding box, one time window, spanning every listed site's
  # observations. A3's coverage check (pricing report) confirmed this
  # combined bbox, padded by the same 0.15 degrees make_sites() already
  # uses, gives >=9km margin beyond the outermost site's own footprint --
  # more than any single site's own domain currently gets, since the
  # cluster's raw bbox is itself larger than any one site's.
  df <- readr::read_csv(OBSERVATIONS_CSV, na = c("", "NA", "N/A"),
                        col_types = COMBINED_COL_TYPES, show_col_types = FALSE) |>
    dplyr::filter(!is.na(Source), !is.na(Area_or_Site), Area_or_Site %in% CLUSTER_SITES) |>
    dplyr::filter(!is.na(lat), !is.na(lon), !is.na(datetime), !is.na(Height_m)) |>
    dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon),
                  datetime = as.POSIXlt(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
  site <- data.frame(
    Site        = TARGET_SITE,
    lat_min     = min(df$lat) - 0.15, lat_max = max(df$lat) + 0.15,
    lon_min     = min(df$lon) - 0.15, lon_max = max(df$lon) + 0.15,
    tme_start   = min(df$datetime),   tme_end = max(df$datetime),
    hObs_min    = min(df$Height_m, na.rm = TRUE), hObs_max = max(df$Height_m, na.rm = TRUE),
    hCanopy_max = suppressWarnings({ v <- max(df$CanopyHeight_m, na.rm = TRUE); if (is.finite(v)) v else NA_real_ })
  )
  # Per-site (unpadded) observation bbox, for .site_footprint_pixels() below
  # -- computed once here from the same filtered df, independent of the
  # shared raster's own (padded) extent.
  CLUSTER_SITE_BBOX <- lapply(CLUSTER_SITES, function(s) {
    d <- df[df$Area_or_Site == s, ]
    list(lat_range = range(d$lat), lon_range = range(d$lon))
  })
  names(CLUSTER_SITE_BBOX) <- CLUSTER_SITES
  log_needed <- TRUE
} else {
  sites         <- make_sites(OBSERVATIONS_CSV, pad = 0.15)
  site          <- sites[sites$Site == TARGET_SITE, ]
  if (nrow(site) == 0) stop(sprintf("Site '%s' not found in sites table.", TARGET_SITE))
  site <- split(site, seq_len(nrow(site)))[[1]]

  # 2026-09-10: ERA5 grid-snap REMOVED. It shifted the model domain (up to
  # 11.5 km) to guarantee >=2 ERA5 cell centres, but the terrain/vegetation
  # rasters (get_dtm() etc.) load from a filename-keyed cache without
  # checking the cached extent, so terrain stayed at the pre-snap bbox
  # while ERA5 was pulled for the shifted domain -- ERA5 cells the snap
  # added fell off the cached DTM, their point models came back NA, and
  # microclimf's spatial run collapsed. All three snapped sites (Mindo
  # Mirador 0.000, MindoTarabita 0.513, Saloya 0.550 finite temp fraction)
  # were degraded; none of the three unsnapped sites were. The snap's own
  # justification is also gone: the independent 12-timestep monthly
  # validity check finds every candidate ERA5 cell valid at every site, so
  # the "1 valid cell" problem the snap solved was an artifact of
  # microclimf::runpointmodela()'s hour-1-only test, not real.
  #
  # Domain is now just the raw-observation bounding box + the fixed 0.15deg
  # buffer from make_sites() -- deterministic, reproducible, centred on the
  # site. A site whose domain falls within a single ERA5 cell is fine: all
  # horizontal microclimatic variation in the model derives from terrain
  # and canopy structure (microclimf), not the weather forcing.
  # `snap` is kept as a null-op record so the manifest field
  # (.era5_grid_snap) and any code testing exists(snap) still resolve.
  snap <- list(lat_min = site$lat_min, lat_max = site$lat_max,
               lon_min = site$lon_min, lon_max = site$lon_max,
               era5_grid_snap = list(axis = "none", shift_deg = 0, shift_km = 0,
                                     removed = TRUE))
  message(sprintf("Domain for %s: raw obs bbox + 0.15deg buffer (no ERA5 grid snap). lat [%.4f, %.4f] lon [%.4f, %.4f]",
                  TARGET_SITE, site$lat_min, site$lat_max, site$lon_min, site$lon_max))
  # CANOPY_GRIDORIGIN_OFFSET_M: retained diagnostic hook only (grid-origin
  # sensitivity test); unset in every normal run.
  test_offset_m <- as.numeric(Sys.getenv("CANOPY_GRIDORIGIN_OFFSET_M", unset = "0"))
  if (test_offset_m != 0) {
    d_lat <- test_offset_m / 111320
    d_lon <- test_offset_m / (111320 * cos(mean(c(site$lat_min, site$lat_max)) * pi / 180))
    site$lat_min <- site$lat_min + d_lat; site$lat_max <- site$lat_max + d_lat
    site$lon_min <- site$lon_min + d_lon; site$lon_max <- site$lon_max + d_lon
    message(sprintf("GRID-ORIGIN SENSITIVITY TEST: additional offset of %.1fm applied (lat+%.6fdeg, lon+%.6fdeg)",
                    test_offset_m, d_lat, d_lon))
  }
  # site$lat_min/lat_max etc. are already padded by 0.15deg (make_sites()) --
  # that's the RASTER domain, not the landscape footprint. Recompute the
  # site's own unpadded observation bbox directly, so solo-site and
  # cluster-site runs feed .site_footprint_pixels() the same kind of input.
  df1 <- readr::read_csv(OBSERVATIONS_CSV, na = c("", "NA", "N/A"),
                         col_types = COMBINED_COL_TYPES, show_col_types = FALSE) |>
    dplyr::filter(!is.na(Source), Area_or_Site == TARGET_SITE,
                  !is.na(lat), !is.na(lon), !is.na(datetime), !is.na(Height_m)) |>
    dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  CLUSTER_SITE_BBOX <- setNames(list(list(lat_range = range(df1$lat), lon_range = range(df1$lon))), TARGET_SITE)
}

ee$Initialize(project = "ee-lizethestevezt")

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR,
  sprintf("microclim_%s_%s.log", TARGET_SITE, format(Sys.time(), "%Y%m%d_%H%M%S")))
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(log_file = log_file)
log_msg(sprintf("run_microclimate_site.R started for %s", TARGET_SITE))

# ── Directories ───────────────────────────────────────────────────────────────
site_dir        <- file.path(RAW_DIR,  site$Site)
site_dtm_dir    <- file.path(site_dir, "dtm")
site_soil_dir   <- file.path(site_dir, "soil")
site_era5_dir   <- file.path(site_dir, "era5")
site_alb_dir    <- file.path(site_dir, "albedo")
site_lai_dir    <- file.path(site_dir, "lai")
site_lcover_dir <- file.path(site_dir, "landcover")
for (d in c(site_dir, site_dtm_dir, site_soil_dir,
            site_era5_dir, site_alb_dir, site_lcover_dir, site_lai_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Resolution-suffixed file names only kick in for a non-default height step,
# so production (0.1m) runs keep their original, unsuffixed file names and
# coarser diagnostic runs never collide with or overwrite them.
res_suffix      <- if (HEIGHT_STEP != 0.1) sprintf("_h%.2f", HEIGHT_STEP) else ""
if (exists("test_offset_m") && test_offset_m != 0) {
  res_suffix <- sprintf("%s_gridoffset%.0fm", res_suffix, test_offset_m)
}
site_env_path   <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site$Site, res_suffix))
site_model_path <- file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", site$Site))

# Per-height temp files go to Lustre scratch (set via CANOPY_SCRATCH in hpc_job.sh).
# Falls back to PROCESSED_DIR for local/interactive runs.
scratch_base <- Sys.getenv("CANOPY_SCRATCH", unset = "")
scratch_base <- if (nchar(scratch_base) > 0 && dir.exists(scratch_base)) scratch_base else PROCESSED_DIR

if (file.exists(site_env_path)) {
  log_msg(sprintf("Microenv already exists for %s — nothing to do.", site$Site))
  quit(status = 0)
}

# ── Spatial / temporal setup ──────────────────────────────────────────────────
raster <- terra::rast(
  nrows = 2, ncols = 2,
  xmin  = site$lon_min, xmax = site$lon_max,
  ymin  = site$lat_min, ymax = site$lat_max,
  crs   = "EPSG:4326"
)
terra::values(raster) <- 1

tme_obs_end <- as.POSIXlt(site$tme_end, tz = "UTC")
tme_to <- as.POSIXlt(format(
  seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
      by = "month", length.out = 2L)[2L] - 3600L,
  "%Y-%m-%d %H:00:00"), tz = "UTC")
tme_from <- as.POSIXlt(format(
  seq(as.POSIXlt(format(tme_obs_end, "%Y-%m-01"), tz = "UTC"),
      by = "-1 month", length.out = N_MONTHS)[N_MONTHS],
  "%Y-%m-01 00:00:00"), tz = "UTC")
tme <- as.POSIXlt(seq(from = tme_from, to = tme_to, by = "hour"), tz = "UTC")
site$tme_start <- tme_from
site$tme_end   <- tme_to

# ── 1. Data acquisition ───────────────────────────────────────────────────────
log_msg("Acquiring weather data...")
weatherdata <- get_weather(site = site, credentials = mycredentials,
                           r = raster, tme = tme, dir = site_era5_dir, output = "grid")

# 2026-09-08: independent, our-own-code validity check -- microclimf's own
# per-cell test (runpointmodela(), quoted in docs/methods_update_report.md)
# only inspects the FIRST hourly timestep (`climdf$temp[1]`), so a cell with
# one missing hour at the very start of the series is discarded even if the
# rest of its year is fine, and (the opposite, more dangerous failure) a
# cell fine at hour 1 but gappy afterward passes silently. This samples one
# timestep per month (12 total) and requires ALL of them finite -- a
# tougher, more representative bar, computed independently of and BEFORE
# runpointmodela() is ever called, so it never depends on or patches that
# package internal. Logged as a cross-check; N_VALID_ERA5_CELLS (below,
# microclimf's own hour-1 count) is still what actually gates which cells
# feed the model, per the decision not to patch runmicro()/runpointmodela().
# get_weather(..., output="grid") returns wrapped (PackedSpatRaster)
# components for safe cross-process return -- unwrap before any terra::
# call, same pattern this script already uses for dtmdata_w/dtmc_w below.
weatherdata_temp <- if (inherits(weatherdata$temp, "PackedSpatRaster")) terra::unwrap(weatherdata$temp) else weatherdata$temp
ERA5_GRID_NROW <- terra::nrow(weatherdata_temp)
ERA5_GRID_NCOL <- terra::ncol(weatherdata_temp)
n_time_avail <- terra::nlyr(weatherdata_temp)
month_sample_idx <- unique(round(seq(1, n_time_avail, length.out = 12)))
temp_monthly_samples <- terra::values(weatherdata_temp[[month_sample_idx]])
valid_monthly <- apply(temp_monthly_samples, 1, function(x) all(is.finite(x)))
N_VALID_ERA5_CELLS_MONTHLY <- sum(valid_monthly)
N_TOTAL_ERA5_CELLS_MONTHLY <- length(valid_monthly)
log_msg(sprintf(
  "Valid grid cells (monthly-sampled, %d timesteps across the year, ALL must be finite): %d / %d",
  length(month_sample_idx), N_VALID_ERA5_CELLS_MONTHLY, N_TOTAL_ERA5_CELLS_MONTHLY))

log_msg("Acquiring DTM...")
dtmdata <- get_dtm(r = raster, dir = site_dtm_dir, mask = FALSE)

log_msg("Acquiring landcover...")
landcoverdata <- get_landcover(site = site, r = raster, out_dir = site_lcover_dir)

log_msg("Acquiring LAI...")
laidata <- get_lai(r = raster, tme = tme, pathout = site_lai_dir, credentials = mycredentials)

log_msg("Acquiring albedo...")
albedodata <- get_albedo(r = raster, tme = tme, pathout = site_alb_dir, credentials = mycredentials)

# ── Terrain-coverage assertion (2026-09-10) ──────────────────────────────
# Every acquired terrain/vegetation raster MUST enclose the requested model
# domain. Before the .cache_covers() guard (lib.R), a filename-keyed cache
# could silently return terrain for an earlier, different domain while ERA5
# was re-pulled for the current one -- MindoMirador's microclimate
# collapsed exactly this way. Fail loudly here, at build time.
.req_ext <- as.vector(terra::ext(raster))
.assert_covers <- function(x, name, hard = FALSE) {
  ce <- as.vector(terra::ext(x)); tol <- 0.003
  m_s <- (ce[3] - .req_ext[3]) * 111.32
  m_n <- (.req_ext[4] - ce[4]) * -111.32
  m_w <- (ce[1] - .req_ext[1]) * 111.32
  m_e <- (.req_ext[2] - ce[2]) * -111.32
  log_msg(sprintf("  terrain coverage %-9s: gap S/N/W/E from domain edge = %+.2f / %+.2f / %+.2f / %+.2f km (negative = raster extends past domain)",
                  name, m_s, m_n, m_w, m_e))
  short <- ce[1] > .req_ext[1] + tol || ce[2] < .req_ext[2] - tol ||
           ce[3] > .req_ext[3] + tol || ce[4] < .req_ext[4] - tol
  if (short) {
    msg <- sprintf("%s raster [%.4f,%.4f,%.4f,%.4f] does NOT enclose the model domain [%.4f,%.4f,%.4f,%.4f] -- stale cache (see .cache_covers()).",
                   name, ce[1], ce[2], ce[3], ce[4], .req_ext[1], .req_ext[2], .req_ext[3], .req_ext[4])
    if (hard) stop(msg) else log_msg(paste("  WARN:", msg))
  }
}
.assert_covers(dtmdata, "DTM", hard = TRUE)
.assert_covers(landcoverdata, "landcover")
.assert_covers(laidata, "LAI")

landcover_dtm <- terra::resample(landcoverdata, dtmdata, method = "near")

log_msg("Computing reflectance...")
refldata <- get_reflectance(lai = laidata, alb = albedodata, landcover = landcover_dtm,
                            cachefile = file.path(site_dir, "reflectance.rds"))

log_msg("Deriving vegetation parameters...")
vegetationdata <- get_vegetation(r = raster, lcover = landcover_dtm, lai = laidata,
                                 refldata = refldata, dir = site_dir, site_name = site$Site)

log_msg("Deriving soil parameters...")
soildata <- get_soil(r = raster, dir = site_soil_dir,
                     landcover = landcover_dtm, refldata = refldata)

# ── 2. Point model ────────────────────────────────────────────────────────────
# Cache is only valid if it's newer than the merged ERA5 file it was built
# from -- otherwise weatherdata (always freshly read above) and model can
# disagree on grid-cell count/structure. This is exactly what happened to
# Saloya on 2026-07-24: a point model cached from an earlier ERA5 download
# got silently reused against a later, re-downloaded Saloya.nc with a
# different valid-cell count, producing "weather[[k]] : subscript out of
# bounds" during the height loop below.
site_era5_path <- file.path(site_era5_dir, sprintf("%s.nc", site$Site))
model_is_stale <- file.exists(site_model_path) && file.exists(site_era5_path) &&
  file.info(site_era5_path)$mtime > file.info(site_model_path)$mtime

if (model_is_stale) {
  log_msg("Cached point model predates the merged ERA5 file -- discarding stale cache.")
}

if (file.exists(site_model_path) && !model_is_stale) {
  log_msg("Loading existing point model...")
  model <- readRDS(site_model_path)
} else {
  log_msg("Running point model...")
  model <- microclimf::runpointmodela(
    climarrayr = weatherdata, tme = tme, reqhgt = 0.05,
    dtm = dtmdata, vegp = vegetationdata, soilc = soildata)
  saveRDS(model, site_model_path)
  log_msg(sprintf("Saved point model to %s", site_model_path))
}

# ── 3. Grid model (parallel heights) ─────────────────────────────────────────

# runpointmodela returns NA for cells where the model fails (e.g. ocean/masked
# cells). runmicro iterates weather[[k]] up to length(dtmc) — the ERA5 cell
# count — so the model list must stay that length.
#
# 2026-09-08: substitution changed from "copy the first valid cell" (scan-
# order in runpointmodela()'s own i,j loop -- row-major, so effectively
# "the northernmost, then westernmost valid cell", with no geographic
# relationship to the masked cell it's filling in for) to "copy the
# NEAREST valid cell" (Euclidean distance in grid-index space -- ERA5 cells
# are on a uniform grid, so index distance is proportional to physical
# distance). The old comment here claimed masked cells' spatial output is
# "masked/unused anyway" -- disproven directly: MindoMirador's 2026-09-07
# re-run (1 valid ERA5 cell) produced a per-pixel array 99.97% NA, i.e. the
# substituted values clearly DO reach and dominate the final per-pixel
# climate output, not just fill in something later discarded. With the
# ERA5 grid-snap (era5_grid_snap(), lib.R) now deliberately choosing a
# domain's shift DIRECTION by a geometric rule (not geography), "first in
# scan order" would have made which cell's weather ends up governing most
# of a domain's masked area depend on which direction the snap happened to
# shift it -- an artifact of the fix itself, not a real forcing difference
# between sites. Nearest-valid substitution removes that dependency: for
# every masked cell, whichever valid cell is geographically closest to IT
# specifically is used, not whichever was scanned first.
n_before  <- length(model)
valid_idx <- which(vapply(model, inherits, logical(1), "micropoint"))
if (length(valid_idx) == 0) stop("All grid cells failed in runpointmodela — check ERA5 coverage and LSM.")
# k = (i-1)*ncol + j (row-major, matching runpointmodela()'s own i,j loop
# order) -> recover each cell's (row, col) grid position.
.k_to_rc <- function(k, ncol) c(row = ((k - 1L) %/% ncol) + 1L, col = ((k - 1L) %% ncol) + 1L)
valid_rc <- t(vapply(valid_idx, .k_to_rc, numeric(2), ncol = ERA5_GRID_NCOL))
model <- lapply(seq_along(model), function(k) {
  if (inherits(model[[k]], "micropoint")) return(model[[k]])
  rc <- .k_to_rc(k, ERA5_GRID_NCOL)
  d2 <- (valid_rc[, "row"] - rc[["row"]])^2 + (valid_rc[, "col"] - rc[["col"]])^2
  model[[valid_idx[which.min(d2)]]]
})
log_msg(sprintf("Valid grid cells: %d / %d (masked cells filled from their own nearest valid cell, not scan-order-first)", length(valid_idx), n_before))
# 2026-09-08: persisted into the manifest below (N_VALID_ERA5_CELLS/
# N_TOTAL_ERA5_CELLS), not just this log line -- sanity_gate.R's new (g)/(h)
# checks need this as a structural fact they can assert on, not something
# that only ever existed as text in a log file. See runpointmodela()
# (microclimf) for exactly what makes a cell "valid": per its own source,
# a cell is valid iff its first ERA5 hourly temperature reading is not NA
# AND the vegetation-height raster (resampled to the ERA5 grid) is not NA
# at that cell -- checked only at the FIRST timestep, not the whole series.
N_VALID_ERA5_CELLS <- length(valid_idx)
N_TOTAL_ERA5_CELLS <- n_before
if (N_VALID_ERA5_CELLS != N_VALID_ERA5_CELLS_MONTHLY) {
  log_msg(sprintf(
    "NOTE: hour-1-only valid count (%d) differs from monthly-sampled count (%d) -- the hour-1 test missed a temporally-gappy cell, or discarded one only briefly missing at hour 1.",
    N_VALID_ERA5_CELLS, N_VALID_ERA5_CELLS_MONTHLY))
}

# ── 2b. Temporal subsample before grid expansion ─────────────────────────────
# The per-height grid expansion is O(pixels x timesteps); the point model
# build above is not. Halving the time axis here (every other day, full
# 24h blocks) halves both runtime and the (npix x ntime) arrays that
# dominate peak RSS, at full spatial and height resolution.
# subsetpointmodela(days=<numeric>) is purely positional (day 1 = this
# model's own first day, whatever the actual calendar date), so it needs
# no adjustment for this site's Apr->Mar window. microclimf's own
# .runmodel2Cpp() detects the shortened `subs` and switches to its
# complete=FALSE path; .sortvegp() re-indexes PAI phenology via `subs`, so
# vegetation seasonality stays aligned to the retained days.
# 2026-08-14: verified via diag_subset_fidelity.R against a real
# full-year production height (Maquipucuna h2.10) before this was turned
# on in production -- see that script for the exact comparison performed.
# Overridable via env var so a pilot run can validate cores/node-exclusion
# in isolation at SUBSET_DAY_STRIDE=1 (i.e. off, full-year) without editing
# this file again for the next, subsample-on production run.
SUBSET_DAY_STRIDE <- as.integer(Sys.getenv("SUBSET_DAY_STRIDE", "2"))

# Preserve the FULL, unsubsetted ERA5 weather record for the manifest
# BEFORE subsetting `model` in place -- get_colonization.R's
# lookup_climate_by_height() needs every hour, not just the retained days.
weather_full <- model[[1]]$weather

n_days_avail <- nrow(weather_full) %/% 24L
subset_days  <- seq.int(1L, n_days_avail, by = SUBSET_DAY_STRIDE)
model        <- microclimf::subsetpointmodela(model, days = subset_days)
stopifnot(all(vapply(model, inherits, logical(1), "micropoint")))

subset_tme <- as.POSIXct(model[[1]]$weather$obs_time, tz = "UTC")
stopifnot(length(subset_tme) == 24L * length(subset_days),
          length(unique(format(subset_tme, "%m"))) == 12L)
log_msg(sprintf(
  "Temporal subsample: every %dth day -> %d days / %d hours (from %d hours).",
  SUBSET_DAY_STRIDE, length(subset_days), length(subset_tme), nrow(weather_full)))

# Only the 5 variables lookup_climate_by_height()/get_clim_voxel() (get_colonization.R)
# actually read -- Tz, relhum, windspeed, Rdirdown, Rdifdown -- are requested
# via `out`, roughly halving storage/serialization (tleaf, soilm, Rlwdown,
# Rswup, Rlwup are computed internally but never written to disk).
OUT_VARS <- c(
  Tz = TRUE, tleaf = FALSE, relhum = TRUE, soilm = FALSE, windspeed = TRUE,
  Rdirdown = TRUE, Rdifdown = TRUE, Rlwdown = FALSE, Rswup = FALSE, Rlwup = FALSE
)

# The model's height ceiling must be the canopy top, not the tallest
# recorded epiphyte observation (site$hObs_max) -- see height_ceiling()
# (lib.R, shared with check_microenv_progress.R) for the preference order
# and rationale (measured CanopyHeight_m > p99 of vhgt.tif > hObs_max).
height_ceiling_m <- height_ceiling(site$hCanopy_max, file.path(site_dir, "vhgt.tif"),
                                   site$hObs_max, log_fn = log_msg)

heights    <- seq(0.1, height_ceiling_m, by = HEIGHT_STEP)
height_dir <- file.path(scratch_base, sprintf("microenv_%s%s_heights", site$Site, res_suffix))
dir.create(height_dir, recursive = TRUE, showWarnings = FALSE)

era5_template <- terra::rast(weatherdata[[1]])[[1]]
dtmc          <- terra::resample(dtmdata, era5_template, method = "bilinear")
dtmdata_w     <- terra::wrap(dtmdata)
dtmc_w        <- terra::wrap(dtmc)

# v7 3D re-run (Phase A1): per-site footprint pixel indices within THIS
# run's raster (shared across all CLUSTER_SITES for a cluster run, or just
# TARGET_SITE itself for a solo run) -- built once here from dtmdata's own
# extent/dims, reused by every height tier's .compute_pixel_means() call
# below. Uses .site_footprint_pixels()/.lonlat_to_pixel() from
# get_colonization.R (already sourced above).
sp_spatial_run <- list(
  ext  = as.vector(terra::ext(dtmdata)),
  nrow = terra::nrow(dtmdata),
  ncol = terra::ncol(dtmdata)
)
site_footprints <- lapply(CLUSTER_SITE_BBOX, function(bb)
  .site_footprint_pixels(bb$lat_range, bb$lon_range, sp_spatial_run))
log_msg(sprintf("Per-pixel-mean footprints (this run's %dx%d raster): %s",
  sp_spatial_run$nrow, sp_spatial_run$ncol,
  paste(sprintf("%s=%d px", names(site_footprints), lengths(site_footprints)), collapse = ", ")))

n_cores   <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
n_heights <- length(heights)
log_msg(sprintf("Launching parallel height loop: %d heights on %d cores...", n_heights, n_cores))

.log_file <- log_file
# with_pid=TRUE (2026-08-10): concurrent workers' lines used to be
# unattributable, and there was no error-visibility at all here -- see
# .run_height()'s tryCatch below, added after Maquipucuna's first
# production run silently dropped heights 17-76 (node-local /tmp filled;
# every worker's runmicro() call failed with no trace anywhere except
# Rscript's auto-printed, never-inspected mclapply() return value).
.wlog <- wlog(.log_file, with_pid = TRUE)

# One height's full body, as a named function so both the main pass and
# the retry pass below can share it. Returns invisible(NULL) on success,
# or the caught error condition on failure -- the caller inspects which.
.run_height <- function(i) {
  h          <- heights[i]
  h_key      <- sprintf("h%.2f", h)
  h_rds_path <- file.path(height_dir, sprintf("%s.rds", h_key))
  if (file.exists(h_rds_path)) {
    .wlog(sprintf("[%d/%d] %.2f m — skipped.", i, n_heights, h))
    return(invisible(NULL))
  }
  .wlog(sprintf("[%d/%d] %.2f m — starting runmicro (full year)...", i, n_heights, h))
  # 2026-08-10: give THIS forked worker its own tempdir subdirectory,
  # keyed by pid, before touching terra at all. Root-caused a real
  # production failure the same day: every mclapply worker inherits the
  # SAME TMPDIR from the parent at fork time (fork() doesn't regenerate
  # tempdir()'s per-session random subdirectory the way a fresh R start
  # would), so with a shared tempdir, one worker's on.exit() cleanup below
  # was deleting OTHER concurrently-running workers' still-in-use temp
  # rasters -- "[mask] file does not exist" failures that struck even at
  # mc.cores=2, not just 4 (Maquipucuna, run_microclimate_26943387,
  # 2026-08-10: 51/76 heights lost to this across the main pass + retry
  # before the run correctly stop()ed rather than save a truncated
  # manifest). Idempotent and cheap to call on every height, not just the
  # worker's first -- terraOptions() just re-sets the same path each time.
  worker_tmp <- file.path(Sys.getenv("TMPDIR", unset = tempdir()), paste0("worker_", Sys.getpid()))
  dir.create(worker_tmp, showWarnings = FALSE, recursive = TRUE)
  terra::terraOptions(tempdir = worker_tmp)
  # Clear THIS worker's own (now properly isolated) scratch rasters on the
  # way out regardless of success/failure, so a many-height run doesn't
  # accumulate temp files unboundedly.
  on.exit(terra::tmpFiles(remove = TRUE), add = TRUE)
  result <- tryCatch({
    dtm_  <- terra::unwrap(dtmdata_w)
    dtmc_ <- terra::unwrap(dtmc_w)
    # method="Cpp": diag_wrap_method_test.R (2026-08-08/09) confirmed this is
    # numerically identical to "R" (max|diff|=0 across Tz/relhum/windspeed/
    # Rdirdown/Rdifdown, ~1.2e9 elements each, identical NA/NaN positions) and
    # ~8% faster; diag_wrap_concurrency.R confirmed no /tmp wrap/unwrap
    # collision recurrence under this method with mc.cores 2 and 4.
    mout <- microclimf::runmicro(micropoint = model, reqhgt = h,
      vegp = vegetationdata, soilc = soildata, dtm = dtm_, dtmc = dtmc_,
      altcorrect = 1, method = "Cpp", out = unname(OUT_VARS))
    # runmicro() always ends with `mout$tme <- as.POSIXct(micropoint[[1]]$
    # tmeorig)` -- the ORIGINAL, unsubsetted timestamps, regardless of the
    # subsetpointmodela() call above. Left alone this claims length(tme)
    # timesteps for arrays that actually have length(subset_tme), silently
    # mis-shaping every matrix(mout$X, ncol=length(mout$tme)) reshape below
    # and .compute_voxel_quantiles()'s month masks. The stopifnot converts
    # any such mismatch into a caught error (.run_height()'s tryCatch logs
    # it), which the n_done < n_heights hard-fail below turns into a loud
    # job failure instead of ever writing a corrupted manifest.
    mout$tme <- subset_tme
    stopifnot(dim(mout$Tz)[3] == length(mout$tme))
    # Reduce in-memory here, before anything is written to disk: raw hourly
    # per-pixel arrays (Tz/relhum/windspeed/Rdirdown/Rdifdown, each
    # (nrow, ncol, ~8760)) filled the shared Lustre scratch workspace
    # (91% used cluster-wide, 2026-08-04) and caused repeated
    # "No space left on device" job failures after only a fraction of a
    # site's heights completed. Two products are kept instead of the raw
    # arrays, computed once here while the raw arrays are still in memory,
    # then the raw arrays are discarded (never saveRDS()'d):
    #   - voxel_quantiles: per-pixel/month/daypart quantiles for EVERY pixel
    #     (.compute_voxel_quantiles(), get_colonization.R) -- feeds the
    #     niche-scoring/per-voxel establishment system (get_clim_voxel()).
    #     Computed over the whole raster, not one run's landscape footprint,
    #     since no specific colonization run/landscape exists yet at this
    #     point -- get_clim_voxel() subsets down to a run's own footprint at
    #     model-run time, a cheap lookup rather than a recomputation.
    #   - *_mean: spatially-averaged (one value per hour, not per pixel)
    #     series for lookup_climate_by_height()/build_clim_cache() -- that
    #     pathway (wind attenuation, precip/survival fallback means) was
    #     already spatially flattened by design, so nothing is lost moving
    #     its averaging step here instead of computing it lazily at read
    #     time from the (now nonexistent) raw array. O(ntime), not
    #     O(pixels*ntime) -- negligible size next to voxel_quantiles.
    nr <- terra::nrow(dtm_)
    nc <- terra::ncol(dtm_)
    voxel_quantiles <- .compute_voxel_quantiles(mout, nr, nc)
    # v7 3D re-run (Phase A1): per-pixel MEANS (not quantiles -- no
    # quantile-to-mean estimator anywhere in this pipeline, per the
    # decision), restricted to each named site's own small footprint
    # (site_footprints, built once above) rather than the whole raster.
    # This is the input the 3D demographic core (Phase B) will read; the
    # existing voxel_quantiles above is untouched and keeps feeding the
    # niche-establishment system exactly as before.
    pixel_means <- .compute_pixel_means(mout, nr, nc, site_footprints)
    saveRDS(
      list(
        tme             = mout$tme,
        voxel_quantiles = voxel_quantiles,
        pixel_means     = pixel_means,
        # 2026-08-10: apply(arr, 3, mean, na.rm=TRUE) -> matrixStats::
        # colMeans2() on the same (npix, ntime) reshape .compute_voxel_
        # quantiles() already uses -- verified numerically identical
        # (incl. all-NaN-column and Inf edge cases) before this change,
        # see /tmp/test_matrixstats_equivalence.R from that verification
        # pass. apply()'s per-slice generic dispatch is slow at this
        # scale; colMeans2() is the same spatial-average-per-hour
        # reduction, just vectorized.
        temp_mean       = matrixStats::colMeans2(matrix(mout$Tz, nrow = nr * nc, ncol = length(mout$tme)), na.rm = TRUE),
        relhum_mean     = matrixStats::colMeans2(matrix(mout$relhum, nrow = nr * nc, ncol = length(mout$tme)), na.rm = TRUE),
        windspeed_mean  = matrixStats::colMeans2(matrix(mout$windspeed, nrow = nr * nc, ncol = length(mout$tme)), na.rm = TRUE),
        Rdirdown_mean   = matrixStats::colMeans2(matrix(mout$Rdirdown, nrow = nr * nc, ncol = length(mout$tme)), na.rm = TRUE),
        Rdifdown_mean   = matrixStats::colMeans2(matrix(mout$Rdifdown, nrow = nr * nc, ncol = length(mout$tme)), na.rm = TRUE)
      ),
      h_rds_path
    )
    .wlog(sprintf("[%d/%d] %.2f m — done. %d timesteps, %s .. %s. Saved: %.1f MB.",
                  i, n_heights, h, length(mout$tme), format(min(mout$tme)), format(max(mout$tme)),
                  file.info(h_rds_path)$size / 1e6))
    invisible(NULL)
  }, error = function(e) e)
  if (inherits(result, "error")) {
    .wlog(sprintf("[%d/%d] %.2f m — FAILED: %s", i, n_heights, h, conditionMessage(result)))
  }
  result
}

# Which of `heights` already has its .rds file on disk right now -- shared
# by the main pass, the retry pass, and the retry's own completion check
# below, which previously recomputed this file.exists() lookup 3x verbatim.
.height_rds_exists <- function(height_dir, heights) {
  file.exists(file.path(height_dir, sprintf("h%.2f.rds", heights)))
}

results <- mclapply(seq_along(heights), .run_height, mc.cores = n_cores)

heights_ok <- heights[.height_rds_exists(height_dir, heights)]
n_done <- length(heights_ok)
log_msg(sprintf("Height loop done: %d/%d complete.", n_done, n_heights))

# 2026-08-10: automatic retry pass for whatever's still missing, before
# deciding pass/fail below. Reduced concurrency (half the cores, minimum
# 1) so a resource-exhaustion failure (e.g. the /tmp exhaustion that
# silently dropped Maquipucuna's heights 17-76 on 2026-08-09) isn't just
# immediately reproduced at the same concurrency that caused it.
if (n_done < n_heights) {
  missing_idx <- which(!.height_rds_exists(height_dir, heights))
  retry_cores <- max(1L, n_cores %/% 2L)
  log_msg(sprintf("Retrying %d/%d failed heights on %d core(s)...",
                   length(missing_idx), n_heights, retry_cores))
  retry_results <- mclapply(missing_idx, .run_height, mc.cores = retry_cores)
  heights_ok <- heights[.height_rds_exists(height_dir, heights)]
  n_done <- length(heights_ok)
  log_msg(sprintf("Retry done: %d/%d complete overall.", n_done, n_heights))
}

# A run where every worker failed (e.g. the 2026-07-24 Saloya "weather[[k]] :
# subscript out of bounds" incident) must not save a manifest at all --
# previously .heights was set to the full INTENDED seq() regardless of what
# actually got written, so a 0/500 run still produced a manifest claiming
# 500 valid heights and exited 0, and downstream code silently treated the
# site as done. Failing loudly here lets the pipeline's own afterok
# dependency chain cancel everything downstream, as designed.
if (n_done == 0) {
  stop(sprintf("Height loop produced 0/%d usable height files in %s -- aborting without saving a manifest.",
               n_heights, height_dir))
}
# 2026-08-10: partial completion now fails just as loudly as total failure
# (previously only a WARNING, and the run still exited 0 with a truncated
# manifest -- exactly what let Maquipucuna's 16/76-height run pass as
# "COMPLETED" and, via microenv_array.sh's --dependency=afterok chaining,
# let Mashpi start immediately behind it carrying the same unaddressed
# risk). A manifest capped at the wrong height ceiling is not a safe
# degraded result: load_height() (get_colonization.R) does a nearest-
# height snap with no tolerance, so every downstream query above the
# truncated ceiling would silently return the top tier's climate instead
# of erroring.
if (n_done < n_heights) {
  stop(sprintf(
    "Height loop produced only %d/%d usable height files in %s after retry -- aborting without saving a manifest.",
    n_done, n_heights, height_dir))
}

# ── 4. Save manifest ─────────────────────────────────────────────────────────
# Each height file is ~6 GB in memory (10 spatial arrays × 374×372×288).
# Loading all heights at once would require hundreds of GB — instead save a
# small manifest that points to the per-height files in scratch. Downstream
# code reads individual heights on demand via load_height().
log_msg("Saving microenv manifest...")
manifest <- list(
  .heights    = heights_ok,
  .height_dir = height_dir,
  # Extent stored as a plain numeric vector (not a raw terra::ext() S4 object)
  # -- a SpatExtent's only slot is a C++ pointer that does not survive
  # saveRDS()/readRDS() across sessions (raises "NULL value passed as symbol
  # address" on reload) unless wrapped with terra::wrap(), which SpatExtent
  # objects don't support. terra::crs() already returns a plain character
  # (WKT) string, so it round-trips fine as-is.
  .spatial    = list(
    ext  = as.vector(terra::ext(dtmdata)),
    nrow = terra::nrow(dtmdata),
    ncol = terra::ncol(dtmdata),
    crs  = terra::crs(dtmdata)
  ),
  # Deliberately the FULL, pre-subsample record captured as weather_full
  # above -- NOT model[[1]]$weather, which was subsampled to every other
  # day for the grid expansion. Downstream precip/winddir
  # (get_colonization.R:240, lookup_climate_by_height()) needs the
  # complete hourly series.
  .weather    = weather_full,
  .era5_grid_snap = if (exists("snap")) snap$era5_grid_snap else NULL,
  .n_valid_era5_cells = N_VALID_ERA5_CELLS,
  .n_valid_era5_cells_monthly = N_VALID_ERA5_CELLS_MONTHLY,
  .n_total_era5_cells = N_TOTAL_ERA5_CELLS
)
saveRDS(manifest, site_env_path)
log_msg(sprintf("Manifest saved to %s  (height_dir=%s)", site_env_path, height_dir))

# v7 3D re-run (Phase A1): also save one manifest PER REAL SITE NAME (not
# just the combined TARGET_SITE, which for a cluster run is a synthetic
# "MindoCluster_..." label) -- same shared .heights/.height_dir/.spatial/
# .weather (the 4 cluster sites deliberately share identical forcing, see
# above), plus this site's own footprint pixel indices so
# build_clim_cache()/build_clim_cache_voxel() can find exactly its own
# pixels inside the shared height files. Existing code that loads
# `microenv_<SiteName>.rds` by name (get_colonization.R's load_height() and
# everything built on it) needs no changes for this to work -- it just
# starts finding a real file where before (for the 3 satellite cluster
# sites) there wasn't one from this run.
if (IS_CLUSTER) {
  for (real_site in CLUSTER_SITES) {
    real_manifest <- manifest
    real_manifest$.footprint_idx <- site_footprints[[real_site]]
    real_manifest$.site_name     <- real_site
    real_manifest$.cluster_of    <- setdiff(CLUSTER_SITES, real_site)
    real_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", real_site, res_suffix))
    saveRDS(real_manifest, real_path)
    log_msg(sprintf("Cluster member manifest saved: %s (%d footprint px, shared height_dir)",
                    real_path, length(site_footprints[[real_site]])))
  }
} else {
  manifest$.footprint_idx <- site_footprints[[TARGET_SITE]]
  manifest$.site_name     <- TARGET_SITE
  saveRDS(manifest, site_env_path)
}
log_msg("Done. Height files remain in scratch — do not ws_release until downstream analysis is complete.")
