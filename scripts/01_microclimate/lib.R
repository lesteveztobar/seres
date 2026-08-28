# lib.R — microclimate acquisition and post-processing library
# Merged from get_climateinputs.R (ERA5/DTM/LAI/albedo/vegetation/soil acquisition)
# and get_microenv.R (niche extraction, canopy grid, climate lookups).
# Sourced together by run_microclimate.R, run_microclimate_site.R,
# regenerate_missing_dtm.R, and scripts/legacy/allsites.R.
# ─────────────────────────────────────────────────────────────────────────────

# ── Worker logging ────────────────────────────────────────────────────────────

# Returns a closure that timestamps `msg` and appends it to `log_file` — used
# by mclapply worker loops to log to a shared file without every worker
# opening its own connection. Set with_pid = TRUE to prefix the worker's PID,
# which matters when several workers write to the same log concurrently
# (e.g. the diag_wrap_*.R scripts); the main production height loop in
# run_microclimate_site.R logs without it. Formerly copy-pasted as a local
# .wlog() in run_microclimate_site.R, diag_wrap_collision.R, and
# diag_wrap_method_test.R — consolidated here.
wlog <- function(log_file, with_pid = FALSE) {
  tag <- if (with_pid) sprintf("[worker pid=%d]", Sys.getpid()) else "[worker]"
  function(msg) {
    stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "]", tag, " ", msg)
    cat(stamped, "\n", file = log_file, append = TRUE)
  }
}

# ── Height ceiling ────────────────────────────────────────────────────────────

# The model's height ceiling must be the canopy top, not the tallest
# recorded epiphyte observation (hobs_max) -- the latter only reflects where
# individuals happened to be found by observers, not how tall the forest
# actually is, and silently truncates the microclimate/landscape grid below
# the real canopy (e.g. Maquipucuna's forest is ~14.7 m tall on average, but
# hObs_max there is only 5.5 m). Preference order:
#   1. `measured` -- measured CanopyHeight_m from combinedv3.csv
#      (make_sites(), this file), when available: real field measurements
#      at the observation points, more trustworthy than a remote-sensing
#      product for this specific forest.
#   2. 99th percentile of the GEE canopy-height raster already downloaded
#      for microclimf's vegetation parameters (vhgt.tif, at `vhgt_path`) --
#      robust to single-pixel outliers, unlike a bare max().
#   3. `hobs_max`, if neither of the above is available.
# Never goes below hobs_max regardless of source (every observed individual
# must remain inside the modelled height range).
#
# Formerly duplicated (once inline in run_microclimate_site.R, once as
# check_microenv_progress.R's own height_ceiling_for()) -- consolidated
# here. `log_fn` is optional: pass log_msg()/wlog() for the descriptive,
# per-branch messages the real pipeline run wants; leave it NULL for a
# quiet call (check_microenv_progress.R's lightweight report use case).
height_ceiling <- function(measured, vhgt_path, hobs_max, log_fn = NULL) {
  log_it <- if (is.null(log_fn)) function(...) invisible(NULL) else log_fn
  if (!is.null(measured) && is.finite(measured)) {
    log_it(sprintf("Canopy height (measured, CanopyHeight_m): %.1f m (hObs_max was %.1f m)",
                   measured, hobs_max))
    return(max(measured, hobs_max))
  }
  if (!requireNamespace("terra", quietly = TRUE) || !file.exists(vhgt_path)) {
    log_it("No measured CanopyHeight_m and no vhgt.tif -- falling back to hObs_max.")
    return(hobs_max)
  }
  vals <- tryCatch(terra::values(terra::rast(vhgt_path), na.rm = TRUE), error = function(e) numeric(0))
  if (length(vals) == 0) {
    log_it("No measured CanopyHeight_m and vhgt.tif has no valid values -- falling back to hObs_max.")
    return(hobs_max)
  }
  ceiling_from_canopy <- as.numeric(quantile(vals, 0.99, na.rm = TRUE))
  log_it(sprintf("No measured CanopyHeight_m -- using p99 of vhgt.tif: %.1f m (hObs_max was %.1f m)",
                 ceiling_from_canopy, hobs_max))
  max(ceiling_from_canopy, hobs_max)
}

# ═══════════════════════════════════════════════════════════════════════════
# Part 1: data acquisition (formerly get_climateinputs.R)
# ═══════════════════════════════════════════════════════════════════════════
# get_climateinputs.R
# Data acquisition pipeline for canopymicroenv
# Lizeth Estévez Tobar — University of Bonn, 2026

# ── Site preparation ──────────────────────────────────────────────────────────

# Reads the combined CSV and derives a per-site summary dataframe with one row
# per field site (Area_or_Site), each with:
# - bounding box padded by pad° to ensure ERA5 cells overlap the study area
# - time window from actual observation datetimes at that site
# - observed height range for the model height sequence
# Use this for the per-site loop in runmicroenv.R; use make_site() when running
# a single AllSites bounding box instead.
make_sites <- function(csv_path, pad = 0.01) {
  message("Reading combined CSV...")
  df <- read_csv(csv_path,
                 na        = c("", "NA", "N/A"),
                 col_types = COMBINED_COL_TYPES) |>
    dplyr::filter(!is.na(Source), !is.na(Area_or_Site))

  df |>
    dplyr::filter(!is.na(lat), !is.na(lon), !is.na(datetime), !is.na(Height_m)) |>
    dplyr::mutate(
      lat      = as.numeric(lat),
      lon      = as.numeric(lon),
      datetime = as.POSIXlt(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    ) |>
    dplyr::group_by(Area_or_Site) |>
    dplyr::summarise(
      lat_min     = min(lat)      - pad,
      lat_max     = max(lat)      + pad,
      lon_min     = min(lon)      - pad,
      lon_max     = max(lon)      + pad,
      tme_start   = min(datetime),
      tme_end     = max(datetime),
      hObs_min    = min(Height_m, na.rm = TRUE),
      hObs_max    = max(Height_m, na.rm = TRUE),
      # Measured canopy height at this site, for the model's height ceiling
      # (run_microclimate_site.R) — NA (via suppressWarnings on an all-NA
      # max()) if no CanopyHeight_m records exist here, in which case that
      # script falls back to the vhgt.tif remote-sensing raster instead.
      hCanopy_max = suppressWarnings(max(CanopyHeight_m, na.rm = TRUE)),
      .groups     = "drop"
    ) |>
    dplyr::mutate(hCanopy_max = ifelse(is.finite(hCanopy_max), hCanopy_max, NA_real_)) |>
    dplyr::rename(Site = Area_or_Site)
}

# ── ERA5 data acquisition ─────────────────────────────────────────────────────

# ERA5 processing utilities. Restored 2026-08-08 after an in-progress cleanup
# deleted these along with make_site()/get_clim()/get_clim_month()/
# get_canopy_grid() but left get_weather() (below) still calling three of
# them -- concat_era5_nc(), fix_lsm(), .download_era5_months() -- which broke
# get_weather() for any site needing a fresh/extended ERA5 download.
# .download_era5_months() also calls merge_era5_steptype_files(), a fourth
# dependency not directly referenced by get_weather() itself, so it is
# restored here too even though it wasn't named in the original deletion
# report. Order below is dependency order: concat_era5_nc and
# merge_era5_steptype_files are leaf utilities; fix_lsm is used by both
# .download_era5_months and get_weather; .download_era5_months is get_weather's
# direct dependency.

# Concatenates per-month ERA5 nc files into one multi-month file along the
# time dimension. All files must share the same spatial grid and variables
# (produced by merge_era5_steptype_files). Used when the requested time
# window spans more than one calendar year.
concat_era5_nc <- function(infiles, outfile) {
  library(ncdf4)
  message("Concatenating ", length(infiles), " ERA5 yearly nc files...")

  ncs <- lapply(infiles, nc_open)
  src <- ncs[[1]]

  # Identify time dimension by name ("valid_time" or "time"); CDS nc files do
  # not mark it unlimited, so we cannot rely on the unlim flag.
  tname <- names(src$dim)[grepl("time", names(src$dim), ignore.case = TRUE)][1]
  if (is.na(tname)) stop("concat_era5_nc: no time dimension found in ", infiles[1])

  # Concatenate time values across all files
  all_tvals <- unlist(lapply(ncs, function(nc) nc$dim[[tname]]$vals))

  # Build output dimensions
  out_dims <- lapply(src$dim, function(d) {
    if (d$name == tname)
      ncdim_def(d$name, d$units, all_tvals, unlim = TRUE)
    else
      ncdim_def(d$name, d$units, d$vals, unlim = FALSE)
  })
  names(out_dims) <- names(src$dim)

  # Build output variables (same structure as source)
  out_vars <- lapply(names(src$var), function(vname) {
    var_meta <- src$var[[vname]]
    vdims <- lapply(var_meta$dim, function(d) out_dims[[d$name]])
    ncvar_def(vname, var_meta$units, vdims, var_meta$missval)
  })
  names(out_vars) <- names(src$var)

  nc_out <- nc_create(outfile, vars = out_vars)
  for (vname in names(src$var)) {
    var_meta <- src$var[[vname]]
    is_tvar <- any(sapply(var_meta$dim, function(d) d$name == tname))
    if (is_tvar) {
      pieces   <- lapply(ncs, function(nc) ncvar_get(nc, vname))
      combined <- abind::abind(pieces, along = length(dim(pieces[[1]])))
      ncvar_put(nc_out, vname, combined)
    } else {
      ncvar_put(nc_out, vname, ncvar_get(src, vname))
    }
  }
  nc_close(nc_out)
  lapply(ncs, nc_close)
  message("Concatenated -> ", outfile)
}

# Merges the three ERA5 stepType netCDF files (accum, avg, instant) that CDS
# delivers separately into one combined file. Renames radiation variables to
# match the names expected by microclimdata::era5_process().
merge_era5_steptype_files <- function(pathin, pathout) {
  library(ncdf4)

  nc_files <- list.files(pathin, pattern = "stepType.*\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) stop("No stepType .nc files found in ", pathin)
  message("Merging ", length(nc_files), " ERA5 stepType files...")

  # CDS renamed two radiation variables between API versions
  rename_map <- c(avg_snlwrf = "msnlwrf", avg_sdlwrf = "msdwlwrf")

  datasets  <- lapply(nc_files, nc_open)
  src       <- datasets[[1]]
  out_dims  <- lapply(src$dim, function(d) ncdim_def(d$name, d$units, d$vals, unlim = d$unlim))
  names(out_dims) <- names(src$dim)

  seen_vars <- names(src$dim)
  out_vars  <- list()
  var_data  <- list()

  for (ds in datasets) {
    for (vname in names(ds$var)) {
      if (vname %in% seen_vars) next
      seen_vars <- c(seen_vars, vname)
      var       <- ds$var[[vname]]
      out_name  <- ifelse(vname %in% names(rename_map), rename_map[vname], vname)
      var_dims  <- lapply(var$dim, function(d) out_dims[[d$name]])
      out_vars[[out_name]] <- ncvar_def(name = out_name, units = var$units,
                                        dim = var_dims, missval = var$missval)
      var_data[[out_name]] <- ncvar_get(ds, vname)
      message("  ", vname, " -> ", out_name)
    }
  }

  nc_out <- nc_create(pathout, vars = out_vars)
  for (vname in names(var_data)) ncvar_put(nc_out, vname, var_data[[vname]])
  nc_close(nc_out)
  lapply(datasets, nc_close)
  message("Merged file written to ", pathout)
  return(pathout)
}

# ERA5 land-sea mask (lsm) sometimes has near-land cells with values just below 1
# (e.g. 0.96) which microclimdata treats as ocean and excludes.
# This fix rounds near-land cells up to 1 so they are included in processing.
fix_lsm <- function(nc_path) {
  nc  <- ncdf4::nc_open(nc_path, write = TRUE)
  lsm <- ncdf4::ncvar_get(nc, "lsm")
  n   <- sum(lsm < 1 & lsm >= 0.95)
  lsm[lsm >= 0.95] <- 1
  ncdf4::ncvar_put(nc, "lsm", lsm)
  ncdf4::nc_close(nc)
  message("LSM fix: ", n, " near-land cells set to 1")
}

# Downloads (or reuses already-cached) per-month ERA5 nc files for the given
# subset of `req` entries (as returned by mcera5::build_era5_request). Returns
# the vector of month_nc file paths, in the same order as `req`.
.download_era5_months <- function(req, site, credentials, dir, overwrite) {
  month_nc_files <- character(length(req))

  for (i in seq_along(req)) {
    req_year  <- req[[i]]$year
    req_month <- req[[i]]$month
    month_nc  <- file.path(dir, sprintf("%s_%s_%s.nc", site$Site, req_year, req_month))

    if (!file.exists(month_nc) || overwrite) {
      month_tmp <- file.path(dir, sprintf("era5_tmp_%s_%s", req_year, req_month))
      dir.create(month_tmp, showWarnings = FALSE)

      message("Downloading ERA5 for ", site$Site,
              " (", req_year, "-", req_month, ")...")
      ecmwfr::wf_request(
        request  = req[[i]],
        user     = credentials$username[credentials$Site == "CDS"],
        transfer = TRUE, path = paste0(month_tmp, "/"), retry = 120, verbose = TRUE
      )
      for (zip_file in list.files(month_tmp, pattern = "\\.zip$", full.names = TRUE)) {
        unzip(zip_file, exdir = month_tmp)
        unlink(zip_file)
      }

      merge_era5_steptype_files(pathin = month_tmp, pathout = month_nc)
      fix_lsm(month_nc)
      unlink(month_tmp, recursive = TRUE)
    } else {
      message("ERA5 for ", site$Site, " ", req_year, "-", req_month,
              " already cached, skipping.")
      fix_lsm(month_nc)
    }
    month_nc_files[i] <- month_nc
  }
  month_nc_files
}

# Downloads ERA5 hourly climate data for the site bounding box and time window.
# Submits one CDS request per month (by_month = TRUE), each to its own temp dir
# to avoid stepType filename collisions, merges them individually, then
# concatenates all months into the final site nc. Resumes cleanly if interrupted.
#
# Incremental append: a sidecar "<site>_months.txt" manifest (one "<year>_<month>"
# token per line, matching the <site>_<year>_<month>.nc per-month naming above)
# records which months are already folded into merged_file. If the requested
# tme window includes months not yet in the manifest, only those missing
# months are downloaded, and the existing merged file is concatenated with
# them into a new merged file — no re-download of months already covered.
# Existing merged files that predate this manifest (no "<site>_months.txt"
# next to them) are left as-is with a warning, since we can't safely tell
# what months they already cover without it — delete or rename them to force
# a clean re-download instead.
#
# months_file lives ONE LEVEL UP from `dir` (the era5/ directory itself), not
# inside it -- `dir` is passed straight through as microclimdata::era5_process()'s
# pathin later in run_microclimate_site.R, which scans every file in that
# directory expecting only per-site .nc files (see the era5_tmp_* cleanup a
# few lines below, which exists for the same reason). A months.txt sitting
# alongside Saloya.nc there broke era5_process() with "NetCDF: Unknown file
# format" trying to nc_open() the manifest itself (2026-07-22, Saloya's
# first-time generation -- every earlier site predates this incremental-
# append feature entirely, so none of them had ever hit this).
get_weather <- function(site, credentials, r, tme, dir, overwrite = FALSE, output = "point") {
  message("")
  merged_file <- file.path(dir, paste0(site$Site, ".nc"))
  months_file <- file.path(dirname(dir), paste0(site$Site, "_months.txt"))

  # Directory-based mutex: dir.create() is atomic on POSIX filesystems
  # (including Lustre), so this safely serializes every call to this
  # function for the SAME site, regardless of which script/job triggers it
  # -- not just the two orchestration scripts that were fixed to chain
  # height-step jobs sequentially. This closes the actual root cause behind
  # the 2026-07-24 Saloya corruption: fix_lsm() opens merged_file with
  # write=TRUE on EVERY call (even the "already covers requested window,
  # nothing to download" branch), and three concurrent Saloya jobs each
  # calling that -- plus one of them concurrently extending/renaming the
  # same merged_file -- corrupted it badly enough that every subsequent
  # read failed downstream ("weather[[k]] : subscript out of bounds",
  # every height tier, even in a later single non-concurrent run -- the
  # damage was already baked into the file, sequencing alone didn't fix it
  # after the fact). Lock released on exit regardless of success/error.
  lock_dir <- file.path(dirname(dir), paste0(".", site$Site, "_era5.lock"))
  waited <- 0
  while (!dir.create(lock_dir, showWarnings = FALSE)) {
    Sys.sleep(5)
    waited <- waited + 5
    if (waited %% 60 == 0) {
      message("Waiting for ERA5 lock on ", site$Site, " (", waited, "s so far, held by another job)...")
    }
    if (waited > 3600) {
      stop("Timed out after 1h waiting for the ERA5 lock on ", site$Site,
           " -- if no other job for this site is actually running, a previous ",
           "run may have crashed without releasing it; remove manually: ", lock_dir)
    }
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)

  req <- mcera5::build_era5_request(
    xmin       = site$lon_min, xmax = site$lon_max,
    ymin       = site$lat_min, ymax = site$lat_max,
    start_time = site$tme_start, end_time = site$tme_end,
    by_month   = TRUE, outfile_name = site$Site
  )
  req_keys <- vapply(req, function(req_item) sprintf("%s_%s", req_item$year, req_item$month), character(1))

  if (!file.exists(merged_file) || overwrite) {
    new_month_files <- .download_era5_months(req, site, credentials, dir, overwrite = TRUE)

    if (length(new_month_files) == 1L) {
      file.rename(new_month_files, merged_file)
    } else {
      concat_era5_nc(new_month_files, merged_file)
      unlink(new_month_files)
    }
    writeLines(req_keys, months_file)
    message("Merged ERA5 file ready: ", merged_file)

  } else if (!file.exists(months_file)) {
    message("ERA5 merged file already exists for ", site$Site,
            " but has no months manifest (predates incremental-append support) — ",
            "leaving it as-is. Delete/rename ", merged_file,
            " to force a clean re-download if it doesn't cover the full requested window.")
    fix_lsm(merged_file)

  } else {
    covered      <- readLines(months_file, warn = FALSE)
    missing_keys <- setdiff(req_keys, covered)

    if (length(missing_keys) == 0) {
      message("ERA5 merged file already covers the full requested window for ",
              site$Site, ", skipping download.")
      fix_lsm(merged_file)
    } else {
      message("Extending existing ERA5 merged file for ", site$Site, " with ",
              length(missing_keys), " new month(s): ", paste(missing_keys, collapse = ", "))
      new_month_files <- .download_era5_months(req[match(missing_keys, req_keys)],
                                                 site, credentials, dir, overwrite = overwrite)

      # PID-suffixed so two concurrent callers for the same site (e.g. two
      # height-tier-spacing jobs both hitting this branch at once) never
      # write to the same temp path -- see 2026-07-24: three Saloya
      # run_microenv.sh jobs (h0.1/h0.5/h1.0) all extended this file within
      # the same second and corrupted it (weather[[k]] subscript-out-of-
      # bounds downstream, every height tier failed). This alone doesn't
      # make concurrent extends fully safe -- the callers themselves must
      # not run concurrently for the same site (see run_full_analysis_
      # pipeline.sh / height_res_array.sh, now chained sequentially per
      # site) -- it just stops two processes clobbering the exact same file
      # mid-write, which was the sharper edge of the race.
      tmp_out <- file.path(dir, paste0(site$Site, "_extended_tmp_", Sys.getpid(), ".nc"))
      concat_era5_nc(c(merged_file, new_month_files), tmp_out)
      unlink(new_month_files)
      file.rename(tmp_out, merged_file)

      writeLines(union(covered, missing_keys), months_file)
      message("Merged ERA5 file extended: ", merged_file)
    }
  }

  # Remove any leftover era5_tmp_* directories that would confuse era5_process()
  stale_tmp <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
  stale_tmp <- stale_tmp[grepl("era5_tmp_", basename(stale_tmp))]
  if (length(stale_tmp) > 0) {
    message("Cleaning up stale temp dirs: ", paste(basename(stale_tmp), collapse = ", "))
    unlink(stale_tmp, recursive = TRUE)
  }

  message("Processing ERA5 data to point climate data frame...")
  weatherdata <- microclimdata::era5_process(
    tme = tme, req = NA, pathin = paste0(dir, "/"), r = r, out = output
  )
  message("Weather data ready.")
  return(weatherdata)
}

# ── Terrain and landcover ─────────────────────────────────────────────────────

# Downloads a digital elevation model for the site extent via elevatr.
# Projects r to UTM first so the downloaded DEM has a sensible metric resolution,
# then reprojects the result back to WGS84 (EPSG:4326) to match all other inputs.
# Caches as dtm.tif to avoid re-downloading.
get_dtm <- function(r, dir, mask = FALSE) {
  message("")
  cache_file <- file.path(dir, "dtm.tif")

  if (file.exists(cache_file)) {
    message("DTM cache found, loading...")
    return(terra::rast(cache_file))
  }

  message("Downloading digital elevation model...")
  # dem_download resamples to the template raster, so a 2×2 template gives a
  # 2×2 DTM. Instead, build a ~90m resolution UTM template to get a usable DEM,
  # then reproject to WGS84. EPSG:32717 = UTM zone 17S covers western Ecuador.
  e_wgs <- terra::ext(r)
  r_utm <- terra::project(terra::rast(e_wgs, crs = "EPSG:4326"), "EPSG:32717")
  e_utm <- terra::ext(r_utm)
  r_utm90 <- terra::rast(
    xmin = e_utm$xmin, xmax = e_utm$xmax,
    ymin = e_utm$ymin, ymax = e_utm$ymax,
    res  = 90, crs = "EPSG:32717"
  )
  terra::values(r_utm90) <- 1
  dtm <- microclimdata::dem_download(r = r_utm90, msk = FALSE)
  dtm <- terra::project(dtm, "EPSG:4326")

  terra::writeRaster(dtm, cache_file)
  message("DTM downloaded and cached: ", nrow(dtm), " x ", ncol(dtm),
          " pixels at ", round(terra::res(dtm)[1] * 111320, 0), "m resolution")
  return(dtm)
}

# Downloads ESA WorldCover 10m landcover via Google Earth Engine, exports to
# Google Drive, and downloads to disk. Checks Drive first to avoid re-exporting.
# Requires rgee initialisation before calling.
get_landcover <- function(site, r, out_dir, type = "ESA",
                          overwrite = FALSE, google_drive_folder = "rgee_backup") {
  message("")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_path    <- file.path(out_dir, paste0(site$Site, "_landcover_", type, ".tif"))
  drive_prefix <- paste0(site$Site, "_ESA_WorldCover")

  if (file.exists(save_path) && !overwrite) {
    message("Landcover already exists on disk, loading...")
    return(terra::rast(save_path))
  }

  message("Checking Google Drive for landcover...")
  googledrive::drive_auth(email = "lizethestevezt@gmail.com", cache = "~/.secrets")
  folder      <- googledrive::drive_find(pattern = google_drive_folder, type = "folder", n_max = 1)
  drive_files <- googledrive::drive_ls(folder)
  drive_file  <- drive_files[grepl(drive_prefix, drive_files$name), ]

  if (nrow(drive_file) == 0 || overwrite) {
    message("Exporting landcover from GEE to Drive (this takes ~15 min)...")
    e          <- terra::ext(r)
    epsg_code  <- paste0("EPSG:", terra::crs(r, describe = TRUE)$code)
    aoi        <- ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
    aoi_coords <- aoi$bounds()$getInfo()$coordinates[[1]]
    img        <- ee$ImageCollection("ESA/WorldCover/v100")$first()
    task <- ee$batch$Export$image$toDrive(
      image = img, description = paste0(site$Site, "_landcover_export"),
      folder = google_drive_folder, fileNamePrefix = drive_prefix,
      region = aoi_coords, scale = 10, crs = epsg_code
    )
    task$start()
    rgee::ee_monitoring(task, max_attempts = 200, quiet = FALSE)
    drive_file <- googledrive::drive_ls(folder) |>
      dplyr::filter(grepl(drive_prefix, name))
  } else {
    message("Landcover found on Google Drive, downloading...")
  }

  googledrive::drive_download(file = drive_file[1, ], path = save_path, overwrite = TRUE)
  message("Landcover saved to ", save_path)
  return(terra::rast(save_path))
}

# ── Vegetation and soil parameters ────────────────────────────────────────────

# Downloads MODIS LAI (500m) for the site extent and time period via NASA
# Earthdata, then mosaics tiles using microclimdata::lai_mosaic().
# Skips download if HDF files already exist in pathout.
get_lai <- function(r, tme, pathout, credentials, reso = 500) {
  message("")
  dir.create(pathout, recursive = TRUE, showWarnings = FALSE)
  if (!reso %in% c(10, 500)) stop("reso must be one of 10 or 500")
  if (tme[length(tme)] < as.POSIXlt("2000-02-18", tz = "UTC"))
    stop("No MODIS data available prior to 2000-02-18")

  lai_file     <- file.path(pathout, "lai_mosaic.tif")
  existing_hdf <- list.files(pathout, pattern = "\\.hdf$", full.names = TRUE)

  if (!file.exists(lai_file)) {
    if (length(existing_hdf) == 0) {
      message("Downloading MODIS LAI...")
      microclimdata::lai_download(
        r           = r,
        tme         = tme,
        reso        = reso,
        pathout     = paste0(pathout, "/"),
        credentials = data.frame(
          username = credentials$username[credentials$Site == "NASA"],
          password = credentials$password[credentials$Site == "NASA"]
        )
      )
    } else {
      message("MODIS LAI files already exist (", length(existing_hdf), " files), skipping download.")
    }
    message("Mosaicing LAI tiles...")
    laidata <- microclimdata::lai_mosaic(r = r, pathin = paste0(pathout, "/"), reso = reso)
    terra::writeRaster(laidata, lai_file)
  } else {
    message("LAI mosaic cache found, loading...")
    laidata <- terra::rast(lai_file)
  }

  message("LAI ready: ", nrow(laidata), "x", ncol(laidata), " @ ", nlyr(laidata), " layers")
  return(laidata)
}

# Downloads and processes MODIS BRDF/albedo (WSA shortwave, band 30) for the
# site. Caches the processed SpatRaster as albedo_processed.rds to avoid
# re-downloading on subsequent runs.
get_albedo <- function(r, tme, pathout, credentials) {
  message("")
  dir.create(pathout, recursive = TRUE, showWarnings = FALSE)
  alb_cache <- file.path(pathout, "albedo_processed.rds")

  if (file.exists(alb_cache)) {
    message("Albedo cache found, loading...")
    return(readRDS(alb_cache))
  }

  message("Downloading MODIS albedo...")
  microclimdata::albedo_download(
    r           = r,
    tme         = tme,
    pathout     = paste0(pathout, "/"),
    credentials = credentials
  )
  message("Processing albedo...")
  albedodata <- microclimdata::albedo_process(r = r, pathin = paste0(pathout, "/"))

  saveRDS(albedodata, alb_cache)
  message("Albedo ready and cached.")
  return(albedodata)
}

# Computes ground and leaf reflectance from LAI, albedo, and the leaf inclination
# coefficient derived from landcover. Aggregates all inputs to a common coarse
# grid before calling reflectance_calc() to avoid geometry mismatch errors.
# Returns a list with $gref (ground reflectance) and $lref (leaf reflectance).
get_reflectance <- function(lai, alb, landcover, cachefile = NULL) {
  message("")

  if (!is.null(cachefile) && file.exists(cachefile)) {
    message("Reflectance cache found, loading...")
    cached <- readRDS(cachefile)
    return(list(gref = terra::unwrap(cached$gref), lref = terra::unwrap(cached$lref)))
  }

  # x_calc maps ESA landcover codes to per-pixel leaf inclination coefficients
  message("Computing leaf inclination coefficients...")
  x_lc <- microclimdata::x_calc(landcover = landcover, lctype = "ESA")

  # Collapse multi-layer LAI to a single mean layer, then align all three inputs
  # to x_lc's grid before passing to reflectance_calc
  message("Aggregating inputs to common grid for reflectance_calc...")
  lai_agg <- terra::resample(terra::app(lai, mean, na.rm = TRUE), x_lc)
  alb_agg <- terra::resample(alb, x_lc)

  message("Computing reflectance...")
  refldata <- microclimdata::reflectance_calc(
    lai          = lai_agg,
    alb          = alb_agg,
    x            = x_lc,
    plotprogress = FALSE
  )
  message("refldata$gref range: ", round(min(terra::values(refldata$gref), na.rm = TRUE), 3),
          " – ", round(max(terra::values(refldata$gref), na.rm = TRUE), 3))
  message("refldata$lref range: ", round(min(terra::values(refldata$lref), na.rm = TRUE), 3),
          " – ", round(max(terra::values(refldata$lref), na.rm = TRUE), 3))

  if (!is.null(cachefile)) {
    saveRDS(list(gref = terra::wrap(refldata$gref), lref = terra::wrap(refldata$lref)), cachefile)
    message("Reflectance cached.")
  }
  return(refldata)
}

# Builds the vegparams object for microclimf using microclimdata::create_veggrid().
# Downloads vegetation height from GEE if not already cached in dir.
# Resamples vhgt and lai to the landcover grid before creating vegp, as
# create_veggrid() requires all inputs to share the same geometry.
get_vegetation <- function(r, lcover, lai, refldata, dir, site_name) {
  message("")
  vhgt_file    <- file.path(dir, "vhgt.tif")
  drive_prefix <- paste0("canopy_height_", site_name)

  if (!file.exists(vhgt_file)) {
    googledrive::drive_auth(email = "lizethestevezt@gmail.com", cache = "~/.secrets")
    folder      <- googledrive::drive_find(pattern = "rgee_backup", type = "folder", n_max = 1)
    drive_files <- googledrive::drive_ls(folder)
    drive_file  <- drive_files[grepl(drive_prefix, drive_files$name), ]

    if (nrow(drive_file) == 0) {
      message("Vegetation height not found on Drive — exporting from GEE for ", site_name, "...")
      # patch vegheight_download to use a site-specific Drive filename
      .orig_vhgt <- get("vegheight_download", envir = getNamespace("microclimdata"))
      assignInNamespace("vegheight_download",
        function(r, GoogleDrivefolder, pathtopython, projectname = NA, silent = FALSE) {
          reticulate::use_python(pathtopython, required = TRUE)
          if (!is.na(projectname)) rgee::ee$Initialize(project = projectname)
          e  <- terra::ext(r)
          r2 <- terra::rast(e); terra::crs(r2) <- terra::crs(r)
          r2 <- terra::project(r2, "EPSG:4326"); e <- terra::ext(r2)
          proj_string <- terra::crs(r, describe = TRUE)
          epsg_code   <- paste0("EPSG:", proj_string$code)
          aoi         <- rgee::ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
          aoi_coords  <- aoi$bounds()$getInfo()$coordinates[[1]]
          canopy_height <- rgee::ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")
          task <- rgee::ee$batch$Export$image$toDrive(
            image          = canopy_height,
            description    = paste0("canopy_height_", site_name),
            folder         = GoogleDrivefolder,
            fileNamePrefix = drive_prefix,
            region         = aoi_coords,
            scale          = 10,
            crs            = epsg_code
          )
          task$start()
          if (!silent) microclimdata:::.monitor_task(task$id)
        },
        ns = "microclimdata"
      )
      microclimdata::vegheight_download(
        r                 = r,
        GoogleDrivefolder = "rgee_backup",
        pathtopython      = Sys.getenv("CANOPY_PYTHON",
                              unset = "/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python3.12"),
        projectname       = "ee-lizethestevezt"
      )
      # re-check Drive after export
      drive_files <- googledrive::drive_ls(folder)
      drive_file  <- drive_files[grepl(drive_prefix, drive_files$name), ]
      if (nrow(drive_file) == 0)
        stop("GEE export completed but ", drive_prefix, " not found on Drive")
    } else {
      message("Vegetation height found on Drive for ", site_name, ", downloading...")
    }

    googledrive::drive_download(file = drive_file[1, ], path = vhgt_file, overwrite = TRUE)
    message("Vegetation height downloaded and cached.")
  }
  vhgt <- terra::rast(vhgt_file)

  # All inputs to create_veggrid must share the same geometry — resample to lcover grid
  message("Resampling vhgt, lai, and lref to landcover grid...")
  vhgt_rs <- terra::resample(vhgt,         lcover,        method = "bilinear")
  lai_rs   <- terra::resample(lai,          lcover,        method = "bilinear")
  lref_rs  <- terra::resample(refldata$lref, lcover,       method = "bilinear")
  gref_rs  <- terra::resample(refldata$gref, lcover,       method = "bilinear")

  message("Creating vegparams...")
  vegetationdata <- microclimdata::create_veggrid(
    landcover = lcover,
    vhgt      = vhgt_rs,
    lai       = lai_rs,
    refldata  = list(gref = gref_rs, lref = lref_rs),
    lctype    = "ESA"
  )
  message("vegparams ready | pai layers: ", nlyr(terra::rast(vegetationdata$pai)))
  return(vegetationdata)
}

# Builds the soilcharac object for microclimf using microclimdata::create_soilgrid().
# Downloads SoilGrids physical properties if not already cached.
# refldata provides ground reflectance (gref); the soil type is derived internally
# by create_soilgrid() from the physical properties.
get_soil <- function(r, dir, landcover, refldata) {
  message("")
  soil_cache <- file.path(dir, "soilproperties.rds")

  if (!file.exists(soil_cache)) {
    message("Downloading SoilGrids data...")
    soil_r <- r
    terra::values(soil_r) <- 1
    soilprops <- microclimdata::soildata_download(
      r           = soil_r,
      pathdir     = paste0(dir, "/"),
      deletefiles = FALSE
    )
    saveRDS(soilprops, soil_cache)
  } else {
    message("Soil cache found, loading...")
    soilprops <- readRDS(soil_cache)
  }

  # gref and lref must match the landcover grid geometry for create_soilgrid
  message("Creating soilcharac...")
  gref_rs <- terra::resample(refldata$gref, landcover, method = "bilinear")
  lref_rs <- terra::resample(refldata$lref, landcover, method = "bilinear")
  soildata <- microclimdata::create_soilgrid(
    soildata  = soilprops,
    refldata  = list(gref = gref_rs, lref = lref_rs),
    landcover = landcover
  )
  # soildata_downscale assigns 0 to water/masked pixels; checkinputs requires 1–11
  st <- terra::rast(soildata$soiltype)
  st[st < 1 | st > 11] <- NA
  soildata$soiltype <- terra::wrap(st)

  message("soilcharac ready")
  return(soildata)
}