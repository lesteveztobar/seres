# sanity_gate.R
# Phase 2.2 sanity gate, per explicit instruction: HARD assertions that
# stop() the run, not warnings printed and ignored. Run this before every
# experiment (source it, or `Rscript scripts/diagnostics/sanity_gate.R`
# as a standalone check) -- a failure here means don't submit the
# downstream run yet.
#
# (a) every site with >=1 determined (identification-confirmed) record has
#     >=1 presence row in the niche cache
# (b) no observation has a blank/NA/default FinalID
# (c) every observation lies within its own site's footprint -- errors if
#     any observation is farther than DIST_THRESHOLD_M from its site's own
#     median lat/lon (catches exactly the Saloya-outlier class of bug)
# (d) every site's footprint area (from its own observation coordinate
#     range) is within a plausible range [AREA_MIN_HA, AREA_MAX_HA]
# (e) site_elev is finite for every site
# (f) niche-cache keys match the taxon-unit roster (post-filter, post-
#     rename FinalID values) by EXACT name -- not a superset/subset check
# (g) hour-1 valid-ERA5-cell count -- NOTE ONLY (2026-09-10). It is
#     microclimf::runpointmodela()'s own hour-1-only test and structurally
#     undercounts for ERA5 accumulated fields (MindoMirador 1/6 at hour 1
#     vs 6/6 across the year). The hard gate is (g2). A hard gate on a
#     known-unreliable metric fires falsely and trains overrides.
# (g2) HARD GATE: every site's manifest has >= MIN_VALID_ERA5_CELLS ERA5
#     grid cells valid across a 12-timestep yearly sample -- the
#     trustworthy count (run_microclimate_site.R computes it directly from
#     the raw ERA5 array, one timestep per calendar month, all finite).
# (h) at height tier 1 (the lowest, closest to ground), the per-pixel
#     voxel_quantiles array is at least MIN_FINITE_FRAC finite for temp --
#     catches the actual downstream symptom of a coverage collapse
#     directly, on every manifest regardless of vintage.
# (i) every height file the manifest names exists on disk, and the manifest
#     file is no older than its newest height file -- catches a Phase A job
#     that finished the per-height loop but silently exited without running
#     the manifest-save step (Saloya/LaElenita, 2026-09-10: COMPLETED 0:0,
#     no error, stale on-disk manifest pointing at an earlier run's data).
# (k) the manifest's terrain extent (.spatial$ext) encloses the model
#     domain (raw obs bbox + 0.15deg buffer), with the margin reported --
#     catches a stale filename-keyed terrain cache serving a DTM for a
#     different domain (the MindoMirador collapse, 2026-09-10).
# (j) the LANDSCAPE's physical extent (from ALL raw observations at the
#     site, as init_colonization() now derives it) and the tree count it
#     implies are within plausible bounds -- catches a site whose landscape
#     has collapsed to a handful of points (MindoMirador 21 raw obs / ~7 ha
#     -> 3 species-filtered / 0.40 ha / 120 trees, which check (d)'s 0.1 ha
#     floor let pass silently).
#
# Usage: source() this file (defines run_sanity_gate(), does not run it),
# or `Rscript scripts/diagnostics/sanity_gate.R [niche_cache_path]` to
# run it standalone against the given (or default) niche cache.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

DIST_THRESHOLD_M <- 2000   # generous vs. a real field plot (tens-hundreds of m); Saloya's outliers were ~7-8km out
AREA_MIN_HA <- 0.1
AREA_MAX_HA <- 200         # generous upper bound; a corrupted-coordinate blowup (Saloya pre-fix: 2805ha) fails this by >10x
# 2026-09-08: MindoMirador's regression was 1 valid cell (of 9) -- every
# healthy site/run seen in this project's own history has had >=2. A count-
# based floor (not a fraction) because the ERA5 grid cell denominator
# itself varies by site (4-9), and 1 valid cell isn't just "low fraction"
# -- it means zero genuine spatial contrast is even possible (every masked
# cell copies the one real point verbatim, runpointmodela()'s own fallback).
#
# LIMITATION (documented, not patched -- see the (g2) check below for the
# mitigation actually taken): this threshold is checked against TWO counts.
# The (g) count (.n_valid_era5_cells) comes from microclimf::
# runpointmodela()'s own per-cell test, which inspects ONLY the FIRST
# hourly ERA5 timestep at that cell (`is.na(climdf$temp[1]) == FALSE`) --
# not the whole year. That means a cell missing data only at hour 1 is
# wrongly discarded, and (more dangerous) a cell fine at hour 1 but gappy
# the rest of the year is wrongly kept. This project does not patch
# runpointmodela()/runmicro() (both compiled into the microclimf package;
# explicit decision, 2026-09-08) -- instead, run_microclimate_site.R
# computes a SECOND, independent count (.n_valid_era5_cells_monthly)
# directly from the raw ERA5 array, sampling one timestep per calendar
# month (12 total) and requiring ALL of them finite. Check (g2) asserts
# the SAME floor against that count. The two counts can differ (logged
# when they do) -- (g2) is the more trustworthy of the two; (g) is kept
# for continuity with what actually gates the point model itself.
MIN_VALID_ERA5_CELLS <- 2
# 2026-09-08: the healthy baseline observed directly (old MindoMirador: 51%
# finite for temp) is far above this; the corrupted run was 0.03%. Set well
# below the healthy baseline so this doesn't false-positive on a genuinely
# small/sparse site, while still catching anything resembling the collapse
# that just happened.
MIN_FINITE_FRAC <- 0.10
# (j) 2026-09-10: the LANDSCAPE's physical extent (init_colonization() ->
# site_landscape_bbox(), from ALL raw observations at the site) and the
# tree count it implies. Check (d) above computes the SPECIES-FILTERED
# footprint area and only floors it at 0.1 ha, so MindoMirador's
# 21-raw-obs / ~7 ha site collapsing to 3 filtered obs / 0.40 ha / 120
# trees passed silently. These bounds are on the real landscape: every
# modelled site here is a forest survey spanning several hectares.
LAND_AREA_MIN_HA <- 1.0
LAND_MIN_TREES   <- 100L
LAND_MIN_RAW_OBS <- 5L

run_sanity_gate <- function(niche_cache_path = NICHE_CACHE_PATH, obs_csv = OBSERVATIONS_CSV, verbose = TRUE) {
  fail <- function(msg) stop("SANITY GATE FAILED: ", msg, call. = FALSE)
  ok <- function(msg) if (verbose) message("  OK: ", msg)

  niches <- load_observations()
  niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                   !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
  sites <- sort(unique(niches$Area_or_Site))

  # (b) no blank/NA/default FinalID
  if (any(is.na(niches$FinalID) | trimws(niches$FinalID) == "")) {
    fail("(b) at least one observation has a blank/NA FinalID after load_observations() -- the unidentified-row filter did not run or was bypassed.")
  }
  if (any(niches$FinalID == "Maxillaria acutifolia" & is.na(niches$Identification))) {
    fail("(b) at least one observation has FinalID='Maxillaria acutifolia' with no raw Identification -- the old blank-defaults-to-acutifolia behavior has returned.")
  }
  ok("(b) no blank/NA/default FinalID")

  # (a)/(f) niche cache coverage and exact-name match
  if (!file.exists(niche_cache_path)) fail(sprintf("(a)/(f) niche cache not found at %s", niche_cache_path))
  niche_cache <- readRDS(niche_cache_path)
  cache_keys <- sort(names(niche_cache))
  roster <- sort(unique(niches$FinalID))
  missing_from_cache <- setdiff(roster, cache_keys)
  extra_in_cache <- setdiff(cache_keys, roster)
  if (length(missing_from_cache) > 0) {
    fail(sprintf("(a) %d taxon(a) with confirmed observations have NO niche-cache entry: %s",
                 length(missing_from_cache), paste(missing_from_cache, collapse = ", ")))
  }
  if (length(extra_in_cache) > 0) {
    fail(sprintf("(f) niche cache has %d key(s) with no matching taxon in the current roster (stale/renamed?): %s",
                 length(extra_in_cache), paste(extra_in_cache, collapse = ", ")))
  }
  ok(sprintf("(a)/(f) niche cache (%d keys) exactly matches the %d-taxon roster", length(cache_keys), length(roster)))

  for (site in sites) {
    site_obs <- niches[niches$Area_or_Site == site, ]

    # (c) every observation within a plausible distance of its own site's centroid
    lat0 <- median(site_obs$lat); lon0 <- median(site_obs$lon)
    dist_m <- sqrt(((site_obs$lat - lat0) * 111000)^2 +
                   ((site_obs$lon - lon0) * 111000 * cos(lat0 * pi / 180))^2)
    if (any(dist_m > DIST_THRESHOLD_M)) {
      bad <- which(dist_m > DIST_THRESHOLD_M)
      fail(sprintf("(c) %s has %d observation(s) > %dm from the site's own median coordinate (max %.0fm) -- likely a mislabeled/corrupted row, same class of bug as the Saloya outliers.",
                   site, length(bad), DIST_THRESHOLD_M, max(dist_m)))
    }

    # (d) footprint area plausibility (same formula as init_colonization()'s own xDim/yDim basis)
    lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
    lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 * cos(mean(site_obs$lat) * pi / 180)
    area_ha <- (lat_range_m * lon_range_m) / 10000
    if (!is.finite(area_ha) || area_ha < AREA_MIN_HA || area_ha > AREA_MAX_HA) {
      fail(sprintf("(d) %s footprint area = %.2f ha, outside the plausible [%.1f, %.1f] ha range.",
                   site, area_ha, AREA_MIN_HA, AREA_MAX_HA))
    }

    # (j) landscape physical extent + implied tree count (the values
    # init_colonization() and build_forest() actually use, post-2026-09-10)
    land_bbox <- site_landscape_bbox(site, obs_csv)
    if (is.null(land_bbox)) {
      fail(sprintf("(j) %s: site_landscape_bbox() returned NULL -- fewer than %d raw observations with finite lat/lon. init_colonization() would fall back to the species-filtered extent (the MindoMirador/Saloya collapse mode).",
                   site, 2L))
    }
    if (land_bbox$n < LAND_MIN_RAW_OBS) {
      fail(sprintf("(j) %s has only %d raw observations to define its landscape extent (minimum %d).",
                   site, land_bbox$n, LAND_MIN_RAW_OBS))
    }
    land_lat_m <- diff(land_bbox$lat) * 111000
    land_lon_m <- diff(land_bbox$lon) * 111000 * cos(mean(land_bbox$lat) * pi / 180)
    land_area_ha <- (land_lat_m * land_lon_m) / 10000
    land_ntree <- max(1L, round(land_area_ha * default_forestparams()$stems_per_ha))
    if (!is.finite(land_area_ha) || land_area_ha < LAND_AREA_MIN_HA || land_area_ha > AREA_MAX_HA) {
      fail(sprintf("(j) %s landscape area = %.2f ha (from %d raw obs), outside the plausible [%.1f, %.1f] ha range -- a real forest survey spans several hectares; this looks truncated.",
                   site, land_area_ha, land_bbox$n, LAND_AREA_MIN_HA, AREA_MAX_HA))
    }
    if (land_ntree < LAND_MIN_TREES) {
      fail(sprintf("(j) %s landscape would generate only %d trees (%.2f ha x %d stems/ha) -- below the minimum of %d.",
                   site, land_ntree, land_area_ha, default_forestparams()$stems_per_ha, LAND_MIN_TREES))
    }
    filt_lat_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
    filt_lon_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 * cos(mean(site_obs$lat) * pi / 180)
    filt_area_ha <- (filt_lat_m * filt_lon_m) / 10000
    if (is.finite(filt_area_ha) && filt_area_ha < 0.25 * land_area_ha) {
      message(sprintf("  NOTE: (j) %s species-filtered footprint (%.2f ha, %d obs) is far smaller than the raw landscape (%.2f ha, %d obs) -- init_colonization() correctly uses the raw extent; noting the divergence.",
                      site, filt_area_ha, nrow(site_obs), land_area_ha, land_bbox$n))
    }

    # (e) site_elev finite
    elev_source <- read.csv(obs_csv)
    elev_source <- elev_source[!is.na(elev_source$lat) & !is.na(elev_source$lon), ]
    elev_source <- augment_elevation(elev_source)
    site_elev <- mean(elev_source$Elevation_final_m[elev_source$Area_or_Site == site], na.rm = TRUE)
    if (!is.finite(site_elev)) {
      fail(sprintf("(e) %s has no finite site_elev (no observation there has lat/lon, or DTM lookup failed for all of them).", site))
    }

    # (g) hour-1 valid-ERA5-cell count -- NOTE ONLY (2026-09-10). This is
    # microclimf::runpointmodela()'s own per-cell test, which inspects only
    # the FIRST hourly ERA5 timestep; for ERA5 accumulated fields
    # (radiation, precipitation) there is no valid accumulation at hour 1,
    # so this count structurally undercounts (MindoMirador: 1/6 at hour 1
    # vs 6/6 on the 12-timestep yearly sample). The hard gate is (g2)
    # below, on the trustworthy monthly-sampled count. A hard gate on (g)
    # fires falsely and trains us to override it.
    microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site))
    if (file.exists(microenv_path)) {
      microenv <- readRDS(microenv_path)
      if (is.null(microenv$.n_valid_era5_cells)) {
        message(sprintf("  WARN: (g) %s's manifest predates .n_valid_era5_cells -- cannot check; re-run to get coverage.", site))
      } else {
        nv1 <- microenv$.n_valid_era5_cells
        nvm <- microenv$.n_valid_era5_cells_monthly
        if (!is.null(nvm) && !is.na(nv1) && nv1 < nvm) {
          message(sprintf("  NOTE: (g) %s hour-1 valid ERA5 cells = %d, but the 12-timestep yearly sample (g2) says %d -- expected divergence for ERA5 accumulated fields; gating on (g2).",
                          site, nv1, nvm))
        } else if (!is.na(nv1) && nv1 < MIN_VALID_ERA5_CELLS) {
          message(sprintf("  NOTE: (g) %s hour-1 valid ERA5 cells = %d (of %d); see (g2).",
                          site, nv1, microenv$.n_total_era5_cells))
        }
      }
      # (g2) same floor, but against the monthly-sampled count (2026-09-08)
      # -- run_microclimate_site.R's own independent check, sampling 12
      # across-the-year timesteps and requiring all finite, rather than
      # trusting microclimf::runpointmodela()'s hour-1-only test (see
      # MIN_VALID_ERA5_CELLS's own comment above). Catches a cell the hour-1
      # test wrongly passed (fine at hour 1, gappy later) that (g) alone
      # would miss.
      if (is.null(microenv$.n_valid_era5_cells_monthly)) {
        message(sprintf("  WARN: (g2) %s's manifest predates .n_valid_era5_cells_monthly -- cannot check; re-run to get coverage.", site))
      } else if (microenv$.n_valid_era5_cells_monthly < MIN_VALID_ERA5_CELLS) {
        fail(sprintf("(g2) %s has only %d ERA5 grid cell(s) valid across a 12-timestep yearly sample (of %d) -- below the minimum of %d, even though microclimf's own hour-1-only test may have passed it.",
                     site, microenv$.n_valid_era5_cells_monthly, microenv$.n_total_era5_cells, MIN_VALID_ERA5_CELLS))
      }

      # (h) per-pixel voxel_quantiles at the lowest height tier is at least
      # MIN_FINITE_FRAC finite for temp -- the direct downstream symptom of
      # (g)'s cause, checkable on every manifest regardless of vintage.
      h0_path <- file.path(microenv$.height_dir, sprintf("h%.2f.rds", microenv$.heights[1]))
      if (file.exists(h0_path)) {
        h0 <- readRDS(h0_path)
        q <- h0$voxel_quantiles$quantiles[["annual_both_temp"]]
        frac <- if (is.null(q)) NA_real_ else mean(is.finite(q))
        if (is.na(frac)) {
          message(sprintf("  WARN: (h) %s's height file has no annual_both_temp quantile entry -- cannot check.", site))
        } else if (frac < MIN_FINITE_FRAC) {
          fail(sprintf("(h) %s's per-pixel temp array is only %.2f%% finite (height tier 1) -- below the minimum of %.0f%%. Same failure signature as MindoMirador's collapse (0.03%% finite, vs. its own healthy-run baseline of 51%%).",
                       site, frac * 100, MIN_FINITE_FRAC * 100))
        }
      } else {
        message(sprintf("  WARN: (h) %s's height_dir has no h%.2f.rds -- cannot check.", site, microenv$.heights[1]))
      }

      # (k) the manifest's spatial extent (.spatial$ext, from the DTM the
      # microclimate run actually used) must ENCLOSE this site's model
      # domain -- the raw-observation bounding box + the standard 0.15deg
      # buffer. 2026-09-10: MindoMirador's microclimate collapsed because a
      # filename-keyed terrain cache returned a DTM for an earlier
      # (pre-grid-snap) domain while ERA5 was pulled for the shifted one;
      # the added ERA5 cells fell off the stale DTM. The grid snap is now
      # removed and get_dtm() etc. have an extent guard, but this asserts
      # the on-disk result directly, with the margin reported.
      if (!is.null(microenv$.spatial) && !is.null(microenv$.spatial$ext)) {
        raw <- read.csv(obs_csv)
        rs <- raw[!is.na(raw$Area_or_Site) & raw$Area_or_Site == site, c("lat", "lon")]
        rs$lat <- suppressWarnings(as.numeric(rs$lat)); rs$lon <- suppressWarnings(as.numeric(rs$lon))
        rs <- rs[is.finite(rs$lat) & is.finite(rs$lon), ]
        if (nrow(rs) >= 2) {
          BUF <- 0.15
          dom <- c(min(rs$lon) - BUF, max(rs$lon) + BUF, min(rs$lat) - BUF, max(rs$lat) + BUF)
          se <- as.numeric(microenv$.spatial$ext)  # xmin xmax ymin ymax
          tol <- 0.01
          gap_km <- c(S = (se[3] - dom[3]) * 111.32, N = (dom[4] - se[4]) * 111.32,
                      W = (se[1] - dom[1]) * 111.32, E = (dom[2] - se[2]) * 111.32)
          message(sprintf("  (k) %s manifest terrain vs domain: gap S/N/W/E = %+.2f/%+.2f/%+.2f/%+.2f km (>0 = terrain short of domain)",
                          site, gap_km["S"], gap_km["N"], gap_km["W"], gap_km["E"]))
          if (se[1] > dom[1] + tol || se[2] < dom[2] - tol ||
              se[3] > dom[3] + tol || se[4] < dom[4] - tol) {
            fail(sprintf("(k) %s's manifest terrain extent [%.4f,%.4f,%.4f,%.4f] does NOT enclose its model domain [%.4f,%.4f,%.4f,%.4f] (raw obs bbox + %.2fdeg). Stale terrain cache -- re-run Phase A.",
                         site, se[1], se[2], se[3], se[4], dom[1], dom[2], dom[3], dom[4], BUF))
          }
        }
      }

      # (i) manifest is not stale relative to its own height files
      # (2026-09-10 -- Saloya and LaElenita's Phase A jobs completed the
      # full per-height loop, then exited `COMPLETED 0:0` WITHOUT running
      # the manifest-save step: no error, exit 0, but the on-disk manifest
      # still pointed at an earlier run's data while the scratch height
      # files had been regenerated. A silent success is worse than a
      # failure -- it survives every content check because the content it
      # points at is internally consistent, just not this run's.). Two
      # assertions: every height file the manifest names exists, and the
      # manifest file is at least as new as the newest height file it
      # indexes. A manifest older than its own height data was written
      # before those files and never refreshed.
      hfiles <- file.path(microenv$.height_dir, sprintf("h%.2f.rds", microenv$.heights))
      missing_h <- hfiles[!file.exists(hfiles)]
      if (length(missing_h) > 0) {
        fail(sprintf("(i) %s's manifest names %d height file(s) that do not exist on disk (e.g. %s) -- manifest and height_dir are out of sync.",
                     site, length(missing_h), basename(missing_h[1])))
      }
      man_mtime <- file.info(microenv_path)$mtime
      newest_h  <- max(file.info(hfiles)$mtime, na.rm = TRUE)
      if (is.finite(newest_h) && man_mtime < newest_h - 60) {
        fail(sprintf("(i) %s's manifest (%s) is OLDER than its newest height file (%s) -- the manifest-save step did not run against this height data (the 2026-09-10 Saloya/LaElenita 'COMPLETED 0:0 without writing manifest' failure mode). Re-run the assembly.",
                     site, format(man_mtime, "%Y-%m-%d %H:%M"), format(newest_h, "%Y-%m-%d %H:%M")))
      }
    } else {
      message(sprintf("  WARN: (g)/(h) no microenv manifest found for %s at %s -- cannot check ERA5 coverage.", site, microenv_path))
    }
  }
  ok(sprintf("(c)/(d)/(e) every one of %d sites: observations within %dm of their own centroid, footprint area in [%.1f,%.1f]ha, finite elevation",
             length(sites), DIST_THRESHOLD_M, AREA_MIN_HA, AREA_MAX_HA))
  ok(sprintf("(g)/(h) every site with a checkable manifest: >=%d valid ERA5 cells, >=%.0f%% finite per-pixel temp array",
             MIN_VALID_ERA5_CELLS, MIN_FINITE_FRAC * 100))

  message("SANITY GATE PASSED.")
  invisible(TRUE)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cache_path <- if (length(args) >= 1) args[1] else NICHE_CACHE_PATH
  run_sanity_gate(cache_path)
}
