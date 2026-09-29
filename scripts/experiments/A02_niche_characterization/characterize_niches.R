# characterize_niches.R
# Precomputes each species' realized climate niche (temp/RH/light) pooled
# across EVERY site where it was observed, not just the one site.rds file
# happens to be modeling — a species with 2-3 records at one site may have
# several more at others. Saves data/processed/species_niches.rds, loaded by
# init_colonization() (get_colonization.R) in preference to the this-site-
# only get_niche_voxel() fallback. Also saves data/processed/niche_background_density.rds
# (the pooled reference distribution every species is scored against), reused
# by check_niche_suitability.R / plot_niche_suitability() for diagnostics.
#
# Presence and background climate both come from the per-voxel quantile
# cache (build_clim_cache_voxel()/voxel_climate_table()/
# voxel_background_table() in get_colonization.R) rather than a spatially-
# flattened per-height mean — real per-pixel resolution activates
# automatically once microenv$.spatial carries a usable plain-vector extent
# (see .spatial_extent_usable()); every existing manifest degrades
# gracefully to pooled (whole-raster) quantiles today, same as everywhere
# else in the per-voxel system.
#
# Each observed individual is augmented with the 12 within-year monthly
# conditions at its recorded lon/lat/height, rather than collapsed to that
# height's single annual value — a species realistically tolerates whatever
# seasonal range its occupied height experiences over a year, not just that
# height's average, and pooling the full monthly record turns each raw
# observation into 12 data points instead of 1 (helpful for the kernel
# density fit, not just the old box), giving a far better-supported niche
# estimate than the handful of field observations alone could. swdown
# (light) uses only daytime quantiles; temp/relhum use the pooled day+night
# quantile — see get_clim_voxel()'s `both_vars` in get_colonization.R.
#
# The background reference for a site is every (footprint pixel x height
# tier x month) combination touched by that site's own field observations
# — "available but not necessarily occupied" conditions at the pixels
# actually surveyed, pooled once per site and then across all sites so
# every species niche is scored against the same yardstick.
#
# Rerun this whenever you add observations to OBSERVATIONS_CSV (paths.R;
# data/csv/combined_with_identification.csv by default, overridable via
# CANOPY_OBS_CSV) — it's the only step that needs to change; every
# colonization run downstream picks up the refined niches automatically
# without needing to recompute anything itself.
#
# Requires microenv_<site>[_h<step>].rds to already exist for every site with
# observations (run_microenv.sh) — reads each site's climate once.
#
# Usage: Rscript scripts/experiments/A02_niche_characterization/characterize_niches.R [height_step]
#   height_step defaults to 0.4 (the production resolution as of
#   2026-08-17, was 0.25) — pass 0.1 to use the unsuffixed manifests
#   instead, if you have those for every site.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")  # .filter_maxillariinae() -- 2026-08-29
source("scripts/02_model/engine/get_colonization.R")

args <- commandArgs(trailingOnly = TRUE)
HEIGHT_STEP <- if (length(args) >= 1) as.numeric(args[1]) else 0.4
manifest_suffix <- if (HEIGHT_STEP != 0.1) sprintf("_h%.2f", HEIGHT_STEP) else ""

VARS <- NICHE_VARS  # c("temp", "relhum", "swdown") — see get_colonization.R

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
# Check 2 (2026-09-02): Saloya excluded from the background/niche pool here --
# coordinate outlier corrupts its contribution to the pooled background
# density. Comparison-only variant, saved to separate CANOPY_NICHE_CACHE/
# CANOPY_NICHE_BACKGROUND paths -- does not touch the production cache.
if (Sys.getenv("EXCLUDE_SALOYA", unset = "0") == "1") {
  niches <- niches[niches$Area_or_Site != "Saloya", ]
  message("EXCLUDE_SALOYA=1 -- Saloya dropped from the background/niche pool")
}

# 2026-09-02 (v7 rebuild, Phase 1.7): exclude the SAME held-out individuals
# runcolonization() (get_colonization.R) validates against
# (get_held_out_split(), shared_helpers.R -- one fixed split, not a fresh
# random draw per caller). Before this, every observation was pooled into
# the niche cache regardless of whether it was later held out for
# validation, so a held-out individual's own presence record could
# contribute to the very niche model it was then scored against (Task 3's
# circularity finding, methods_update_report.md). This also automatically
# excludes held-out individuals' own pixels from each site's background
# footprint below (`obs_site` is derived from this already-filtered
# `niches`), not just from the presence table.
held_out <- get_held_out_split(niches)
message(sprintf(
  "Held-out validation split: excluding %d of %d observations from niche estimation (train_frac=0.70, matches runcolonization()'s split).",
  sum(held_out), length(held_out)))
niches <- niches[!held_out, ]

sites <- sort(unique(niches$Area_or_Site))

# ── Pass 1: per-observation ("presence") climate via voxel_climate_table(),
# and this site's background via voxel_background_table() (both
# get_colonization.R) ─────────────────────────────────────────────────────
obs_rows <- list()
bg_rows  <- list()

for (s in sites) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", s, manifest_suffix))
  if (!file.exists(microenv_path)) {
    message("Skipping ", s, " -- no microenv at ", microenv_path, " (run_microenv.sh first)")
    next
  }
  message("Reading climate for ", s, "...")
  microenv <- readRDS(microenv_path)
  obs_site <- niches[niches$Area_or_Site == s, ]

  # Footprint for the BACKGROUND (availability) distribution and the voxel
  # cache: the model's ACTUAL simulated landscape -- the raw-observation
  # bounding box rasterised (site_landscape_bbox(), the SAME extent
  # init_colonization() derives xDim/yDim from), NOT the handful of pixels
  # touched by identification-confirmed observations.
  #
  # 2026-09-10 fix. The niche model is a presence-vs-background kernel-
  # density ratio (niche_density_model(), methods.tex Eq. nicheratio).
  # Previously `footprint` came from `obs_site` (species-filtered,
  # held-out-excluded observations) -- 3-8 raster pixels per site -- so the
  # availability distribution was centred on presence locations, the
  # presence/background ratio was compressed toward 1, every niche came out
  # less selective than it is, and it was not the availability set the
  # model actually samples (run_pass2_establish() scores every voxel in the
  # 7-33 ha landscape against this niche). Same species-filter leak as the
  # landscape-bbox and canopy-ceiling bugs. See
  # docs/methods_update_report.md "Audit: site-level physical quantities".
  #
  # Presence rows below still come from `obs_site` -- presence IS the
  # confirmed-observation set, by definition. Only the background/footprint
  # changes. Falls back to NULL (pooled whole-raster quantiles) when
  # microenv$.spatial isn't a usable plain-vector extent.
  footprint <- if (.spatial_extent_usable(microenv)) {
    bb <- site_landscape_bbox(s)
    if (is.null(bb)) bb <- list(lat = range(obs_site$lat), lon = range(obs_site$lon))
    corners <- .lonlat_to_pixel(rep(bb$lon, 2L), c(bb$lat, rev(bb$lat)), microenv$.spatial)
    row_lo <- max(1L, min(corners$row) - 2L)
    row_hi <- min(microenv$.spatial$nrow, max(corners$row) + 2L)
    col_lo <- max(1L, min(corners$col) - 2L)
    col_hi <- min(microenv$.spatial$ncol, max(corners$col) + 2L)
    expand.grid(row = row_lo:row_hi, col = col_lo:col_hi)
  } else {
    NULL
  }

  cc_voxel <- build_clim_cache_voxel(microenv, footprint = footprint)

  # 2026-09-06 (item 2, background pool weighting -- CORRECTED AGAIN, per
  # the author's explicit instruction to retract the pooled-per-site fix
  # below and replace it): two earlier fixes at this exact spot are both
  # retired now --
  #   (a) an original fix assumed footprint-pixel rows were numerically
  #       IDENTICAL duplicates and tried de-duplicating them -- WRONG,
  #       disproven directly (unique(bg_df) removed 0 of 25,788 rows: once
  #       .spatial_extent_usable() returns TRUE, voxel_climate_table()'s
  #       per-pixel lookup is genuinely spatially resolved).
  #   (b) a second fix requested the background in POOLED mode
  #       (footprint=NULL) for every site, so every site contributed
  #       exactly zDim x 12 rows regardless of its own pixel count -- this
  #       DID equalize site weight, but at the cost of replacing real
  #       per-pixel values with one whole-raster spatial average per
  #       height/month, narrowing the background distribution and
  #       systematically inflating every taxon's apparent selectivity.
  # Final fix: keep the REAL per-pixel background rows (cc_voxel, this
  # site's own real footprint, built above) -- full spatial variance
  # preserved -- and weight each row by 1/(this site's own footprint pixel
  # count), so every site still contributes equal TOTAL weight to the
  # pooled density (build_background_density()'s new `weight_col`) without
  # discarding per-pixel variance to get there. A site with 3 footprint
  # pixels and a site with 8 each end up mattering equally to the pooled
  # background; within a site, all pixels are weighted equally to each
  # other (no additional per-pixel ecological reason to prefer one).
  # 2026-09-08 BUGFIX: weighting each row by 1/n_fp (pixel count alone)
  # does NOT give equal total site weight -- each site contributes
  # n_fp*zDim*12 rows (one per pixel x height tier x month), so summing
  # 1/n_fp over all of them gives zDim*12, not 1 -- and zDim varies by
  # site (48-125 in this project), so a tall-canopy site still ended up
  # weighted ~2.6x a short one (confirmed directly: printed effective
  # weights of 576/828/1236/900/1500/1200 for the 6 sites -- exactly
  # zDim*12 each, not equal). Fixed by weighting each row by 1/(this
  # site's own total row count) instead, so every site's weights sum to
  # exactly 1 regardless of pixel count OR height-tier count.
  site_bg_table <- as.data.frame(voxel_background_table(cc_voxel, microenv, footprint, VARS))
  bg_rows[[s]] <- data.frame(.site = s, .weight = 1 / nrow(site_bg_table), site_bg_table)

  obs <- data.frame(lon = obs_site$lon, lat = obs_site$lat, height = obs_site$Height_m)
  presence <- voxel_climate_table(cc_voxel, microenv, obs, VARS)
  # voxel_climate_table()'s default months = 1:12 means each observation
  # contributes 12 consecutive monthly rows in order -- rep(..., each = 12)
  # tags them with the matching species.
  obs_rows[[s]] <- cbind(species = rep(obs_site$FinalID, each = 12), as.data.frame(presence))

  message(sprintf(
    "  %s: %d observations, %d footprint pixel(s)", s, nrow(obs_site),
    if (is.null(footprint)) 0L else nrow(footprint)
  ))
}

if (length(obs_rows) == 0) stop("No usable observations across any site with a microenv file.")
if (length(bg_rows) == 0) stop("No usable background climate across any site with a microenv file.")

obs_df <- do.call(rbind, obs_rows)
obs_df$species <- as.character(obs_df$species)

bg_df <- do.call(rbind, bg_rows)

# 2026-09-06 (item 2, background pool weighting -- FINAL fix, see the note
# where bg_rows[[s]] is built, above, for the full history/rationale):
# real per-pixel background rows are kept (no pooling-to-one-value-per-
# site), each carrying a `.weight` of 1/(that row's site's own footprint
# pixel count) so every site still contributes equal total weight.
site_weights <- tapply(bg_df$.weight, bg_df$.site, function(w) sum(w))
message(sprintf(
  "Background pool: %d real per-pixel rows across %d site(s) -- effective site weights (should be ~equal): %s",
  nrow(bg_df), length(bg_rows),
  paste(sprintf("%s=%.3f", names(site_weights), site_weights), collapse = ", ")))

bg_density <- build_background_density(bg_df, VARS, weight_col = ".weight")
bg_out_path <- NICHE_BACKGROUND_PATH
saveRDS(bg_density, bg_out_path)
message("Saved background density to ", bg_out_path)

# ── Pass 2: pooled niche model per species, across every site it was
# observed at ────────────────────────────────────────────────────────────────
species_ids <- sort(unique(obs_df$species))
niche_cache <- lapply(species_ids, function(sp) {
  clim_vals <- as.matrix(obs_df[obs_df$species == sp, VARS])
  n_obs     <- nrow(clim_vals) / 12L  # each field observation contributes 12 monthly rows
  n_sites   <- length(unique(niches$Area_or_Site[niches$FinalID == sp]))
  message(sprintf("%-30s %3d observations (%4d monthly data points) across %d site(s)",
                  sp, n_obs, nrow(clim_vals), n_sites))
  niche_density_model(clim_vals, bg_density, VARS)
})
names(niche_cache) <- species_ids

out_path <- NICHE_CACHE_PATH
saveRDS(niche_cache, out_path)
message(sprintf("\nSaved %d species niches to %s", length(niche_cache), out_path))
