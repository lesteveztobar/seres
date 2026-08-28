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
# Usage: Rscript scripts/02_model/setup/characterize_niches.R [height_step]
#   height_step defaults to 0.4 (the production resolution as of
#   2026-08-17, was 0.25) — pass 0.1 to use the unsuffixed manifests
#   instead, if you have those for every site.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

args <- commandArgs(trailingOnly = TRUE)
HEIGHT_STEP <- if (length(args) >= 1) as.numeric(args[1]) else 0.4
manifest_suffix <- if (HEIGHT_STEP != 0.1) sprintf("_h%.2f", HEIGHT_STEP) else ""

VARS <- NICHE_VARS  # c("temp", "relhum", "swdown") — see get_colonization.R

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
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

  # Footprint: raster pixels actually touched by this site's own field
  # observations -- the "landscape" this script has to work with, since
  # (unlike init_colonization()) there's no simulated grid here. Falls back
  # to NULL (pooled whole-raster quantiles) when microenv$.spatial isn't a
  # usable plain-vector extent -- the current state of every existing
  # manifest (see .spatial_extent_usable()).
  footprint <- if (.spatial_extent_usable(microenv)) {
    px <- .lonlat_to_pixel(obs_site$lon, obs_site$lat, microenv$.spatial)
    unique(data.frame(row = px$row, col = px$col))
  } else {
    NULL
  }

  cc_voxel <- build_clim_cache_voxel(microenv, footprint = footprint)

  bg_rows[[s]] <- voxel_background_table(cc_voxel, microenv, footprint, VARS)

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
message(sprintf("Background pool: %d pixel-height-months across %d site(s)", nrow(bg_df), length(bg_rows)))

bg_density <- build_background_density(bg_df, VARS)
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
