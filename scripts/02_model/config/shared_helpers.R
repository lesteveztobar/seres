# shared_helpers.R
# Small helpers factored out of duplicated logic across scripts/02_model/,
# found while extending the FOL-style refactor of get_colonization.R to the
# rest of the pipeline (2026-08). Each consolidates code that was previously
# copy-pasted verbatim (or near-verbatim) across 2-4 files.
# Sourced by: run/run_colonization.R, diagnostics/check_colonization_run.R,
# analysis/summarize_all_results.R, analysis/competition_analysis.R,
# plots/plot_functions.R, resolution/resolution_diagnostics.R,
# resolution/height_resolution_experiment.R.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

# ── Result-shape classification ────────────────────────────────────────────
# A saved colonization result (colonization_<site>_<tag>.rds) is one of two
# shapes: a sweep/factorial data.frame (one row per replicate x year x swept-
# parameter combination; produced by run_experiment()/run_factorial_
# experiment()), or a run_replicated() list(runs=, summary=). This computes
# the final-year combination/persistence summary each shape needs, without
# prescribing how the caller prints it -- check_colonization_run.R and
# summarize_all_results.R want visibly different text (one verbose, one
# compact-per-line for a multi-result survey), but were duplicating the
# exact same swept-column detection + aggregate() + persisted logic to get
# there. Callers format their own text from the fields returned here.
.classify_result_shape <- function(result) {
  if (is.data.frame(result)) {
    fixed_cols <- c("rep", "t", "totalS", "totalJ", "totalA", "total", "extinct")
    swept <- setdiff(names(result), fixed_cols)
    if (length(swept) == 0 || nrow(result) == 0) {
      return(list(shape = "sweep", degenerate = TRUE, swept = swept))
    }
    t_max <- max(result$t)
    final <- result[result$t == t_max, ]
    form  <- as.formula(paste("cbind(extinct, total) ~", paste(swept, collapse = " + ")))
    combo <- aggregate(form, data = final, FUN = mean)
    combo$persisted <- !as.logical(round(combo$extinct))
    list(shape = "sweep", degenerate = FALSE, swept = swept, t_max = t_max,
         combo = combo, persisting = combo[combo$persisted, ])
  } else {
    runs <- if (is.list(result) && !is.null(result$runs)) result$runs else list(result)
    n_ok <- sum(!vapply(runs, is.null, logical(1)))
    if (n_ok == 0 || is.null(result$summary)) {
      return(list(shape = "replicated", n_runs = length(runs), n_ok = n_ok, summary = NULL))
    }
    final_t <- max(result$summary$t)
    final <- result$summary[result$summary$t == final_t, ]
    # n_recruiting: how many replicates saw ANY seedling/juvenile presence
    # in the second half of the run (`recruited`, get_colonization.R,
    # added 2026-08-28) -- distinct from n_ok/persisted, which only track
    # whether a replicate ran and whether adults survived. A replicate can
    # be "persisting" (adults present) with n_recruiting excluding it
    # entirely -- i.e. it's just the founder cohort decaying with zero
    # replacement, not a self-sustaining population.
    n_recruiting <- if ("recruited" %in% names(result$summary)) {
      length(unique(result$summary$rep[result$summary$recruited]))
    } else {
      NA_integer_
    }
    list(shape = "replicated", n_runs = length(runs), n_ok = n_ok,
         n_recruiting = n_recruiting,
         final_t = final_t, final = final, summary = result$summary)
  }
}

# ── Forest structure defaults ──────────────────────────────────────────────
# Literature-calibrated forestparams (Myster 2017), previously copy-pasted
# identically in run/run_colonization.R, resolution/resolution_diagnostics.R,
# and resolution/height_resolution_experiment.R.
default_forestparams <- function() {
  list(
    stems_per_ha          = 298,    # Myster (2017) Table 5, mean of 4 primary MR plots
    mean_hgt              = 8.4,    # m — mean canopy height at this elevation
    sd_hgt                = 3.5,
    mean_crown_r          = 2.0,    # m — crown radius
    sd_crown_r            = 0.8,
    trunk_r               = 0.114,  # m — mean dsh 22.7 cm / 2 (Myster 2017)
    branch_density        = 3.0,    # m² branch surface per m² projected crown area
    epiphyte_footprint_m2 = 0.02    # m² bark area per Maxillariinae individual
  )
}

# ── Significance stars ─────────────────────────────────────────────────────
# Standard 3-level significance-star formatting for an (adjusted) p-value,
# previously duplicated in plot_functions.R, resolution/height_resolution_
# experiment.R (dunn_res$P.adjusted), and analysis/competition_analysis.R
# (result_df$ks_p, which additionally maps NA to "" rather than leaving it
# NA -- preserved here via `na_str`, defaulting to NA_character_ to match
# the other two call sites' original behavior exactly).
.sig_stars <- function(p, na_str = NA_character_) {
  ifelse(is.na(p), na_str,
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,  "*", "ns"))))
}

# ── Species color palette ──────────────────────────────────────────────────
# Per-species plotting colors: a single mid-palette color for 1 species (so
# it doesn't default to the palette's extreme end), else one color per
# species across the full palette range. Previously duplicated 4x in
# plot_functions.R (.plot_live, plot_abundance, plot_3d_abundance,
# plot_3d_abundance_animated).
.species_colors <- function(n_species) {
  if (n_species == 1) {
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  } else {
    scico::scico(n_species, palette = "lipari", begin = 0.2, end = 0.8)
  }
}

# ── Confirmed (non-default-fallback) species x site membership ─────────────
# load_observations() (paths.R) copies Identification into FinalID, then
# hardcodes any blank/NA FinalID to "Maxillaria acutifolia" -- a real
# identification and a blank-ID default are indistinguishable in FinalID
# alone. load_observations()'s own output still carries the RAW,
# unmodified Identification column though (only FinalID gets the
# substitution) -- so this filters on that raw column instead, to find
# sites where a species was actually, explicitly identified rather than
# defaulted. Added 2026-08-25 for plot_species_across_sites() (plot_
# functions.R), which must not compare a site's "acutifolia" curve against
# another site's if one of them is really just "nobody identified this."
.confirmed_species_sites <- function(species_name, niches = load_observations()) {
  hit <- !is.na(niches$Identification) & nzchar(trimws(niches$Identification)) &
    niches$Identification == species_name
  sort(unique(niches$Area_or_Site[hit]))
}

# ── Elevation helpers ───────────────────────────────────────────────────────
# Formerly elevation_helpers.R (merged into climate_variation_test.R
# 2026-07-28), relocated here 2026-08-27 once a second script (climate_
# variation_between_sites.R) needed them too -- this is exactly the
# duplication shared_helpers.R exists to catch.
#
# Per-observation elevation, combining field-recorded Elevation_m (present
# for Maquipucuna/Mashpi/MindoMirador/LaElenita/Saloya) with digital-
# elevation-model extraction at each observation's exact lat/lon (fills
# MindoTarabita/Yanayacu, which have NO recorded elevation at all).
#
# Reuses the per-site DTM already downloaded for microclimate modelling
# (data/raw/<site>/dtm/dtm.tif, cached by get_dtm() in lib.R, already
# reprojected to EPSG:4326 there) -- no new data acquisition needed, every
# site already has one.
#
# Does NOT modify combinedv3.csv (a manually curated, merged dataset -- see
# README's pipeline diagram) -- augment_elevation() returns an augmented
# copy of whatever data frame you pass it. Field-recorded values are kept
# as-is (not overwritten by the DTM) since they may reflect a field
# GPS/altimeter reading more precise than a resampled DEM pixel; the DTM is
# only used to fill in what the field record doesn't have.

# DTM-extracted elevation (m) at each row's (lon, lat) for one site.
.dtm_elevation <- function(site_name, lon, lat, raw_dir = RAW_DIR) {
  dtm_path <- file.path(raw_dir, site_name, "dtm", "dtm.tif")
  if (!file.exists(dtm_path)) {
    warning("No DTM for ", site_name, " at ", dtm_path, " -- returning NA elevation")
    return(rep(NA_real_, length(lon)))
  }
  dtm <- terra::rast(dtm_path)
  pts <- terra::vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
  if (!terra::same.crs(pts, dtm)) pts <- terra::project(pts, terra::crs(dtm))
  vals <- terra::extract(dtm, pts)
  as.numeric(vals[[2]])  # column 1 is the auto ID, column 2 is the DTM's single band
}

# Adds Elevation_dtm_m (always DTM-derived, for auditing/comparison against
# the field record) and Elevation_final_m (field-recorded Elevation_m where
# present and numeric, DTM-derived otherwise -- the column downstream
# analyses should actually use) to a niches-style data frame. Requires
# Area_or_Site, lon, lat columns; Elevation_m is optional (treated as
# entirely missing if absent).
augment_elevation <- function(niches_df, raw_dir = RAW_DIR) {
  niches_df$Elevation_dtm_m <- NA_real_
  for (site_name in unique(niches_df$Area_or_Site)) {
    idx <- which(niches_df$Area_or_Site == site_name)
    if (length(idx) == 0) next
    niches_df$Elevation_dtm_m[idx] <- .dtm_elevation(
      site_name, niches_df$lon[idx], niches_df$lat[idx], raw_dir)
  }

  field_elev <- if ("Elevation_m" %in% names(niches_df)) {
    suppressWarnings(as.numeric(niches_df$Elevation_m))
  } else {
    rep(NA_real_, nrow(niches_df))
  }
  niches_df$Elevation_final_m <- ifelse(!is.na(field_elev), field_elev, niches_df$Elevation_dtm_m)
  niches_df
}

# ── Per-site, per-height climate series (shared build step) ────────────────
# The data-loading core of climate_variation_test.R's within-site test --
# extracted 2026-08-27 so climate_variation_between_sites.R (cross-site,
# elevation-matched comparisons) can reuse it instead of re-reading every
# microenv_<site>_h0.40.rds from scratch. Returns, per site with a usable
# microenv on disk: its full long-format per-height-tier hourly climate
# series (site_pixels[[site]]: one row per height tier x hourly timestep,
# columns height/temp/relhum/swdown -- not actually per-pixel despite the
# name inherited from the original script; lookup_climate_by_height()
# already spatially averages each height's raster at write time, see
# get_colonization.R) and its mean recorded/DTM-filled elevation
# (site_elev[site], from augment_elevation()'s Elevation_final_m).
#
# Requires get_colonization.R (microenv_heights(), lookup_climate_by_height())
# to already be sourced by the caller.
.build_site_climate_series <- function(sites = NULL, processed_dir = PROCESSED_DIR,
                                        niches = NULL) {
  if (is.null(niches)) {
    niches <- load_observations()
    niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                     !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
    niches <- augment_elevation(niches)
  }
  if (is.null(sites)) sites <- sort(unique(niches$Area_or_Site))

  site_pixels <- list()
  site_elev   <- numeric(0)
  for (site in sites) {
    microenv_path <- file.path(processed_dir, sprintf("microenv_%s_h0.40.rds", site))
    if (!file.exists(microenv_path)) {
      message("Skipping ", site, " -- no ", microenv_path)
      next
    }
    message("Reading ", site, "...")
    microenv <- readRDS(microenv_path)
    heights  <- microenv_heights(microenv)

    n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))
    px <- do.call(rbind, parallel::mclapply(heights, function(h) {
      cl <- lookup_climate_by_height(h, microenv)
      if (is.null(cl)) return(NULL)
      data.frame(height = h, temp = cl$temp, relhum = cl$relhum, swdown = cl$swdown)
    }, mc.cores = n_cores))
    site_pixels[[site]] <- px
    site_elev[site] <- mean(niches$Elevation_final_m[niches$Area_or_Site == site], na.rm = TRUE)
  }
  list(site_pixels = site_pixels, site_elev = site_elev)
}
