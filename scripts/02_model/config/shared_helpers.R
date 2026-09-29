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
  fp <- .default_forestparams_base()
  # K calibration against literature stand density (2026-09-29): written by
  # scripts/diagnostics/calibrate_k.R; overrides the two K terms when present.
  # Target = 2,800 orchid stands/ha (Alzate-Q et al. 2019, Veracruz) x 271/1,348
  # (Ecuadorian Maxillariinae spp. / Mexican orchid spp., i.e. Mexican
  # Maxillariinae density 41/1,348 scaled up by 271/41) = 562.9 per ha at Maquipucuna.
  base <- if (exists("BASE_DIR")) BASE_DIR else "/home/s38leste_hpc/seres"
  cal_path <- Sys.getenv("CANOPY_K_CALIBRATION", unset = file.path(base, "data", "params", "k_calibration.rds"))
  if (nzchar(cal_path) && file.exists(cal_path)) {
    cal <- readRDS(cal_path)
    fp$occupiable_bark_fraction <- cal$occupiable_bark_fraction
    fp$maxillariinae_community_share <- cal$maxillariinae_community_share
  }
  fp
}

.default_forestparams_base <- function() {
  list(
    stems_per_ha          = 298,    # Myster (2017) Table 5, mean of 4 primary MR plots
    mean_hgt              = 8.4,    # m — mean canopy height at this elevation
    sd_hgt                = 3.5,
    mean_crown_r          = 2.0,    # m — crown radius
    sd_crown_r            = 0.8,
    trunk_r               = 0.114,  # m — mean dsh 22.7 cm / 2 (Myster 2017)
    branch_density        = 3.0,    # m² branch surface per m² projected crown area
    epiphyte_footprint_m2 = 0.02,   # m² bark area per Maxillariinae individual
    # ── Carrying-capacity terms added 2026-09-10 (K recalibration) ─────────
    # See build_forest() (get_colonization.R) and
    # docs/methods_update_report.md "Carrying capacity audit". The pre-2026-
    # 09-10 formula gave total K of 1.6-92 million individuals per site
    # (vs. 78 field observations across all 7 sites) for three reasons, all
    # now addressed: (1) the crown bark-area term was dimensionally wrong
    # (m³ not m²) and evaluated once per vertical tier instead of once per
    # tree — fixed in build_forest(); (2) no term for how much woody
    # surface is actually colonisable; (3) `total_occ` sums over species so
    # K represented the whole epiphyte community while the model contains
    # one subtribe.
    #
    # occupiable_bark_fraction: proportion of woody surface area actually
    # colonisable by an establishing epiphyte — excludes smooth bark,
    # sun-exposed faces, and wood too young/thin to retain a seedling.
    # An explicit modelling assumption (treated like sigma / p_germ):
    # stated as such in Methods and swept across an order of magnitude in
    # the sensitivity design. Central value 0.02.
    occupiable_bark_fraction   = 0.02,
    # maxillariinae_community_share: fraction of the vascular-epiphyte
    # community (by which `total_occ` capacity is implicitly shared) that
    # is Maxillariinae. Derivation: our POWO/WCVP parse gives 271
    # Maxillariinae of 4,355 Ecuadorian orchid species (6.22%); scaled by
    # the orchid share of vascular-epiphyte species (Zotz 2013,
    # Bot. J. Linn. Soc. 171: ~68% of ~27,600 vascular epiphyte species are
    # Orchidaceae) -> 0.0622 * 0.68 = 0.0423. This is a species-richness
    # proxy for the abundance share — a stated assumption, flagged in
    # Methods. Swept alongside occupiable_bark_fraction.
    maxillariinae_community_share = 0.0423
  )
}

# ── Site canopy ceiling: ONE shared definition (v7 rebuild, Phase 1.6) ─────
# Both the landscape's own height-tier count (zDim, set in
# run_microclimate_site.R from MAX measured CanopyHeight_m via
# height_ceiling()) and the dispersal-kernel wind-decay reference height
# (canopy_z, disperse()/get_colonization.R) need "this site's canopy top."
# Before this fix, that same MAX-measured-CanopyHeight_m-with-fallback
# formula was independently duplicated in 9 places (run_colonization.R plus
# 8 analysis scripts, several of which had drifted to MEAN instead of MAX
# before Phase 1.6) -- one shared implementation now, called everywhere.
# Physical justification for MAX specifically (not an arbitrary match to
# zDim): the Cionco (1972) exponential wind-decay profile disperse()
# implements is defined relative to the height at which wind is
# UNDECAYED -- the canopy top -- not a within-stand average, so MAX is the
# physically correct choice for canopy_z on its own terms, independent of
# also matching zDim's basis.
# fallback_max_height: caller-supplied (this site's own microenv/landscape
# max height tier), used only when no observation at this site has a
# measured CanopyHeight_m.
site_canopy_ceiling <- function(site_name, niches, fallback_max_height,
                                obs_csv = OBSERVATIONS_CSV) {
  # 2026-09-10: the canopy ceiling is a PHYSICAL property of the site --
  # the tallest tree recorded there -- not of which epiphytes were
  # confirmed. Previously this took max(CanopyHeight_m) over the
  # species-filtered `niches`, so at MindoTarabita (3 raw CanopyHeight_m
  # measurements, 21-30 m; only the one at 21 m sits on a
  # confirmed-Maxillariinae row) the ceiling came out 21.0 m instead of
  # 30 m -- capping mean_hgt / mean_crown_r / zDim-consistent geometry
  # well below the site's real canopy and below its own microclimate
  # landscape's height ceiling (29.7 m). Same root cause as the
  # landscape-bbox bug (see site_landscape_bbox()). Now: max over ALL raw
  # CanopyHeight_m at the site, and never below fallback_max_height (the
  # microclimate manifest's own height ceiling, itself the site's measured
  # or p99-vhgt canopy top).
  raw_ceiling <- suppressWarnings(max(niches$CanopyHeight_m[niches$Area_or_Site == site_name], na.rm = TRUE))
  csv <- tryCatch(utils::read.csv(obs_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(csv) && all(c("CanopyHeight_m", "Area_or_Site") %in% names(csv))) {
    ch <- suppressWarnings(as.numeric(csv$CanopyHeight_m[csv$Area_or_Site == site_name]))
    ch <- ch[is.finite(ch)]
    if (length(ch) > 0) raw_ceiling <- max(raw_ceiling, max(ch), na.rm = TRUE)
  }
  ceiling <- if (is.finite(raw_ceiling)) raw_ceiling else fallback_max_height
  if (is.finite(fallback_max_height)) ceiling <- max(ceiling, fallback_max_height)
  ceiling
}

# 2026-09-02 (v7 rebuild, Phase 1.5): site_forestparams() is now the
# UNCONDITIONAL default (the CANOPY_SITE_TREE_HEIGHT opt-in flag from the
# original Task 6 version is removed -- every site now gets its own
# mean_hgt/sd_hgt, always).
#
# mean_hgt/sd_hgt calibration (Phase 1.4): Lang et al.'s raster measures
# TOP-OF-CANOPY per pixel (a canopy height model); Myster (2017)'s 8.4m is
# the MEAN of every stem in a forest-plot inventory -- a stand-level
# statistic, not a canopy-top statistic. These are not the same quantity
# (see methods_update_report.md's "Task 6" item: every site's raw Lang
# pixel mean came out 3-4x Myster's own value, even at Maquipucuna itself,
# the site Myster actually measured). Calibrated at Maquipucuna, the one
# site with both quantities available: ratio = 8.4 / 13.7 (Myster's value /
# Maquipucuna's own measured CanopyHeight_m mean). That ratio is then
# applied to EVERY site's own canopy-top statistic: measured CanopyHeight_m
# mean where field data exists (Maquipucuna, Mashpi, MindoTarabita), else
# the TIGHT-buffer Lang et al. raster mean (fetch_vhgt_site.R -- median-
# centroid + 300m buffer, saved as vhgt_tight_forestparams.tif, NOT the
# wide-bbox production vhgt.tif) for MindoMirador/Yanayacu/Saloya/
# LaElenita.
#
# 2026-09-02 REVISED: the Lang-based calibration above was tried and
# rejected -- it doesn't hold up against the only check available (the 3
# sites with an actual field CanopyHeight_m mean). Measured/Lang-tight
# ratios come out 0.41 (Maquipucuna), 0.53 (Mashpi), 0.77 (MindoTarabita) --
# not a stable conversion factor, and the Lang route doesn't even reproduce
# the measured ORDERING: Lang gives Maquipucuna/Mashpi/MindoTarabita a
# 2m spread (33.19/35.19/35.25) against the measured sites' actual 13.3m
# spread (13.7/18.8/27.0), and it makes MindoMirador's forest taller than
# Saloya's while the max-based canopy ceiling says the opposite. The tight-
# buffer Lang extraction and its summary table (`output/vhgt_site_summary_
# tight.csv`) are KEPT as a reported diagnostic documenting why the raster
# was not used for parameterisation -- not deleted, just not the basis for
# mean_hgt/sd_hgt below.
#
# REPLACEMENT -- a ceiling-scaled rule, calibrated at Maquipucuna (the one
# site with both a real stand inventory AND a well-supported canopy
# ceiling): ratio = 8.4 (Myster's stand mean) / 19.0 (Maquipucuna's own
# MAX-based canopy ceiling, site_canopy_ceiling()) = 0.442. Applied to
# EVERY site's own ceiling -- the SAME per-site quantity that already
# drives zDim and canopy_z (Phase 1.6) -- so mean_hgt, sd_hgt, canopy_z,
# and zDim all derive from one shared number per site, not three
# independently-sourced ones. sd_hgt uses Myster's own coefficient of
# variation (3.5/8.4 = 0.417) applied to the resulting mean_hgt, rather
# than an independently-sourced spread.
#
# `ceiling`: this site's own site_canopy_ceiling() value -- passed in by
# the caller (already computed there for zDim/canopy_z) rather than
# recomputed here, so there is exactly one ceiling computation per call
# site, not two.
# ── Landscape physical extent (2026-09-10) ────────────────────────────────
# The colonization landscape is a PHYSICAL PLACE; its horizontal extent
# must come from every observation recorded at the site, NOT from the
# species-filtered (Maxillariinae + identified) subset that
# init_colonization() uses for the species list and niche lookups. At sites
# where few Maxillariinae were confirmed, deriving extent from the filtered
# set collapsed the landscape to the bounding box of a handful of points --
# MindoMirador: 21 raw observations spanning ~7 ha -> 3 filtered -> 0.40 ha,
# 120 trees, K=221; Saloya: 17 -> 5 -> 0.47 ha. That is a species-modelling
# choice leaking into the landscape geometry. This helper returns the raw
# lat/lon bounding box (matches canopy_audit.R's long-standing choice).
site_landscape_bbox <- function(site_name, obs_csv = OBSERVATIONS_CSV) {
  raw <- tryCatch(utils::read.csv(obs_csv, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (is.null(raw) || !all(c("lat", "lon", "Area_or_Site") %in% names(raw))) return(NULL)
  d <- raw[!is.na(raw$Area_or_Site) & raw$Area_or_Site == site_name, c("lat", "lon"), drop = FALSE]
  d$lat <- suppressWarnings(as.numeric(d$lat))
  d$lon <- suppressWarnings(as.numeric(d$lon))
  d <- d[is.finite(d$lat) & is.finite(d$lon), , drop = FALSE]
  if (nrow(d) < 2) return(NULL)
  list(lat = range(d$lat), lon = range(d$lon), n = nrow(d))
}

site_forestparams <- function(site_name, ceiling) {
  fp <- default_forestparams()
  if (!is.finite(ceiling)) {
    message("site_forestparams(): no finite canopy ceiling for ", site_name,
            " -- keeping literature default mean_hgt/sd_hgt (", fp$mean_hgt, "m/", fp$sd_hgt, "m).")
    return(fp)
  }
  ratio <- 8.4 / 19.0    # Myster stand mean / Maquipucuna's own ceiling
  cv    <- 3.5 / 8.4     # Myster's own sd/mean, applied to the scaled mean_hgt
  fp$mean_hgt <- ratio * ceiling
  fp$sd_hgt   <- cv * fp$mean_hgt
  # 2026-09-10: crown radius now scales with stature too. Previously
  # mean_crown_r was the fixed literature default (2.0 m) at every site
  # regardless of canopy ceiling (19-50 m), so crown bark area
  # (pi*r^2*branch_density) -- and therefore carrying capacity -- carried
  # no between-site signal (K/tree was 1.82-1.86 everywhere). Crown radius
  # scales roughly linearly with tree height in tropical forest allometry;
  # apply Myster's own crown-radius / stand-mean-height ratio (2.0/8.4) to
  # the site-scaled mean_hgt, so crown area ~ h^2. sd_crown_r keeps Myster's
  # own crown-radius CV (0.8/2.0).
  crown_ratio <- 2.0 / 8.4
  crown_cv    <- 0.8 / 2.0
  fp$mean_crown_r <- crown_ratio * fp$mean_hgt
  fp$sd_crown_r   <- crown_cv * fp$mean_crown_r
  message(sprintf(
    "site_forestparams(): %s ceiling=%.2fm -> mean_hgt=%.2fm sd_hgt=%.2fm mean_crown_r=%.2fm sd_crown_r=%.2fm",
    site_name, ceiling, fp$mean_hgt, fp$sd_hgt, fp$mean_crown_r, fp$sd_crown_r))
  fp
}

# ── Held-out validation split: ONE fixed split, shared by every consumer ──
# (v7 rebuild, Phase 1.7). Previously, runcolonization() (get_colonization.R)
# drew its own random 70/30 train/val split fresh on EVERY call -- since
# each replicate passes a different `seed`, every replicate of the same
# design actually validated against a DIFFERENT held-out set, and
# characterize_niches.R never excluded anything from the niche cache at
# all, so held-out individuals were routinely pooled into the very niche
# model they were later scored against (methods_update_report.md's Task 3
# circularity finding). This function computes ONE deterministic split
# (fixed seed, independent of any caller's own stochastic `seed`), split
# per (site, species) group so every species keeps its own train_frac
# fraction -- both characterize_niches.R and runcolonization() call this
# and get the identical answer for the same input `niches`, and it does
# not touch the caller's own RNG stream (restores .Random.seed on exit, so
# a stochastic colonization run's own random draws are unaffected by
# having called this first).
# Returns a logical vector aligned with `niches`'s own row order (TRUE =
# held out for validation, FALSE = used for spin-up/niche estimation).
get_held_out_split <- function(niches, train_frac = 0.70, seed = 42) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv))
  # This set.seed() no longer determines the split itself (see 2026-09-06
  # note below, in the loop) -- it exists only to guarantee .Random.seed is
  # present in .GlobalEnv before the loop's first get() call.
  set.seed(seed)
  held_out <- rep(FALSE, nrow(niches))
  key <- interaction(niches$Area_or_Site, niches$FinalID, drop = TRUE)
  groups <- split(seq_len(nrow(niches)), key)

  # 2026-09-06: real bug, root-caused and fixed -- a SINGLE set.seed(seed)
  # call before this loop meant every group's sample() draw consumed from
  # the SAME continuing RNG stream, in whatever order split() produced the
  # groups (alphabetical by the site.species interaction factor level).
  # This made the split depend on GROUP ORDER, not just on `seed` and group
  # membership: renaming "Sudamerlycaste sp." -> "Ida sp." moved that
  # group's alphabetical position (S... near the end -> I... near the
  # start), which shifted every subsequent group's position in the RNG
  # stream and changed WHICH individuals got held out in groups that were
  # never touched by the rename itself (confirmed: Maquipucuna's and
  # Mashpi's footprint pixel counts each changed by 1 between two v7
  # rebuilds with an identical held-out COUNT, 19/78 both times -- same
  # count, different set). Fixed by seeding EACH group independently from a
  # hash of its own (site, species) key, not from a shared incrementing
  # stream -- the split for any one group is now invariant to every other
  # group's existence, name, or position, so renaming/adding/reordering
  # taxa can never again change an unrelated group's held-out individuals.
  for (g in names(groups)) {
    idx <- groups[[g]]
    n_train <- max(1L, round(length(idx) * train_frac))
    group_seed <- seed + strtoi(substr(digest_key(g), 1, 7), base = 16L)
    old_seed2 <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(group_seed)
    # See the 2026-09-05 note (still applicable): sample(x, n) on a
    # length-1 numeric x samples from 1:x, not "the one candidate" -- guard
    # singleton groups explicitly.
    train_idx <- if (length(idx) == 1) idx else sample(idx, n_train)
    assign(".Random.seed", old_seed2, envir = .GlobalEnv)
    held_out[setdiff(idx, train_idx)] <- TRUE
  }
  held_out
}

# Deterministic, dependency-free string hash (R's own -- no need for
# digest:: package just for a stable per-group seed offset). Sums character
# codes with a simple multiplicative mix; not cryptographic, just needs to
# be stable and well-distributed across (site, species) key strings.
digest_key <- function(s) {
  codes <- utf8ToInt(s)
  # Reduce in double precision (not integer) to avoid 32-bit overflow, then
  # floor to an integer that fits set.seed()'s valid range at the end.
  h <- Reduce(function(acc, c) (acc * 31 + c) %% 2147483647, codes, 7)
  sprintf("%07x", as.integer(h %% 16777215))
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

  # 2026-09-05: site_elev is a physical property of the SITE, unrelated to
  # whether any individual observation there was ever identified to
  # species -- but it was being read off `niches`, which every caller
  # filters to identification-confirmed rows first (`!is.na(FinalID)`).
  # LaElenita has zero rows surviving that filter (Task 0), so
  # `mean(niches$Elevation_final_m[Area_or_Site=="LaElenita"], na.rm=TRUE)`
  # silently returned NaN -- LaElenita computed a perfectly good climate
  # profile here (site_pixels comes from the microenv, independent of
  # `niches` entirely) but got dropped from any elevation-based test
  # downstream (e.g. Task 1 Part B) purely because of this, not because it
  # genuinely lacks a location. Fixed by deriving site_elev from a
  # SEPARATELY, minimally-filtered (lat/lon present only, no
  # identification/height requirement) elevation source, independent of
  # whatever filtering the caller applied to `niches` for its own
  # (species-level) purposes.
  # NOTE: load_observations() itself already applies .filter_maxillariinae()
  # (drops unidentified rows) before this function ever sees the result --
  # LaElenita's raw rows are ALL unidentified, so load_observations() alone
  # returns zero LaElenita rows, before any lat/lon filtering even happens.
  # Must read the RAW csv directly to actually decouple elevation from
  # identification status.
  elev_source <- read.csv(OBSERVATIONS_CSV)
  elev_source <- elev_source[!is.na(elev_source$lat) & !is.na(elev_source$lon), ]
  elev_source <- augment_elevation(elev_source)

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
    site_elev[site] <- mean(elev_source$Elevation_final_m[elev_source$Area_or_Site == site], na.rm = TRUE)
  }
  list(site_pixels = site_pixels, site_elev = site_elev)
}

# ── Footprint-restricted climate series (2026-09-08, Phase D) ─────────────────
# Repoints Task 1/E2-style analyses at the SAME per-pixel footprint climate
# the demographic model actually uses (voxel_quantiles, via
# build_clim_cache_voxel()/get_clim_voxel()), instead of
# .build_site_climate_series()'s whole-raster spatial pool
# (lookup_climate_by_height()) -- see Item 0's finding, docs/
# methods_update_report.md: these are two genuinely different climate
# systems, and this project's own analysis scripts were reading the wrong
# one for anything claiming to describe "this site's own microclimate."
#
# Per (site, height, footprint pixel): the pixel's own per-height MEDIAN
# (the same deterministic `.sample_quantile(..., stochastic=FALSE)` value
# every other consumer in this pipeline uses under default settings -- not
# a new statistic). temp/relhum use the "both" (day+night pooled) daypart;
# swdown uses "day" ONLY (`swdown>0` hours) -- radiation averaged over the
# full day/night cycle is dominated by the (uninformative) zero-radiation
# night half; every other per-pixel radiation consumer in this codebase
# already reads the day-only bucket (get_clim_voxel()/niche scoring), this
# just makes analysis consistent with that.
#
# Returns, per site: $pixels (site x height x pixel x variable, both
# relative height [height/ceiling] and absolute height carried as columns),
# $footprint_centroid (lon/lat, mean of the site's own footprint pixel
# centers, via .pixel_to_lonlat()), $footprint_elevation (min/mean/max from
# the DTM, sampled at every footprint pixel's own center -- replaces
# site_elev's single observation-coordinate-derived value).
.build_site_climate_series_footprint <- function(sites, niches, processed_dir = PROCESSED_DIR,
                                                  ceiling_fn = site_canopy_ceiling) {
  VARS <- c("temp", "relhum", "swdown")
  daypart_by_var <- c(temp = "both", relhum = "both", swdown = "day")

  out_pixels <- list()
  out_centroid <- list()
  out_elev <- list()
  for (site in sites) {
    microenv_path <- file.path(processed_dir, sprintf("microenv_%s_h0.40.rds", site))
    if (!file.exists(microenv_path)) { message("Skipping ", site, " -- no ", microenv_path); next }
    microenv <- readRDS(microenv_path)
    if (!.spatial_extent_usable(microenv)) { message("Skipping ", site, " -- no usable per-pixel spatial extent"); next }

    obs_site <- niches[niches$Area_or_Site == site, ]
    if (nrow(obs_site) == 0) { message("Skipping ", site, " -- no observations to build a footprint from"); next }
    px_idx <- .lonlat_to_pixel(obs_site$lon, obs_site$lat, microenv$.spatial)
    footprint <- unique(data.frame(row = px_idx$row, col = px_idx$col))

    cc_voxel <- build_clim_cache_voxel(microenv, footprint = footprint)
    heights <- microenv_heights(microenv)
    ceiling <- ceiling_fn(site, niches, max(heights))

    n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))
    rows <- do.call(rbind, parallel::mclapply(seq_along(heights), function(zi) {
      h <- heights[zi]
      per_var <- lapply(VARS, function(v) {
        .clim_voxel_slice(list(clim_cache_voxel = cc_voxel, clim_mode = "voxel",
                               clim_pixel_row = footprint$row, clim_pixel_col = footprint$col,
                               xDim = nrow(footprint), yDim = 1),
                          zi, "annual", daypart_by_var[[v]], v, stochastic = FALSE)
      })
      names(per_var) <- VARS
      n_px <- nrow(footprint)
      data.frame(height = h, rel_height = h / ceiling, pixel = seq_len(n_px),
                temp = as.vector(per_var$temp)[seq_len(n_px)],
                relhum = as.vector(per_var$relhum)[seq_len(n_px)],
                swdown = as.vector(per_var$swdown)[seq_len(n_px)])
    }, mc.cores = n_cores))
    out_pixels[[site]] <- rows

    centers <- .pixel_to_lonlat(footprint$row, footprint$col, microenv$.spatial)
    out_centroid[[site]] <- c(lon = mean(centers$lon), lat = mean(centers$lat))
    elevs <- .dtm_elevation(site, centers$lon, centers$lat)
    out_elev[[site]] <- c(min = min(elevs, na.rm = TRUE), mean = mean(elevs, na.rm = TRUE), max = max(elevs, na.rm = TRUE))
  }
  list(pixels = out_pixels, footprint_centroid = out_centroid, footprint_elevation = out_elev)
}

# Great-circle (haversine) distance in km between two (lon,lat) pairs.
.haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# ── Maxillariinae subtribe filter ────────────────────────────────────────────
# 2026-08-29: OBSERVATIONS_CSV has always carried a handful of non-target
# records -- confirmed by genus-tokenizing Identification directly: Scaphy-
# glottis (4), Brassia (3), Dichaea (1), Stelis (1), Epidendrum (1),
# Elleanthus (1) -- none of them Maxillariinae, none of them ever excluded
# by any prior filter (.confirmed_species_sites() only checks whether an ID
# is genuine, not whether the genus belongs in scope). This is the first
# actual taxonomic scope filter in the pipeline.
#
# Two steps, in order, on the raw Identification column (not the post-
# default FinalID -- this must run before that default is ever applied, so
# it belongs at the very top of load_observations(), before anything else
# touches the data):
#   1. Rename synonym genera to their accepted name (e.g. "Ida sp." ->
#      "Sudamerlycaste sp.") -- these ARE Maxillariinae, just under an
#      outdated generic name, so they're corrected in place rather than
#      dropped.
#   2. Drop any row whose (possibly renamed) genus isn't in the accepted
#      Maxillariinae list at all -- these are real other-subtribe/family
#      misidentifications or contaminant records, not a naming quirk.
# 2026-08-30: blank/NA Identification rows are now DROPPED here too (not
# left for a species-level default) -- load_observations() previously
# defaulted every blank/unidentified row to "Maxillaria acutifolia" as a
# stopgap, but that fabricates a specific species identity for an
# individual nobody actually identified, which is exactly the contamination
# mechanism this whole filter exists to stop (see the module-level comment
# above: LaElenita/MindoMirador/Saloya's niche figures looked convincing
# only because they were borrowing OTHER sites' real acutifolia sightings
# via this same default). Per explicit decision: no placeholder name of any
# kind for an unidentified individual -- drop it, don't guess.
#
# Accepted names per the current Maxillariinae circumscription (WCVP,
# checked 2026-08-29) -- genus only, author citations dropped.
#
# 2026-09-05 CORRECTION: "Sudamerlycaste" (below) was accepted here in
# error. POWO treats Sudamerlycaste Archila as an illegitimate, superfluous
# name for Ida A.Ryan & Oakeley -- Ida is the correct accepted genus, not a
# synonym to be renamed away from. This was backwards: the rename map below
# used to send Ida -> Sudamerlycaste (the wrong direction). Fixed to Ida as
# the accepted name, Sudamerlycaste as the synonym renamed to it -- affects
# the "Sudamerlycaste sp." niche-cache key (now "Ida sp."; same underlying
# observations, same counts, name only).
MAXILLARIINAE_GENERA <- c(
  "Anguloa", "Bifrenaria", "Brasiliorchis", "Cryptocentrum", "Cyrtidiorchis",
  "Guanchezia", "Horvatia", "Hylaeorchis", "Ida", "Lycaste", "Maxillaria",
  "Mormolyca", "Neomoorea", "Pityphyllum", "Rudolfiella", "Scuticaria",
  "Teuscheria", "Trigonidium", "Xylobium"
)

# Synonym genus -> accepted genus. Bifrenaria's three synonyms (Adipe,
# Cydoniorchis, Stenocoryne), the two explicitly called out as "included
# in Maxillaria", and Sudamerlycaste -> Ida (see correction above).
.MAXILLARIINAE_RENAME <- c(
  Anthosiphon    = "Maxillaria",
  Chrysocycnis   = "Maxillaria",
  Sudamerlycaste = "Ida",
  Adipe          = "Bifrenaria",
  Cydoniorchis   = "Bifrenaria",
  Stenocoryne    = "Bifrenaria"
)

# Applies both steps to a niches-style data frame's Identification column.
# Called from load_observations() (paths.R) before anything else, so every
# downstream consumer (FinalID default, characterize_niches.R, all of this
# session's site/species helpers) sees already-scoped, already-renamed data
# without needing its own copy of this logic.
.filter_maxillariinae <- function(niches_df) {
  ident <- niches_df$Identification
  has_id <- !is.na(ident) & nzchar(trimws(ident))
  genus <- rep(NA_character_, length(ident))
  genus[has_id] <- vapply(strsplit(trimws(ident[has_id]), "\\s+"), `[`, character(1), 1)

  # Step 1: rename synonym genera in place, keeping the rest of the
  # Identification string (species epithet, "sp.", etc.) unchanged.
  needs_rename <- has_id & genus %in% names(.MAXILLARIINAE_RENAME)
  if (any(needs_rename)) {
    rest <- sub("^\\S+\\s*", "", trimws(ident[needs_rename]))
    new_genus <- .MAXILLARIINAE_RENAME[genus[needs_rename]]
    niches_df$Identification[needs_rename] <- trimws(paste(new_genus, rest))
    genus[needs_rename] <- unname(new_genus)
  }

  # Step 2: drop rows with a genus that's asserted but out of scope, AND
  # (2026-08-30) rows with no Identification at all -- blank/unidentified is
  # no longer defaulted to any species, so there's nothing left for such a
  # row to usefully contribute to species-level niche characterization.
  out_of_scope <- has_id & !(genus %in% MAXILLARIINAE_GENERA)
  unidentified <- !has_id
  if (any(out_of_scope)) {
    message("Excluding ", sum(out_of_scope), " non-Maxillariinae row(s): ",
            paste(sort(unique(ident[out_of_scope])), collapse = ", "))
  }
  if (any(unidentified)) {
    message("Excluding ", sum(unidentified), " unidentified (blank/NA) row(s) -- ",
            "no species-level placeholder is assigned.")
  }
  niches_df[!out_of_scope & !unidentified, ]
}
