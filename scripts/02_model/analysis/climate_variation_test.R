# climate_variation_test.R
# Statistical test of whether temp/relhum/swdown vary (a) by height tier
# within a site, and (b) by elevation across sites -- the two questions
# behind the "resolution test" discussion in the thesis: does coarser height
# sampling risk missing real vertical structure, and does that structure
# itself shift along the elevational gradient the sites span? SITES below is
# derived dynamically (like everywhere else in this pipeline), so this scales
# automatically as more sites are added -- currently seven (2026-07-25).
#
# Generalizes plot_temperature_profile()'s existing Kruskal-Wallis + Dunn
# post-hoc (currently Maquipucuna-only, temperature-only -- see
# plot_functions.R) to all three niche-relevant variables (temp, relhum,
# swdown -- see NICHE_VARS, get_colonization.R) and every site, plus a new
# cross-site elevation comparison that a single site's data can't support at
# all (essentially no elevation range within one site).
#
# Two separate analyses, because "height tier" and "elevation" don't mix
# cleanly into one test -- a given height in meters means something
# different in a 15m-canopy site than a 30m-canopy one:
#   1. WITHIN-SITE, per variable: does the variable differ across height
#      tiers? (Kruskal-Wallis omnibus + Dunn pairwise, Holm-adjusted --
#      same method as plot_temperature_profile(), just for all 3 variables
#      and every site instead of 1 and 1.)
#   2. ACROSS-SITE (elevation), per variable: using each site's own
#      mean climate value (averaged across its full height profile) against
#      that site's mean elevation (see elevation_helpers.R) -- the standard
#      "does climate follow a lapse-rate-like trend with elevation" check.
#      With only a handful of sites, this is too small a sample for a formal
#      test, so it's reported as a Spearman correlation (rank-based, robust
#      to the small n and to any non-linearity) rather than claiming
#      inferential power it doesn't have -- read this as descriptive/
#      indicative, not a confirmatory test, regardless of how many sites
#      are currently in the data.
#
# Usage: Rscript scripts/02_model/analysis/climate_variation_test.R
# Output: output/climate_variation_by_height.csv (test 1, one row per
#   site x variable), output/climate_variation_by_elevation.csv (test 2,
#   one row per variable), printed summary of both.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")
library(terra)

# ── Elevation helpers (formerly elevation_helpers.R; merged here 2026-07-28 —
# ── this was its only sourcer) ──────────────────────────────────────────────
# Per-observation elevation, combining field-recorded Elevation_m (present
# for Maquipucuna/Mashpi/MiradorMindo) with digital-elevation-model
# extraction at each observation's exact lat/lon (fills MindoTarabita/
# Yanayacu, which have NO recorded elevation at all -- audit 2026-07-18:
# 29/57 combinedv3.csv rows had Elevation_m == NA, entirely concentrated in
# those two sites; the other 3 sites are 100% populated already).
#
# Reuses the per-site DTM already downloaded for microclimate modelling
# (data/raw/<site>/dtm/dtm.tif, cached by get_dtm() in lib.R,
# already reprojected to EPSG:4326 there) -- no new data acquisition needed,
# every site already has one.
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
# ── End elevation helpers ───────────────────────────────────────────────────

VARS <- NICHE_VARS  # c("temp", "relhum", "swdown") -- get_colonization.R

if (!requireNamespace("dunn.test", quietly = TRUE)) {
  stop("dunn.test package not installed. Install once with:\n",
       '  Rscript -e \'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("dunn.test")\'')
}

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
niches <- augment_elevation(niches)

# Sites derived from OBSERVATIONS_CSV itself, not hardcoded -- a newly added
# site is picked up automatically the next time this runs (its own
# per-site loop below already skips gracefully if that site's microenv
# doesn't exist yet).
SITES <- sort(unique(niches$Area_or_Site))

# ── Per-site, per-height, per-variable raw climate (every spatial pixel at
# every height, not just the mean) -- pixel-level values give the
# within-site Kruskal-Wallis/Dunn test real statistical power, the same way
# plot_temperature_profile()'s existing test uses every pixel rather than
# per-height means. ──────────────────────────────────────────────────────────
site_pixels   <- list()   # site -> data.frame(height, temp, relhum, swdown)
site_elev     <- numeric(0)

for (site in SITES) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.25.rds", site))
  if (!file.exists(microenv_path)) {
    message("Skipping ", site, " -- no ", microenv_path)
    next
  }
  message("Reading ", site, "...")
  microenv <- readRDS(microenv_path)
  heights  <- microenv_heights(microenv)

  # swdown: keep every row, including night/zero -- unlike the niche axis
  # (which deliberately drops night hours, see .niche_var_mean() in
  # get_colonization.R), here we want the full diurnal distribution actually
  # sampled at each height, not a daytime-only summary.
  # Parallel, not sequential: each lookup_climate_by_height() call reads one
  # independent per-height file from scratch (see load_height()), same as
  # build_clim_cache() (get_colonization.R) -- this loop hit its own SLURM
  # time limit on 2026-07-27 running single-threaded on sites with 100-400+
  # height tiers.
  n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))
  px <- do.call(rbind, parallel::mclapply(heights, function(h) {
    cl <- lookup_climate_by_height(h, microenv)
    if (is.null(cl)) return(NULL)
    data.frame(height = h, temp = cl$temp, relhum = cl$relhum, swdown = cl$swdown)
  }, mc.cores = n_cores))
  site_pixels[[site]] <- px
  site_elev[site] <- mean(niches$Elevation_final_m[niches$Area_or_Site == site], na.rm = TRUE)
}

if (length(site_pixels) == 0) stop("No sites had a usable microenv_<site>_h0.25.rds.")

# ── Test 1: within-site, per variable -- Kruskal-Wallis + Dunn pairwise ────
kw_dunn_one <- function(df, var) {
  x <- df[[var]]; g <- factor(df$height)
  ok <- is.finite(x)
  x <- x[ok]; g <- droplevels(g[ok])
  if (nlevels(g) < 2 || length(x) < 4) return(NULL)
  kt <- tryCatch(kruskal.test(x, g), error = function(e) NULL)
  if (is.null(kt)) return(NULL)
  dunn_res <- tryCatch(dunn.test::dunn.test(x, g, method = "holm", kw = FALSE, label = FALSE),
                       error = function(e) NULL)
  n_sig <- if (!is.null(dunn_res)) sum(dunn_res$P.adjusted < 0.05) else NA_integer_
  n_pairs <- if (!is.null(dunn_res)) length(dunn_res$P.adjusted) else NA_integer_
  data.frame(kw_chisq = unname(kt$statistic), kw_df = unname(kt$parameter),
             kw_p = kt$p.value, n_height_tiers = nlevels(g),
             n_sig_pairs = n_sig, n_pairs = n_pairs)
}

within_rows <- list()
for (site in names(site_pixels)) {
  for (v in VARS) {
    res <- kw_dunn_one(site_pixels[[site]], v)
    if (!is.null(res)) within_rows[[length(within_rows) + 1]] <- cbind(site = site, variable = v, res)
  }
}
within_df <- do.call(rbind, within_rows)

cat("\n========================================\n")
cat("1. Within-site variation by height tier (Kruskal-Wallis + Dunn pairwise)\n")
cat("========================================\n")
print(within_df, row.names = FALSE)
out1 <- file.path(OUTPUT_DIR, "climate_variation_by_height.csv")
write.csv(within_df, out1, row.names = FALSE)
cat("\nSaved: ", out1, "\n")

# ── Test 2: across-site (elevation) -- site-mean climate vs. site elevation ──
site_means <- do.call(rbind, lapply(names(site_pixels), function(site) {
  px <- site_pixels[[site]]
  row <- c(site = site, elevation_m = unname(site_elev[site]),
          setNames(sapply(VARS, function(v) mean(px[[v]], na.rm = TRUE)), VARS))
  as.data.frame(as.list(row), stringsAsFactors = FALSE)
}))
for (v in c("elevation_m", VARS)) site_means[[v]] <- as.numeric(site_means[[v]])

cat("\n========================================\n")
cat("2. Site-mean climate vs. elevation (n=", nrow(site_means), " sites)\n", sep = "")
cat("========================================\n")
print(site_means, row.names = FALSE)

elev_rows <- lapply(VARS, function(v) {
  ok <- is.finite(site_means$elevation_m) & is.finite(site_means[[v]])
  if (sum(ok) < 3) return(data.frame(variable = v, n = sum(ok), rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(site_means$elevation_m[ok], site_means[[v]][ok], method = "spearman"))
  data.frame(variable = v, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
})
elev_df <- do.call(rbind, elev_rows)
cat("\nSpearman correlation with elevation (indicative only, n=", nrow(site_means), " sites):\n", sep = "")
print(elev_df, row.names = FALSE)
out2 <- file.path(OUTPUT_DIR, "climate_variation_by_elevation.csv")
write.csv(elev_df, out2, row.names = FALSE)
cat("\nSaved: ", out2, "\n")
