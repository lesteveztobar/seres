# Task 1 / "E2" -- vertical-profile slope (climate vs height) per variable,
# then whether that slope covaries with elevation across sites. REWRITTEN
# 2026-09-08 (Phase D) to read the same per-pixel FOOTPRINT climate the
# demographic model actually uses (.build_site_climate_series_footprint(),
# shared_helpers.R -- voxel_quantiles via build_clim_cache_voxel()), not
# the whole-raster-pooled lookup_climate_by_height() the old version used
# (see docs/methods_update_report.md, Item 0 -- these are two genuinely
# different climate systems). Also: swdown is daylight-hours only (built
# into the new function); height reported BOTH as absolute 5m bins and
# relative-height deciles, side by side; rel_height<0.3 split dropped
# entirely (was never well-justified -- see the report's own flag);
# elevation is now each site's own FOOTPRINT elevation (DTM-sampled at
# every footprint pixel, not observation coordinates); the elevation~slope
# correlation is reported at LOCALITY level (Mindo cluster -- MindoMirador/
# MindoTarabita/Saloya/LaElenita -- collapsed to one point, n=4) as
# PRIMARY, and at site level (n=6/7) as a sensitivity check, per the
# confirmed finding that those 4 sites share a single ERA5 forcing cell
# and substantially overlapping terrain -- treating them as 4 independent
# samples was never defensible.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
library(MASS)  # rlm() -- robust (IRLS) regression, used when residuals misbehave

VARS <- c("temp", "relhum", "swdown")
SITES <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu", "LaElenita")
# CANOPY_SKIP_SITES (comma/space-separated, default empty): drop sites,
# e.g. CANOPY_SKIP_SITES=Saloya -- same convention as sensitivity_run.R.
.skip <- strsplit(Sys.getenv("CANOPY_SKIP_SITES", unset = ""), "[ ,]+")[[1]]
if (length(.skip) && any(nzchar(.skip))) {
  message(sprintf("CANOPY_SKIP_SITES=%s: dropping from SITES", paste(.skip, collapse = ",")))
  SITES <- setdiff(SITES, .skip)
}
MINDO_CLUSTER <- c("MindoMirador", "MindoTarabita", "Saloya", "LaElenita")

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
niches <- augment_elevation(niches)

fc <- .build_site_climate_series_footprint(SITES, niches)
site_pixels <- fc$pixels
cat("Sites entering the analysis:", length(site_pixels), "of", length(SITES), "requested\n")
cat("Sites used:", paste(names(site_pixels), collapse = ", "), "\n\n")

ceiling_for <- function(site) site_canopy_ceiling(site, niches, max(site_pixels[[site]]$height))
site_ceiling <- setNames(sapply(names(site_pixels), ceiling_for), names(site_pixels))

# ── Great-circle distances between footprint centroids ─────────────────────
cent <- fc$footprint_centroid
site_names <- names(cent)
dist_mat <- matrix(NA_real_, length(site_names), length(site_names), dimnames = list(site_names, site_names))
for (i in site_names) for (j in site_names) {
  dist_mat[i, j] <- .haversine_km(cent[[i]]["lon"], cent[[i]]["lat"], cent[[j]]["lon"], cent[[j]]["lat"])
}
cat("=== Pairwise great-circle distances between footprint centroids (km) ===\n")
print(round(dist_mat, 2))
write.csv(dist_mat, resolve_output_path("climate_variation_relative_height", "footprint_centroid_distances_km.csv",
  legacy_path = file.path(OUTPUT_DIR, "footprint_centroid_distances_km.csv")))

# ── Per-site footprint elevation (replaces site_elev) ───────────────────────
elev_df <- do.call(rbind, lapply(names(fc$footprint_elevation), function(s) {
  e <- fc$footprint_elevation[[s]]
  data.frame(site = s, elev_min = e["min"], elev_mean = e["mean"], elev_max = e["max"])
}))
rownames(elev_df) <- NULL
cat("\n=== Per-site footprint elevation (DTM-sampled at every footprint pixel) ===\n")
print(elev_df, row.names = FALSE)
write.csv(elev_df, resolve_output_path("climate_variation_relative_height", "footprint_elevation.csv",
  legacy_path = file.path(OUTPUT_DIR, "footprint_elevation.csv")), row.names = FALSE)
site_elev_mean <- setNames(elev_df$elev_mean, elev_df$site)

# ── Part A: within-site slope, BOTH absolute-height bins and relative-height deciles ──
fit_slope <- function(yy, xx) {
  ols <- lm(yy ~ xx)
  resid_ok <- tryCatch({
    rs <- resid(ols)
    samp <- if (length(rs) > 4000) sample(rs, 4000) else rs
    shapiro.test(samp)$p.value
  }, error = function(e) NA_real_)
  use_robust <- is.na(resid_ok) || resid_ok < 0.05
  if (use_robust) {
    rfit <- tryCatch(MASS::rlm(yy ~ xx, maxit = 100), error = function(e) NULL)
    if (!is.null(rfit)) {
      co <- summary(rfit)$coefficients
      return(list(slope = co["xx", "Value"], se = co["xx", "Std. Error"],
                  method = "rlm", resid_p = resid_ok))
    }
  }
  co <- summary(ols)$coefficients
  ci <- confint(ols)["xx", ]
  list(slope = co["xx", "Estimate"], se = co["xx", "Std. Error"],
       method = if (use_robust) "OLS (rlm failed, fell back)" else "OLS", resid_p = resid_ok)
}

rows_abs <- list(); rows_rel <- list()
for (site in names(site_pixels)) {
  px <- site_pixels[[site]]
  for (v in VARS) {
    y <- px[[v]]
    ok_abs <- is.finite(y) & is.finite(px$height)
    ok_rel <- is.finite(y) & is.finite(px$rel_height)
    if (sum(ok_abs) >= 10) {
      f <- fit_slope(y[ok_abs], px$height[ok_abs])
      rows_abs[[length(rows_abs) + 1]] <- data.frame(
        site = site, variable = v, n = sum(ok_abs), binning = "absolute_height_m",
        slope = f$slope, se = f$se, method = f$method,
        spearman_rho = suppressWarnings(cor(px$height[ok_abs], y[ok_abs], method = "spearman")))
    }
    if (sum(ok_rel) >= 10) {
      f <- fit_slope(y[ok_rel], px$rel_height[ok_rel])
      rows_rel[[length(rows_rel) + 1]] <- data.frame(
        site = site, variable = v, n = sum(ok_rel), binning = "relative_height_decile",
        slope = f$slope, se = f$se, method = f$method,
        spearman_rho = suppressWarnings(cor(px$rel_height[ok_rel], y[ok_rel], method = "spearman")))
    }
  }
}
slope_df <- rbind(do.call(rbind, rows_abs), do.call(rbind, rows_rel))
cat("\n=== Part A: per-site, per-variable slope, BOTH bandings side by side ===\n")
print(slope_df, row.names = FALSE)
write.csv(slope_df, resolve_output_path("climate_variation_relative_height", "elevation_test_relative_height_slopes.csv",
  legacy_path = file.path(OUTPUT_DIR, "elevation_test_relative_height_slopes.csv")), row.names = FALSE)

# ── Part B: does the per-site slope vary with elevation? Locality primary, site sensitivity ──
elev_test <- function(sites_subset, elev_lookup, label, slope_binning) {
  sub_slopes <- slope_df[slope_df$binning == slope_binning & slope_df$site %in% sites_subset, ]
  rows <- lapply(VARS, function(v) {
    sub <- sub_slopes[sub_slopes$variable == v, ]
    sub$elev <- elev_lookup[sub$site]
    ok <- is.finite(sub$elev) & is.finite(sub$slope)
    if (sum(ok) < 3) return(data.frame(set = label, variable = v, n = sum(ok), rho = NA, p = NA))
    ct <- suppressWarnings(cor.test(sub$elev[ok], sub$slope[ok], method = "spearman"))
    data.frame(set = label, variable = v, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
  })
  do.call(rbind, rows)
}

# Locality collapse: Mindo cluster (4 sites sharing one ERA5 cell and
# substantially overlapping terrain -- confirmed this session) averaged to
# ONE locality point; Maquipucuna, Mashpi, Yanayacu each their own
# locality (no other site shares their ERA5 cell). n=4 localities (or 3 if
# LaElenita/the cluster isn't fully present).
localities <- list(
  Mindo_cluster = intersect(MINDO_CLUSTER, names(site_pixels)),
  Maquipucuna = "Maquipucuna", Mashpi = "Mashpi", Yanayacu = "Yanayacu")
locality_slope <- do.call(rbind, lapply(names(localities), function(loc) {
  sites_in <- intersect(localities[[loc]], names(site_pixels))
  if (length(sites_in) == 0) return(NULL)
  sub <- slope_df[slope_df$binning == "relative_height_decile" & slope_df$site %in% sites_in, ]
  do.call(rbind, lapply(VARS, function(v) {
    s <- sub$slope[sub$variable == v]
    if (length(s) == 0) return(NULL)
    data.frame(locality = loc, variable = v, slope = mean(s), n_sites_averaged = length(s))
  }))
}))
locality_elev <- setNames(sapply(names(localities), function(loc) {
  sites_in <- intersect(localities[[loc]], names(site_elev_mean))
  if (length(sites_in) == 0) return(NA_real_)
  mean(site_elev_mean[sites_in])
}), names(localities))

cat("\n=== Part B PRIMARY: locality-level (Mindo cluster collapsed), n=", length(localities), " ===\n", sep = "")
loc_rows <- lapply(VARS, function(v) {
  sub <- locality_slope[locality_slope$variable == v, ]
  sub$elev <- locality_elev[sub$locality]
  ok <- is.finite(sub$elev) & is.finite(sub$slope)
  if (sum(ok) < 3) return(data.frame(variable = v, n = sum(ok), rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(sub$elev[ok], sub$slope[ok], method = "spearman"))
  data.frame(variable = v, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
})
locality_elev_df <- do.call(rbind, loc_rows)
print(locality_elev_df, row.names = FALSE)
write.csv(locality_slope, resolve_output_path("climate_variation_relative_height", "elevation_test_locality_slopes.csv",
  legacy_path = file.path(OUTPUT_DIR, "elevation_test_locality_slopes.csv")), row.names = FALSE)
write.csv(locality_elev_df, resolve_output_path("climate_variation_relative_height", "elevation_test_locality_vs_elevation.csv",
  legacy_path = file.path(OUTPUT_DIR, "elevation_test_locality_vs_elevation.csv")), row.names = FALSE)

cat("\n=== Part B SENSITIVITY: site-level (n=", length(site_pixels), "), for comparison only ===\n", sep = "")
site_elev_df <- elev_test(names(site_pixels), site_elev_mean, "site_level", "relative_height_decile")
print(site_elev_df, row.names = FALSE)
write.csv(site_elev_df, resolve_output_path("climate_variation_relative_height", "elevation_test_slope_vs_elevation.csv",
  legacy_path = file.path(OUTPUT_DIR, "elevation_test_slope_vs_elevation.csv")), row.names = FALSE)

cat("\n(All correlations descriptive, not confirmatory -- n<=7 sites / n=4 localities.)\n")
cat("Done.\n")
