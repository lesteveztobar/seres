# climate_variation_between_sites.R
# Two BETWEEN-SITE climate comparisons, complementing climate_variation_
# test.R's two existing analyses (within-site by height tier, and the
# 7-site Spearman correlation of site-mean climate vs. site-mean elevation).
# Both of those either stay inside one site or collapse each site to a
# single mean -- neither one asks "do two different sites actually feel
# different, at a comparable canopy position or a comparable absolute
# elevation?" This script answers that, two ways:
#
#   A. SAME RELATIVE HEIGHT TIER, across all sites -- for a handful of
#      reference heights-above-ground (1/2/5/10/20 m), find each site's
#      nearest available height tier and compare (Kruskal-Wallis + Dunn
#      pairwise, same method as Test 1) across ALL sites present at that
#      reference height, then call out the highest-vs-lowest-elevation
#      pair by name so "does the high site feel different from the low
#      site at the same canopy position" has a direct answer.
#   B. SAME ABSOLUTE ELEVATION (masl), across sites -- height tiers aren't
#      pinned to a location, so "meters above ground" means a different
#      absolute elevation at every site (site's own mean elevation +
#      height). Binning by that combined absolute elevation and comparing
#      sites within each bin asks the sharper question: if you stand at
#      the SAME point in space (masl-wise), does which site you're at
#      still matter? Only sites whose absolute-elevation range actually
#      overlaps another site's can appear together in a bin -- see the
#      printed overlap summary; some sites (this dataset's elevational
#      extremes) may have no valid partner at all, which is a real
#      property of the site selection, not a bug in the test.
#
# Both reuse .build_site_climate_series() (shared_helpers.R), the same
# per-height hourly climate series climate_variation_test.R's Test 1 uses --
# not literally per-pixel (see that function's own comment): each height
# tier's raster is already spatially averaged to one value per hourly
# timestep at write time, so cross-site absolute-elevation matching here is
# "site's mean recorded/DTM elevation + height tier", not a true per-pixel
# DTM lookup -- the per-pixel spatial detail needed for that was already
# discarded upstream (run_microclimate_site.R's write-time reduction), so
# this is the finest resolution available without re-deriving the raw
# microclimf rasters.
#
# Usage: Rscript scripts/experiments/A08_climate_variation_between_sites/climate_variation_between_sites.R
# Output: output/climate_variation_between_sites_by_height.csv (analysis A),
#   output/climate_variation_between_sites_by_elevation.csv (analysis B),
#   printed summary of both plus the elevation-overlap table.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/config/shared_helpers.R")
library(terra)

if (!requireNamespace("dunn.test", quietly = TRUE)) {
  stop("dunn.test package not installed. Install once with:\n",
       '  Rscript -e \'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("dunn.test")\'')
}

VARS <- NICHE_VARS  # c("temp", "relhum", "swdown") -- get_colonization.R
REF_HEIGHTS <- c(1, 2, 5, 10, 20)  # m above ground -- analysis A's reference tiers
HEIGHT_TOL  <- 0.5                 # m -- max distance from a ref height to count as "present"
ELEV_BIN_M  <- 25                  # m -- analysis B's absolute-elevation bin width

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
niches <- augment_elevation(niches)
SITES  <- sort(unique(niches$Area_or_Site))

sc <- .build_site_climate_series(SITES, niches = niches)
site_pixels <- sc$site_pixels
site_elev   <- sc$site_elev
if (length(site_pixels) < 2) stop("Need at least 2 sites with a usable microenv_<site>_h0.40.rds.")

kw_dunn_between <- function(long_df, extreme_lo, extreme_hi) {
  # long_df: data.frame(site, value). Returns one summary row (or NULL if
  # fewer than 2 sites are present).
  #
  # `extreme_pair`/`extreme_p_adj`/`extreme_direction` describe the pair of
  # sites ACTUALLY PRESENT in this bin/height with the largest median
  # difference -- NOT the two sites with the highest/lowest elevation
  # overall (`extreme_lo`/`extreme_hi`, still used below for the separate
  # elev_extreme_* columns). 2026-09-28 fix: the previous version wrote
  # `extreme_pair = paste0(extreme_hi, " vs ", extreme_lo)` unconditionally,
  # using the two GLOBAL elevational-extreme site names regardless of
  # whether either was actually a member of `g` for this specific
  # bin/height -- e.g. it printed "Yanayacu vs Mashpi" for elevation bins
  # containing neither site, because those two happened to be the dataset's
  # global elevation extremes. Per spec: NA when only 2 sites are present in
  # this bin/height (nothing to single out -- the KW/Dunn row already IS
  # that one pair), otherwise the in-bin pair with the largest |median
  # difference|.
  g <- factor(long_df$site)
  x <- long_df$value
  ok <- is.finite(x)
  x <- x[ok]; g <- droplevels(g[ok])
  if (nlevels(g) < 2 || length(x) < 4) return(NULL)
  kt <- tryCatch(kruskal.test(x, g), error = function(e) NULL)
  if (is.null(kt)) return(NULL)
  # See climate_variation_test.R's kw_dunn_one() -- dunn.test() always
  # emits its full pairwise table via rlang::inform() (a message condition,
  # not cat()/print()), regardless of `kw=`/`label=`; capture.output() does
  # NOT catch it (only redirects stdout) -- suppressMessages() does, and
  # that emission dominates runtime once a group has 20+ levels.
  dunn_res <- tryCatch(
    suppressMessages(dunn.test::dunn.test(x, g, method = "holm", kw = FALSE, label = TRUE)),
    error = function(e) NULL)
  eps_sq <- max(0, (unname(kt$statistic) - nlevels(g) + 1) / (length(x) - nlevels(g)))

  # ---- elevational-extreme pair (only when BOTH are members of this bin) --
  elev_extreme_p <- NA_real_; elev_extreme_dir <- NA_character_
  if (!is.null(dunn_res) && extreme_lo %in% levels(g) && extreme_hi %in% levels(g)) {
    lbl1 <- paste(extreme_hi, "-", extreme_lo)
    lbl2 <- paste(extreme_lo, "-", extreme_hi)
    hit <- which(dunn_res$comparisons %in% c(lbl1, lbl2))
    if (length(hit) == 1) {
      elev_extreme_p <- dunn_res$P.adjusted[hit]
      med_hi <- median(x[g == extreme_hi]); med_lo <- median(x[g == extreme_lo])
      elev_extreme_dir <- if (med_hi > med_lo) paste0(extreme_hi, ">", extreme_lo)
                          else if (med_hi < med_lo) paste0(extreme_lo, ">", extreme_hi)
                          else "equal"
    }
  }

  # ---- in-bin extreme pair (largest median difference among sites actually
  # present here) -- NA when there are only 2 sites (nothing to single out).
  extreme_pair <- NA_character_; extreme_p <- NA_real_; extreme_dir <- NA_character_
  if (nlevels(g) > 2) {
    meds <- tapply(x, g, median)
    site_pairs <- combn(levels(g), 2, simplify = FALSE)
    diffs <- vapply(site_pairs, function(p) abs(meds[[p[1]]] - meds[[p[2]]]), numeric(1))
    best <- site_pairs[[which.max(diffs)]]
    med1 <- meds[[best[1]]]; med2 <- meds[[best[2]]]
    hi <- if (med1 >= med2) best[1] else best[2]
    lo <- if (med1 >= med2) best[2] else best[1]
    extreme_pair <- paste0(hi, " vs ", lo)
    extreme_dir  <- paste0(hi, ">", lo)
    if (!is.null(dunn_res)) {
      lbl1 <- paste(hi, "-", lo); lbl2 <- paste(lo, "-", hi)
      hit <- which(dunn_res$comparisons %in% c(lbl1, lbl2))
      if (length(hit) == 1) extreme_p <- dunn_res$P.adjusted[hit]
    }
  }

  data.frame(kw_chisq = unname(kt$statistic), kw_df = unname(kt$parameter),
             kw_p = kt$p.value, eps_sq = eps_sq, n_sites = nlevels(g), n = length(x),
             sites = paste(levels(g), collapse = ","),
             extreme_pair = extreme_pair,
             extreme_p_adj = extreme_p, extreme_direction = extreme_dir,
             elev_extreme_p_adj = elev_extreme_p, elev_extreme_direction = elev_extreme_dir)
}

# Elevational extremes, by each site's own mean elevation (site_elev) --
# named explicitly in every printed/saved row so "highest vs lowest" always
# refers to real site names, not just an index.
site_extreme_hi <- names(which.max(site_elev[names(site_pixels)]))
site_extreme_lo <- names(which.min(site_elev[names(site_pixels)]))
cat(sprintf("\nElevational extremes: highest = %s (%.0f m), lowest = %s (%.0f m)\n",
            site_extreme_hi, site_elev[site_extreme_hi], site_extreme_lo, site_elev[site_extreme_lo]))

# ── Analysis A: same relative height tier, across sites ─────────────────────
cat("\n========================================\n")
cat("A. Between-site comparison at matched height-above-ground tiers\n")
cat("========================================\n")

rows_a <- list()
for (ref_h in REF_HEIGHTS) {
  for (v in VARS) {
    long_rows <- lapply(names(site_pixels), function(site) {
      px <- site_pixels[[site]]
      idx <- which.min(abs(px$height - ref_h))
      if (length(idx) == 0 || abs(px$height[idx] - ref_h) > HEIGHT_TOL) return(NULL)
      nearest_h <- px$height[idx]
      data.frame(site = site, value = px[[v]][px$height == nearest_h])
    })
    long_df <- do.call(rbind, long_rows)
    if (is.null(long_df)) next
    res <- kw_dunn_between(long_df, site_extreme_lo, site_extreme_hi)
    if (!is.null(res)) {
      rows_a[[length(rows_a) + 1]] <- cbind(ref_height_m = ref_h, variable = v, res)
    }
  }
}
between_height_df <- do.call(rbind, rows_a)
if (is.null(between_height_df)) {
  message("No reference height had 2+ sites within tolerance -- nothing to report for analysis A.")
} else {
  print(between_height_df, row.names = FALSE)
}
out_a <- file.path(OUTPUT_DIR, "climate_variation_between_sites_by_height.csv")
write.csv(between_height_df, out_a, row.names = FALSE)
cat("\nSaved: ", out_a, "\n")

# ── Analysis B: same absolute elevation (masl), across sites ────────────────
# First, report which sites CAN be masl-matched at all -- a site's usable
# absolute-elevation range is [site_elev, site_elev + max(height)]; two
# sites can only ever land in the same bin if these ranges overlap.
cat("\n========================================\n")
cat("Site absolute-elevation ranges (site elevation + height-tier range)\n")
cat("========================================\n")
elev_ranges <- do.call(rbind, lapply(names(site_pixels), function(site) {
  h <- site_pixels[[site]]$height
  data.frame(site = site, elev_lo = site_elev[site], elev_hi = site_elev[site] + max(h))
}))
print(elev_ranges, row.names = FALSE)

# Overlap must be checked in the SAME bin-floor terms the comparison below
# actually uses (floor(x / ELEV_BIN_M) * ELEV_BIN_M), not raw continuous
# ranges -- a continuous-range check can say "no overlap" for two ranges
# with, say, a 15 m real gap between them, while ELEV_BIN_M=25 m bins
# still land both sites in the same bin (their nearest boundary points
# floor to the same multiple of 25). Checking continuous ranges here would
# silently contradict analysis B's own bin table below.
elev_ranges$bin_lo <- floor(elev_ranges$elev_lo / ELEV_BIN_M) * ELEV_BIN_M
elev_ranges$bin_hi <- floor(elev_ranges$elev_hi / ELEV_BIN_M) * ELEV_BIN_M
overlap_any <- FALSE
for (i in seq_len(nrow(elev_ranges))) {
  for (j in seq_len(nrow(elev_ranges))) {
    if (j <= i) next
    lo <- max(elev_ranges$bin_lo[i], elev_ranges$bin_lo[j])
    hi <- min(elev_ranges$bin_hi[i], elev_ranges$bin_hi[j])
    if (lo <= hi) {
      overlap_any <- TRUE
      cat(sprintf("  %s <-> %s share elevation bin(s) [%.0f, %.0f] m\n",
                  elev_ranges$site[i], elev_ranges$site[j], lo, hi + ELEV_BIN_M))
    }
  }
}
if (!overlap_any) {
  message("No two sites share an elevation bin at this bin width (",
          ELEV_BIN_M, " m) -- analysis B has nothing to compare (every ",
          "site occupies a disjoint elevational band; this is a property ",
          "of the site selection, not a bug). A wider ELEV_BIN_M would ",
          "bring elevationally-adjacent sites together at the cost of a ",
          "coarser match.")
}

cat("\n========================================\n")
cat(sprintf("B. Between-site comparison at matched absolute elevation (%d m bins)\n", ELEV_BIN_M))
cat("========================================\n")

all_long <- do.call(rbind, lapply(names(site_pixels), function(site) {
  px <- site_pixels[[site]]
  data.frame(site = site, absolute_elev = site_elev[site] + px$height,
             temp = px$temp, relhum = px$relhum, swdown = px$swdown)
}))
all_long$elev_bin <- floor(all_long$absolute_elev / ELEV_BIN_M) * ELEV_BIN_M

rows_b <- list()
for (bin in sort(unique(all_long$elev_bin))) {
  bin_df <- all_long[all_long$elev_bin == bin, ]
  if (length(unique(bin_df$site)) < 2) next
  for (v in VARS) {
    long_df <- data.frame(site = bin_df$site, value = bin_df[[v]])
    res <- kw_dunn_between(long_df, site_extreme_lo, site_extreme_hi)
    if (!is.null(res)) {
      rows_b[[length(rows_b) + 1]] <- cbind(
        elev_bin_lo = bin, elev_bin_hi = bin + ELEV_BIN_M, variable = v, res)
    }
  }
}
between_elev_df <- do.call(rbind, rows_b)
if (is.null(between_elev_df)) {
  message("No elevation bin had 2+ sites -- nothing to report for analysis B ",
          "(see the overlap table above).")
} else {
  print(between_elev_df, row.names = FALSE)
}
out_b <- file.path(OUTPUT_DIR, "climate_variation_between_sites_by_elevation.csv")
write.csv(between_elev_df, out_b, row.names = FALSE)
cat("\nSaved: ", out_b, "\n")
