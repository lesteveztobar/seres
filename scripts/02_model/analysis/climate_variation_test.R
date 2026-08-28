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
source("scripts/02_model/config/shared_helpers.R")
library(terra)

# Elevation helpers (augment_elevation(), .dtm_elevation()) and the
# per-site/per-height climate series builder (.build_site_climate_series())
# used below now live in shared_helpers.R -- relocated there 2026-08-27 once
# climate_variation_between_sites.R needed them too.

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

# ── Per-site, per-height, per-variable raw climate (every hourly timestep at
# every height, not just the mean) -- this gives the within-site
# Kruskal-Wallis/Dunn test real statistical power, the same way
# plot_temperature_profile()'s existing test uses every pixel rather than
# per-height means. Shared with climate_variation_between_sites.R -- see
# .build_site_climate_series() (shared_helpers.R). ──────────────────────────
sc <- .build_site_climate_series(SITES, niches = niches)
site_pixels <- sc$site_pixels   # site -> data.frame(height, temp, relhum, swdown)
site_elev   <- sc$site_elev

if (length(site_pixels) == 0) stop("No sites had a usable microenv_<site>_h0.40.rds.")

# ── Test 1: within-site, per variable -- Kruskal-Wallis + Dunn pairwise ────
# eps_sq (epsilon-squared, Tomczak & Tomczak 2014) is Kruskal-Wallis's
# standard companion effect size: (H - k + 1) / (n - k), clamped at 0 --
# added 2026-08-27 because at this test's sample sizes (every hourly
# timestep x every height tier -> tens of thousands of rows per site x
# variable) kw_p underflows to a literal 0 almost everywhere (real vertical
# structure + huge n saturates the chi-squared tail past double-precision
# range), so p alone can't distinguish "barely detectable" from "enormous"
# vertical structure -- eps_sq (0-1, roughly R^2-like) can.
kw_dunn_one <- function(df, var) {
  x <- df[[var]]; g <- factor(df$height)
  ok <- is.finite(x)
  x <- x[ok]; g <- droplevels(g[ok])
  if (nlevels(g) < 2 || length(x) < 4) return(NULL)
  kt <- tryCatch(kruskal.test(x, g), error = function(e) NULL)
  if (is.null(kt)) return(NULL)
  # dunn.test() always emits its full pairwise comparison table regardless
  # of `kw=`/`label=` -- for a site with 100+ height tiers that's a
  # 100x100+ matrix (up to ~5000 pairs) on every one of the 21 site x
  # variable calls this script makes. It emits via rlang::inform(), i.e. a
  # message condition (stderr), not cat()/print() -- capture.output() (which
  # only redirects stdout) does NOT catch it; suppressMessages() does.
  # (2026-08-27: confirmed via `deparse(dunn.test::dunn.test)` after
  # capture.output() alone left runtime completely unchanged.)
  dunn_res <- tryCatch(
    suppressMessages(dunn.test::dunn.test(x, g, method = "holm", kw = FALSE, label = FALSE)),
    error = function(e) NULL)
  n_sig <- if (!is.null(dunn_res)) sum(dunn_res$P.adjusted < 0.05) else NA_integer_
  n_pairs <- if (!is.null(dunn_res)) length(dunn_res$P.adjusted) else NA_integer_
  eps_sq <- max(0, (unname(kt$statistic) - nlevels(g) + 1) / (length(x) - nlevels(g)))
  data.frame(kw_chisq = unname(kt$statistic), kw_df = unname(kt$parameter),
             kw_p = kt$p.value, eps_sq = eps_sq, n_height_tiers = nlevels(g),
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
