# step2_elevation_bins_v2.R -- 2026-09-28. Elevation-matched between-site
# comparison, redone per the Step-0 audit + user's confirmed methodology:
#   - tiers aggregated to ONE VALUE PER SITE PER TIMESTAMP (median across
#     the height tiers of that site landing in the 25 m absolute-elevation
#     bin) -- fixes the old n=114,192/263,520 tier pseudoreplication.
#   - PRIMARY comparison is CLIMATOLOGICAL: each site's per-timestamp series
#     is collapsed to (month x hour-of-day) cell medians, sites are paired
#     by (month,hour) cell -- exact calendar-date pairing is impossible for
#     12/21 site pairs (Step 0 finding: only 9/21 pairs share ANY exact
#     timestamp). A ROBUSTNESS check re-runs the same test on exact-date-
#     matched raw timestamps for bin/variable combinations where every site
#     present belongs to a mutually-overlapping pair.
#   - swdown: day-only (>10 W/m2), already clipped to >=0 upstream.
#   - VPD computed per-hour in step1 (not from summarised T/RH).
suppressPackageStartupMessages({})
source("scripts/02_model/config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- if (length(args) >= 1) args[1] else stop("Usage: step2_elevation_bins_v2.R <outdir from step1>")

cache <- readRDS(file.path(PROCESSED_DIR, "site_climate_cache_v2.rds"))
site_pixels <- readRDS(file.path(OUTDIR, "site_pixels_with_vpd.rds"))  # has vpd_kPa
site_elev   <- cache$site_elev
tme_by_site <- cache$tme_by_site
SITES <- names(site_pixels)
ELEV_BIN_M <- 25
VARS <- c("temp", "relhum", "vpd_kPa", "swdown")
DAY_THRESH <- 10
set.seed(20260928)

# ---- Per-site: per-timestamp series (one row per hour), aggregating tiers
# within each 25 m absolute-elevation bin by MEDIAN, per variable. swdown
# uses ONLY daytime hours (its own within-bin tier median computed only
# over daytime tiers at that hour) -- temp/relhum/vpd use all hours.
build_site_bin_series <- function(site) {
  px <- site_pixels[[site]]
  tme <- tme_by_site[[site]]
  n_hr <- length(tme)
  px$hour_idx <- ave(seq_len(nrow(px)), px$height, FUN = seq_along)
  px$abs_elev <- site_elev[site] + px$height
  px$elev_bin <- floor(px$abs_elev / ELEV_BIN_M) * ELEV_BIN_M
  px$is_day <- px$swdown > DAY_THRESH
  list(px = px, tme = tme, n_hr = n_hr)
}
site_bin <- lapply(SITES, build_site_bin_series); names(site_bin) <- SITES

all_bins <- sort(unique(unlist(lapply(SITES, function(s) unique(site_bin[[s]]$px$elev_bin)))))

# For a given (site, bin, variable): one value per hour (1..n_hr), median
# across in-bin tiers; NA where no in-bin tier is daytime (swdown) or no
# in-bin tier exists at all for that site.
site_bin_hourly <- function(site, bin, var) {
  d <- site_bin[[site]]$px
  sub <- d[d$elev_bin == bin, ]
  if (nrow(sub) == 0) return(NULL)
  n_hr <- site_bin[[site]]$n_hr
  if (var == "swdown") sub <- sub[sub$is_day, ]
  if (nrow(sub) == 0) return(rep(NA_real_, n_hr))
  agg <- tapply(sub[[var]], sub$hour_idx, median, na.rm = TRUE)
  out <- rep(NA_real_, n_hr)
  out[as.integer(names(agg))] <- as.numeric(agg)
  out
}

# ---- climatology: (month, hour-of-day) cell median, from a per-timestamp
# vector (length n_hr, aligned to that site's own tme).
climatology <- function(site, x) {
  tme <- tme_by_site[[site]]
  mon <- as.integer(format(tme, "%m")); hr <- as.integer(format(tme, "%H"))
  cell <- sprintf("%02d_%02d", mon, hr)
  ok <- is.finite(x)
  if (!any(ok)) return(NULL)
  m <- tapply(x[ok], cell[ok], median, na.rm = TRUE)
  data.frame(cell = names(m), value = as.numeric(m))
}

# ---- day-block bootstrap CI on a paired-difference vector, resampling
# SAMPLED DAYS (each site's own 24h blocks; every-other-day design -- see
# Step 0.1) with a given block length (in sampled days), 2000 reps.
block_bootstrap_ci <- function(diffs, day_id, block_len = 1, reps = 2000, seed = 20260928) {
  set.seed(seed)
  days <- unique(day_id)
  n_days <- length(days)
  if (n_days < block_len * 2) return(c(lo = NA_real_, hi = NA_real_, eff_n = NA_real_))
  starts <- seq_len(n_days - block_len + 1)
  boot_meds <- numeric(reps)
  for (r in seq_len(reps)) {
    n_blocks <- ceiling(n_days / block_len)
    s_idx <- sample(starts, n_blocks, replace = TRUE)
    day_sel <- unlist(lapply(s_idx, function(i) days[i:min(i + block_len - 1, n_days)]))
    sel <- day_id %in% day_sel
    boot_meds[r] <- median(diffs[sel], na.rm = TRUE)
  }
  # effective sample size via lag-1 autocorrelation of per-day mean diffs
  day_means <- tapply(diffs, day_id, mean, na.rm = TRUE)
  ac1 <- suppressWarnings(cor(day_means[-length(day_means)], day_means[-1], use = "complete.obs"))
  ac1 <- if (is.finite(ac1)) ac1 else 0
  eff_n <- n_days * (1 - ac1) / (1 + ac1)
  c(lo = unname(quantile(boot_meds, 0.025, na.rm = TRUE)),
    hi = unname(quantile(boot_meds, 0.975, na.rm = TRUE)),
    eff_n = eff_n)
}

# ---- paired comparison on a matched matrix (rows = matched units, cols =
# sites), climatology-cell or exact-timestamp level.
paired_compare <- function(mat, sites_present, day_id_for_bootstrap = NULL) {
  n_sites <- length(sites_present)
  if (n_sites < 2) return(NULL)
  if (n_sites == 2) {
    x <- mat[, 1]; y <- mat[, 2]
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < 4) return(NULL)
    wt <- tryCatch(wilcox.test(x, y, paired = TRUE, conf.int = TRUE, exact = FALSE), error = function(e) NULL)
    if (is.null(wt)) return(NULL)
    diffs <- x - y
    med_diff <- median(diffs)
    hl <- unname(wt$estimate)
    V <- wt$statistic
    n <- length(diffs)
    rbc <- (2 * unname(V) / (n * (n + 1) / 2)) - 1  # matched-pairs rank-biserial from V (sum of positive ranks)
    p <- wt$p.value
    dir <- if (med_diff > 0) paste0(sites_present[1], ">", sites_present[2])
           else if (med_diff < 0) paste0(sites_present[2], ">", sites_present[1]) else "equal"
    boot <- if (!is.null(day_id_for_bootstrap)) block_bootstrap_ci(diffs, day_id_for_bootstrap[ok], block_len = 1) else c(lo = NA, hi = NA, eff_n = NA)
    data.frame(test = "wilcoxon_signed_rank", n_sites = 2, n = n, statistic = unname(V), p_value = p,
               kendall_w = NA_real_, signed_effect_median_diff = med_diff, hl_shift = hl,
               rank_biserial = rbc, direction = dir, extreme_pair = NA_character_,
               ci_low_medboot = boot["lo"], ci_high_medboot = boot["hi"], eff_n_boot = boot["eff_n"])
  } else {
    ok <- complete.cases(mat)
    m <- mat[ok, , drop = FALSE]
    if (nrow(m) < 4) return(NULL)
    ft <- tryCatch(friedman.test(m), error = function(e) NULL)
    if (is.null(ft)) return(NULL)
    k <- ncol(m); N <- nrow(m)
    W <- unname(ft$statistic) / (N * (k - 1))  # Kendall's W from Friedman chi-sq
    meds <- apply(m, 2, median)
    site_pairs <- combn(sites_present, 2, simplify = FALSE)
    dd <- vapply(site_pairs, function(p) abs(meds[[which(sites_present == p[1])]] - meds[[which(sites_present == p[2])]]), numeric(1))
    best <- site_pairs[[which.max(dd)]]
    i1 <- which(sites_present == best[1]); i2 <- which(sites_present == best[2])
    hi <- if (meds[i1] >= meds[i2]) best[1] else best[2]
    lo <- if (meds[i1] >= meds[i2]) best[2] else best[1]
    data.frame(test = "friedman", n_sites = k, n = N, statistic = unname(ft$statistic), p_value = ft$p.value,
               kendall_w = W, signed_effect_median_diff = NA_real_, hl_shift = NA_real_,
               rank_biserial = NA_real_, direction = paste0(hi, ">", lo),
               extreme_pair = paste0(hi, " vs ", lo),
               ci_low_medboot = NA_real_, ci_high_medboot = NA_real_, eff_n_boot = NA_real_)
  }
}

rows <- list()
for (bin in all_bins) {
  sites_here <- SITES[sapply(SITES, function(s) any(site_bin[[s]]$px$elev_bin == bin))]
  if (length(sites_here) < 2) next
  for (v in VARS) {
    # ---- PRIMARY: climatology, matched by (month,hour) cell ----
    clims <- lapply(sites_here, function(s) climatology(s, site_bin_hourly(s, bin, v)))
    names(clims) <- sites_here
    present <- sites_here[!sapply(clims, is.null)]
    if (length(present) < 2) next
    common_cells <- Reduce(intersect, lapply(clims[present], function(d) d$cell))
    if (length(common_cells) < 4) next
    mat <- sapply(present, function(s) {
      d <- clims[[s]]; d$value[match(common_cells, d$cell)]
    })
    day_id_cells <- substr(common_cells, 1, 2)  # month, as a coarse block unit for the climatology-level bootstrap
    res <- paired_compare(mat, present, day_id_for_bootstrap = day_id_cells)
    if (!is.null(res)) {
      rows[[length(rows) + 1]] <- cbind(elev_bin_lo = bin, elev_bin_hi = bin + ELEV_BIN_M, variable = v,
                                         pairing = "climatology_month_hour", sites = paste(present, collapse = ","), res)
    }
  }
}
between_elev_df <- do.call(rbind, rows)
write.csv(between_elev_df, file.path(OUTDIR, "step2_between_sites_by_elevation_climatology.csv"), row.names = FALSE)
cat("\n== Step 2 (PRIMARY, climatology-paired): elevation-bin comparison ==\n")
print(between_elev_df, row.names = FALSE)
cat("\nSaved:", file.path(OUTDIR, "step2_between_sites_by_elevation_climatology.csv"), "\n")
cat("\nElevation bins present (>=2 sites):\n")
for (bin in all_bins) {
  sites_here <- SITES[sapply(SITES, function(s) any(site_bin[[s]]$px$elev_bin == bin))]
  if (length(sites_here) >= 2) cat(sprintf("  [%d,%d): %s\n", bin, bin + ELEV_BIN_M, paste(sites_here, collapse = ",")))
}
