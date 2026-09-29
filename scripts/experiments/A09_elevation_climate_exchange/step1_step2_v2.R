# step1_step2_v2.R -- 2026-09-28. Derived variables (VPD) + elevation-
# matched between-site comparison, redone per user's confirmed methodology
# after the Step-0 audit (see docs/microclimate_reanalysis_v2/):
#   - site-level hourly (spatially-pooled) series -- no true per-pixel
#     hourly data exists in this pipeline (Step 0.0 finding).
#   - elevation bins: tiers aggregated to ONE VALUE PER SITE PER TIMESTAMP
#     (median across the tiers landing in that bin), fixing the tier-count
#     pseudoreplication in the old n=114,192/263,520.
#   - between-site pairing: PRIMARY = climatology (month x hour-of-day cell
#     medians per site), since only 9/21 site pairs share any exact
#     calendar timestamp (Step 0 finding); exact-date pairing is a
#     robustness check restricted to those 9 pairs.
#   - swdown clipped to >=0 (negative-swdown audit reported separately).
suppressPackageStartupMessages({})
source("scripts/02_model/config/paths.R")

RUN_TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(OUTPUT_DIR, paste0("microclimate_reanalysis_", RUN_TS))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
cat("Writing to:", OUTDIR, "\n")

cache <- readRDS(file.path(PROCESSED_DIR, "site_climate_cache_v2.rds"))
site_pixels <- cache$site_pixels   # list(site -> df(height, temp, relhum, swdown)), swdown already clipped >=0
site_elev   <- cache$site_elev
tme_by_site <- cache$tme_by_site
SITES <- names(site_pixels)

# ---- Step 1.1: VPD per timestamp, per height tier, per site --------------
# es = 0.6108*exp(17.27*T/(T+237.3)) kPa; VPD = es*(1-RH/100), RH clipped [0,100].
rh_clip_log <- list()
for (s in SITES) {
  px <- site_pixels[[s]]
  rh_raw <- px$relhum
  n_clip_hi <- sum(rh_raw > 100, na.rm = TRUE)
  n_clip_lo <- sum(rh_raw < 0, na.rm = TRUE)
  rh <- pmin(pmax(rh_raw, 0), 100)
  es <- 0.6108 * exp(17.27 * px$temp / (px$temp + 237.3))
  px$vpd_kPa <- es * (1 - rh / 100)
  site_pixels[[s]] <- px
  rh_clip_log[[s]] <- data.frame(site = s, n_clipped_hi = n_clip_hi, n_clipped_lo = n_clip_lo, n_total = nrow(px))
}
rh_clip_df <- do.call(rbind, rh_clip_log)
write.csv(rh_clip_df, file.path(OUTDIR, "step1_1_rh_clip_audit.csv"), row.names = FALSE)
cat("\n== Step 1.1: RH clip audit ==\n"); print(rh_clip_df, row.names = FALSE)

# ---- Step 1.2: daytime mask, sensitivity to threshold ---------------------
# site_pixels rows are (height x hour) long-format; each site has its own
# tme recycled once per height tier (same 4392-length vector repeated).
day_thresholds <- c(1, 10, 50)
day_sens <- do.call(rbind, lapply(SITES, function(s) {
  px <- site_pixels[[s]]
  do.call(rbind, lapply(day_thresholds, function(th) {
    data.frame(site = s, threshold = th, frac_day = mean(px$swdown > th, na.rm = TRUE))
  }))
}))
write.csv(day_sens, file.path(OUTDIR, "step1_2_daytime_threshold_sensitivity.csv"), row.names = FALSE)
cat("\n== Step 1.2: daytime-mask sensitivity ==\n"); print(day_sens, row.names = FALSE)
DAY_THRESH <- 10  # primary, per spec

# ---- Step 1.3: relative transmittance vs topmost tier, daytime only ------
transmit_rows <- list()
for (s in SITES) {
  px <- site_pixels[[s]]
  top_h <- max(px$height)
  n_hr <- nrow(px) / length(unique(px$height))
  # px is stacked height-major (rbind over heights, each block length 4392);
  # recover hour index within each height block.
  px$hour_idx <- ave(seq_len(nrow(px)), px$height, FUN = seq_along)
  top_sw <- px$swdown[px$height == top_h]
  names(top_sw) <- seq_along(top_sw)
  px$sw_top <- top_sw[px$hour_idx]
  px_day <- px[px$sw_top > DAY_THRESH, ]
  px_day$transmittance <- ifelse(px_day$sw_top > 0, px_day$swdown / px_day$sw_top, NA_real_)
  agg <- aggregate(transmittance ~ height, data = px_day, FUN = function(x) median(x, na.rm = TRUE))
  agg$site <- s; agg$top_height_m <- top_h
  transmit_rows[[s]] <- agg
}
transmit_df <- do.call(rbind, transmit_rows)
write.csv(transmit_df, file.path(OUTDIR, "step1_3_relative_transmittance.csv"), row.names = FALSE)
cat("\n== Step 1.3: relative transmittance (median per height, daytime only) -- written, n rows =", nrow(transmit_df), "==\n")

saveRDS(site_pixels, file.path(OUTDIR, "site_pixels_with_vpd.rds"))
cat("\nStep 1 done. site_pixels (with vpd_kPa) saved to", file.path(OUTDIR, "site_pixels_with_vpd.rds"), "\n")
