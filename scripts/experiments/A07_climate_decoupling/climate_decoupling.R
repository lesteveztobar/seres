# Task 2 -- microclimate decoupling analysis: does radiation covary with
# humidity/temperature through the vertical profile? Tests the assumption
# in Petter et al. (2021) that humidity can be treated as covarying with
# light.
#
# METHODOLOGICAL NOTE, RECONSIDERED 2026-09-09 (Phase D) -- NOT repointed
# at the per-pixel footprint system, and here is why that was checked, not
# assumed: this test needs JOINT (temp, relhum, swdown) values at the SAME
# pixel-hour to say anything about covariation. Neither the existing
# per-pixel quantile cache (voxel_quantiles) NOR the new per-pixel MEAN
# cache (pixel_means, Phase A) preserves that -- both store each variable
# as an INDEPENDENT marginal reduction, keyed by (pixel, month, daypart),
# with no joint structure across variables at the same original hour. Only
# `lookup_climate_by_height()`'s whole-raster-pooled hourly table (used
# here) keeps temp/relhum/swdown as columns of the same row, genuinely
# joint in time -- so THIS script keeps using it, not because repointing
# was skipped, but because the alternative would silently discard the one
# property the test needs. Fixing this properly would mean retaining real
# per-pixel HOURLY records at write time (feasible per the project's own
# Q3 pricing estimate, ~450-770MB/site -- but not what Phase A actually
# implemented; Phase A wrote per-pixel MEANS, not raw hourly arrays).
# Flagging as a genuine follow-up, not doing it here.
#
# What WAS fixed here (2026-09-09, Phase D): radiation restricted to
# daylight hours (swdown>0) EVERYWHERE, not just as one optional breakdown
# -- night-time swdown=0 is a degenerate value that was diluting Part b/c's
# "overall"/"by month" correlations before this fix. Height reported BOTH
# as absolute 5m bins and relative-height deciles; the unsourced
# rel_height<0.3 crown-base split is dropped entirely. Great-circle
# footprint-centroid distances and DTM-sampled footprint elevation added
# (shared with climate_variation_relative_height.R via the same
# .build_site_climate_series_footprint() call, used here only for that
# metadata -- not for the covariation test itself, per the note above).
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

EXCLUDE_SALOYA <- Sys.getenv("EXCLUDE_SALOYA", unset = "0") == "1"
out_suffix <- if (EXCLUDE_SALOYA) "_noSaloya" else ""
SITES <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu", "LaElenita")
if (EXCLUDE_SALOYA) SITES <- setdiff(SITES, "Saloya")
niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
niches <- augment_elevation(niches)

ceiling_for <- function(site, heights) site_canopy_ceiling(site, niches, max(heights))

# ── Great-circle distances + footprint elevation (metadata, shared helper) ──
fc <- .build_site_climate_series_footprint(SITES, niches)
cent <- fc$footprint_centroid
site_names_fp <- names(cent)
dist_mat <- matrix(NA_real_, length(site_names_fp), length(site_names_fp),
                   dimnames = list(site_names_fp, site_names_fp))
for (i in site_names_fp) for (j in site_names_fp) {
  dist_mat[i, j] <- .haversine_km(cent[[i]]["lon"], cent[[i]]["lat"], cent[[j]]["lon"], cent[[j]]["lat"])
}
cat("=== Footprint centroid distances (km) -- shared with climate_variation_relative_height.R ===\n")
print(round(dist_mat, 2))
elev_df <- do.call(rbind, lapply(names(fc$footprint_elevation), function(s) {
  e <- fc$footprint_elevation[[s]]
  data.frame(site = s, elev_min = e["min"], elev_mean = e["mean"], elev_max = e["max"])
}))
rownames(elev_df) <- NULL
cat("\n=== Footprint elevation (replaces site_elev) ===\n")
print(elev_df, row.names = FALSE)

# ── Joint (temp, relhum, swdown) table -- pooled, per the note above ───────
all_rows <- list()
for (site in SITES) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site))
  if (!file.exists(microenv_path)) { message("Skipping ", site); next }
  message("Reading ", site, "...")
  microenv <- readRDS(microenv_path)
  heights <- microenv_heights(microenv)
  ceiling <- ceiling_for(site, heights)

  n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))
  px <- do.call(rbind, parallel::mclapply(heights, function(h) {
    cl <- lookup_climate_by_height(h, microenv)
    if (is.null(cl)) return(NULL)
    data.frame(height = h, rel_height = h / ceiling, month = cl$month,
              temp = cl$temp, relhum = cl$relhum, swdown = cl$swdown)
  }, mc.cores = n_cores))
  px$site <- site
  px$daypart <- ifelse(px$swdown > 0, "day", "night")
  # Height reported both ways -- absolute 5m bins AND relative-height
  # deciles, side by side; the old rel_height<0.3 split is gone.
  px$height_band_abs <- cut(px$height, breaks = seq(0, ceiling(max(px$height) / 5) * 5, by = 5), include.lowest = TRUE)
  px$height_band_rel <- cut(px$rel_height, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  all_rows[[site]] <- px
}
all_df <- do.call(rbind, all_rows)
saveRDS(all_df, file.path(PROCESSED_DIR, sprintf("climate_decoupling_table%s.rds", out_suffix)))
cat("\nTotal rows:", nrow(all_df), "across", length(all_rows), "sites\n\n")

# 2026-09-09: DAYLIGHT-ONLY radiation everywhere from here on -- swdown=0 at
# night is not "low radiation," it's a structurally different (zero,
# uninformative) regime that was diluting every "overall"/"by month"
# correlation below before this fix. n therefore drops to roughly half of
# the previous (all-hours) version -- reported explicitly per pair below.
day_df <- all_df[all_df$daypart == "day", ]

spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10 || length(unique(x[ok])) < 3 || length(unique(y[ok])) < 3) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

cat("=== Part b: per-site swdown~relhum and swdown~temp (DAYLIGHT hours only, all months) ===\n")
site_overall <- do.call(rbind, lapply(SITES, function(s) {
  d <- day_df[day_df$site == s, ]
  if (nrow(d) == 0) return(NULL)
  r1 <- spearman_safe(d$swdown, d$relhum)
  r2 <- spearman_safe(d$swdown, d$temp)
  data.frame(site = s, pair = c("swdown~relhum", "swdown~temp"),
             rho = c(r1["rho"], r2["rho"]), p = c(r1["p"], r2["p"]), n = c(r1["n"], r2["n"]))
}))
print(site_overall, row.names = FALSE)
write.csv(site_overall, file.path(OUTPUT_DIR, sprintf("climate_decoupling_overall%s.csv", out_suffix)), row.names = FALSE)

cat("\n=== Part c: breakdown by month (daylight only) ===\n")
by_month <- do.call(rbind, lapply(SITES, function(s) {
  d <- day_df[day_df$site == s, ]
  do.call(rbind, lapply(sort(unique(d$month)), function(m) {
    dm <- d[d$month == m, ]
    r1 <- spearman_safe(dm$swdown, dm$relhum)
    r2 <- spearman_safe(dm$swdown, dm$temp)
    data.frame(site = s, month = m, pair = c("swdown~relhum", "swdown~temp"),
               rho = c(r1["rho"], r2["rho"]), p = c(r1["p"], r2["p"]), n = c(r1["n"], r2["n"]))
  }))
}))
print(by_month, row.names = FALSE)
write.csv(by_month, file.path(OUTPUT_DIR, sprintf("climate_decoupling_by_month%s.csv", out_suffix)), row.names = FALSE)

cat("\n=== Part c: breakdown by ABSOLUTE height band (5m bins, daylight only) ===\n")
by_band_abs <- do.call(rbind, lapply(SITES, function(s) {
  d <- day_df[day_df$site == s, ]
  do.call(rbind, lapply(levels(d$height_band_abs), function(b) {
    db <- d[!is.na(d$height_band_abs) & d$height_band_abs == b, ]
    r1 <- spearman_safe(db$swdown, db$relhum)
    r2 <- spearman_safe(db$swdown, db$temp)
    data.frame(site = s, height_band = b, binning = "absolute_5m",
               pair = c("swdown~relhum", "swdown~temp"),
               rho = c(r1["rho"], r2["rho"]), p = c(r1["p"], r2["p"]), n = c(r1["n"], r2["n"]))
  }))
}))
cat("\n=== Part c: breakdown by RELATIVE height decile (daylight only) ===\n")
by_band_rel <- do.call(rbind, lapply(SITES, function(s) {
  d <- day_df[day_df$site == s, ]
  do.call(rbind, lapply(levels(d$height_band_rel), function(b) {
    db <- d[!is.na(d$height_band_rel) & d$height_band_rel == b, ]
    r1 <- spearman_safe(db$swdown, db$relhum)
    r2 <- spearman_safe(db$swdown, db$temp)
    data.frame(site = s, height_band = b, binning = "relative_decile",
               pair = c("swdown~relhum", "swdown~temp"),
               rho = c(r1["rho"], r2["rho"]), p = c(r1["p"], r2["p"]), n = c(r1["n"], r2["n"]))
  }))
}))
by_band <- rbind(by_band_abs, by_band_rel)
print(by_band, row.names = FALSE)
write.csv(by_band, file.path(OUTPUT_DIR, sprintf("climate_decoupling_by_height_band%s.csv", out_suffix)), row.names = FALSE)

# ── Figures: one per site, one cross-site summary ───────────────────────────
library(ggplot2)
for (s in SITES) {
  d <- day_df[day_df$site == s, ]
  if (nrow(d) == 0) next
  p1 <- ggplot(d, aes(x = swdown, y = relhum, colour = height_band_rel)) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(title = paste0(s, " -- swdown vs relhum (daylight)"), x = "swdown", y = "relhum (%)", colour = "rel. height decile") +
    theme_minimal(base_size = 11)
  p2 <- ggplot(d, aes(x = swdown, y = temp, colour = height_band_rel)) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(title = paste0(s, " -- swdown vs temp (daylight)"), x = "swdown", y = "temp (C)", colour = "rel. height decile") +
    theme_minimal(base_size = 11)
  p <- patchwork::wrap_plots(p1, p2, nrow = 1) + patchwork::plot_layout(guides = "collect")
  out_path <- file.path(OUTPUT_DIR, sprintf("climate_decoupling_%s%s.png", s, out_suffix))
  ggsave(out_path, plot = p, width = 12, height = 5, dpi = 300, bg = "white")
  message("Saved: ", out_path)
}

summary_df <- site_overall
p_summary <- ggplot(summary_df, aes(x = site, y = rho, fill = pair)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Cross-site summary: swdown~relhum / swdown~temp (daylight, Spearman rho)",
       y = "Spearman rho", x = NULL) +
  theme_minimal(base_size = 12)
out_summary <- file.path(OUTPUT_DIR, sprintf("climate_decoupling_cross_site_summary%s.png", out_suffix))
ggsave(out_summary, plot = p_summary, width = 9, height = 6, dpi = 300, bg = "white")
message("Saved: ", out_summary)
cat("\nDone.\n")
