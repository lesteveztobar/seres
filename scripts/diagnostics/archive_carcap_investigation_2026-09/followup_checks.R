source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
library(MASS)

SITES <- c("Maquipucuna","Mashpi","MindoMirador","MindoTarabita","Saloya","Yanayacu")
niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
niches <- augment_elevation(niches)
sc <- .build_site_climate_series(SITES, niches = niches)
site_pixels <- sc$site_pixels
site_elev <- sc$site_elev

# v7 (Phase 1.6): shared ceiling formula, site_canopy_ceiling() (shared_helpers.R) --
# was an independently-duplicated mean(...)-based formula here.
ceiling_for <- function(site, heights) site_canopy_ceiling(site, niches, max(heights))

cat("#################### CHECK 3: Task 2 height confound ####################\n")
for (s in names(site_pixels)) {
  px <- site_pixels[[s]]
  px$daypart <- ifelse(px$swdown > 0, "day", "night")
  cat(sprintf("\n=== %s ===\n", s))
  for (dp in c("day","night")) {
    d <- px[px$daypart == dp, ]
    heights_here <- sort(unique(d$height))
    per_tier_rho <- sapply(heights_here, function(h) {
      dd <- d[d$height == h, ]
      if (nrow(dd) < 10 || length(unique(dd$swdown)) < 3 || length(unique(dd$relhum)) < 3) return(NA_real_)
      suppressWarnings(cor(dd$swdown, dd$relhum, method = "spearman"))
    })
    per_tier_rho <- per_tier_rho[is.finite(per_tier_rho)]
    if (length(per_tier_rho) == 0) {
      cat(sprintf("  %s: no usable height tiers\n", dp)); next
    }
    # Partial Spearman (swdown~relhum | height): rank-transform all 3, then
    # partial correlation via residuals of linear fits on ranked height.
    ok <- is.finite(d$swdown) & is.finite(d$relhum) & is.finite(d$height)
    dd <- d[ok, ]
    if (nrow(dd) > 20 && length(unique(dd$height)) > 2) {
      r_sw <- rank(dd$swdown); r_rh <- rank(dd$relhum); r_h <- rank(dd$height)
      resid_sw <- resid(lm(r_sw ~ r_h))
      resid_rh <- resid(lm(r_rh ~ r_h))
      partial_rho <- cor(resid_sw, resid_rh)
    } else partial_rho <- NA_real_
    cat(sprintf("  %s: within-tier rho median=%.3f IQR=[%.3f,%.3f] range=[%.3f,%.3f] (n_tiers=%d) | partial rho (height-controlled)=%.3f\n",
        dp, median(per_tier_rho), quantile(per_tier_rho,0.25), quantile(per_tier_rho,0.75),
        min(per_tier_rho), max(per_tier_rho), length(per_tier_rho), partial_rho))
  }
}

cat("\n\n#################### ADDITIONAL: Task 1 slopes split by daypart ####################\n")
rows <- list()
for (s in names(site_pixels)) {
  px <- site_pixels[[s]]
  ceiling <- ceiling_for(s, px$height)
  px$rel_h <- px$height / ceiling
  px$daypart <- ifelse(px$swdown > 0, "day", "night")
  for (dp in c("day","night")) {
    d <- px[px$daypart == dp, ]
    for (v in c("temp","relhum","swdown")) {
      if (dp == "night" && v == "swdown") next  # swdown=0 by construction at night
      y <- d[[v]]; x <- d$rel_h
      ok <- is.finite(y) & is.finite(x)
      yy <- y[ok]; xx <- x[ok]
      if (length(yy) < 20 || length(unique(xx)) < 3) next
      rfit <- tryCatch(MASS::rlm(yy ~ xx, maxit=100), error=function(e) NULL)
      if (is.null(rfit)) next
      co <- summary(rfit)$coefficients
      rho <- suppressWarnings(cor(xx, yy, method="spearman"))
      rows[[length(rows)+1]] <- data.frame(site=s, daypart=dp, variable=v, n=length(yy),
        slope=co["xx","Value"], se=co["xx","Std. Error"], spearman_rho=rho)
    }
  }
}
daypart_slopes <- do.call(rbind, rows)
print(daypart_slopes, row.names = FALSE)
write.csv(daypart_slopes, file.path(OUTPUT_DIR, "elevation_test_relative_height_slopes_by_daypart.csv"), row.names=FALSE)
cat("\ntemp slope sign by site, day vs night:\n")
temp_rows <- daypart_slopes[daypart_slopes$variable=="temp",]
print(reshape(temp_rows[,c("site","daypart","slope")], idvar="site", timevar="daypart", direction="wide"), row.names=FALSE)

cat("\n\n#################### ADDITIONAL: VPD-based Task 4b ####################\n")
crossing_height_vpd <- function(px, target, ceiling) {
  px$es <- 0.6108 * exp(17.27 * px$temp / (px$temp + 237.3))
  px$vpd <- px$es * (1 - px$relhum / 100)
  agg <- aggregate(vpd ~ height, data = px, FUN = median)
  agg <- agg[order(agg$height), ]
  d <- agg$vpd - target
  sign_change <- which(diff(sign(d)) != 0)
  if (length(sign_change) == 0) return(c(abs = NA_real_, rel = NA_real_))
  i <- sign_change[1]
  h1 <- agg$height[i]; h2 <- agg$height[i+1]; d1 <- d[i]; d2 <- d[i+1]
  h_cross <- h1 + (0 - d1) * (h2 - h1) / (d2 - d1)
  c(abs = h_cross, rel = h_cross / ceiling)
}
all_vpd <- unlist(lapply(names(site_pixels), function(s) {
  px <- site_pixels[[s]]
  es <- 0.6108 * exp(17.27 * px$temp / (px$temp + 237.3))
  es * (1 - px$relhum/100)
}))
target_vpd <- median(all_vpd, na.rm=TRUE)
cat("Target VPD (pooled median across all sites):", round(target_vpd,4), "kPa\n\n")
vpd_rows <- lapply(names(site_pixels), function(s) {
  ceiling <- ceiling_for(s, site_pixels[[s]]$height)
  cr <- crossing_height_vpd(site_pixels[[s]], target_vpd, ceiling)
  data.frame(site=s, elevation=site_elev[s], ceiling=ceiling,
             crossing_height_abs=cr["abs"], crossing_height_rel=cr["rel"])
})
vpd_df <- do.call(rbind, vpd_rows)
n_usable <- sum(!is.na(vpd_df$crossing_height_abs))
cat("n sites with a usable VPD crossing height:", n_usable, "of", nrow(vpd_df), "\n")
print(vpd_df, row.names=FALSE)
write.csv(vpd_df, file.path(OUTPUT_DIR, "elevation_canopy_crossing_height_vpd.csv"), row.names=FALSE)
if (n_usable >= 4) {
  d <- vpd_df[!is.na(vpd_df$crossing_height_abs),]
  cat("\nSpearman (elevation vs VPD-crossing-height):\n")
  print(cor.test(d$elevation, d$crossing_height_abs, method="spearman"))
} else {
  cat("n<4 -- not fitting a regression, per the same rule as the relhum version.\n")
}
cat("\nDone.\n")
