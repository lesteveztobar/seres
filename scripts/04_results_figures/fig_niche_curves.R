# fig_niche_curves.R -- climate-niche suitability against height above
# ground, one line per site, one panel per taxon unit, so between-site
# differences in where the suitability peak sits read directly.
# Suitability is the after-ceiling-rescale score (0-100) used by the model
# (same rows as niche_suitability_by_site*.png), smoothed with a centred
# running mean over SMOOTH_M metres of height tiers; the dot marks the
# height of the smoothed peak. No in-figure caption.
# Outputs:
#   fig_niche_curves.{png,pdf}      taxa modelled at 2+ sites (main text)
#   fig_niche_curves_all.{png,pdf}  every taxon unit (appendix)
#   Rscript scripts/04_results_figures/fig_niche_curves.R
source("scripts/04_results_figures/fig_common.R")
suppressMessages(source("scripts/02_model/plots/plot_functions.R"))

SMOOTH_M <- 2.0
built <- .niche_score_rows_all_sites_cached(PROCESSED_DIR, 0.4)
rows <- built$rows
cat(sprintf("Niche score rows: %d (%d sites, %d taxa)\n", nrow(rows),
            length(unique(rows$site)), length(unique(rows$species))))

df <- data.frame(site = as.character(rows$site), species = as.character(rows$species),
                 height = rows$height, score = rows$after)
df <- df[is.finite(df$score), ]

smooth_one <- function(d) {
  d <- d[order(d$height), ]
  half <- SMOOTH_M / 2
  d$smooth <- vapply(d$height, function(h) mean(d$score[abs(d$height - h) <= half + 1e-9]), numeric(1))
  d
}
df <- do.call(rbind, lapply(split(df, list(df$site, df$species), drop = TRUE), smooth_one))
peaks <- do.call(rbind, lapply(split(df, list(df$site, df$species), drop = TRUE),
                               function(d) d[which.max(d$smooth), ]))
# A curve that never rises above PEAK_MIN has no meaningful peak (taxon
# unsuitable at every height at that site): line drawn, no peak marker.
PEAK_MIN <- 5
cat("No peak (max <", PEAK_MIN, "):", paste(peaks$species[peaks$smooth < PEAK_MIN], peaks$site[peaks$smooth < PEAK_MIN], sep = " @ ", collapse = "; "), "\n")
peaks_all <- peaks
peaks <- peaks[peaks$smooth >= PEAK_MIN, ]

sp_order <- .niche_species_order(rows)
n_sites <- tapply(df$site, df$species, function(s) length(unique(s)))
lab <- function(x) ifelse(grepl("^[A-Z][a-z]+ sp", x), x, .niche_sp_abbr(x))  # keep "Maxillaria sp. B" whole

make_plot <- function(keep, ncol) {
  d <- df[df$species %in% keep, ]; pk <- peaks[peaks$species %in% keep, ]
  lev <- unique(taxon_expr(sp_order[sp_order %in% keep], abbreviate = TRUE))
  d$species <- factor(taxon_expr(d$species, abbreviate = TRUE), levels = lev)
  pk$species <- factor(taxon_expr(pk$species, abbreviate = TRUE), levels = lev)
  d$site <- factor(d$site, levels = intersect(SITE_LEVELS, d$site))
  pk$site <- factor(pk$site, levels = levels(d$site))
  ggplot(d, aes(x = height, y = smooth, colour = site)) +
    geom_line(linewidth = 0.6) +
    geom_segment(data = pk, aes(x = height, xend = height, y = 0, yend = smooth),
                 linetype = "dotted", linewidth = 0.3) +
    geom_point(data = pk, aes(fill = site), shape = 21, colour = "black", size = 1.8, stroke = STROKE_W) +
    scale_colour_manual(values = SITE_COL[levels(d$site)], name = "Site") +
    scale_fill_manual(values = SITE_COL[levels(d$site)], guide = "none") +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    facet_wrap(~species, ncol = ncol, labeller = label_parsed) +
    labs(x = "Height above ground (m)", y = "Climate-niche suitability (0–100)") +
    theme_seres() +
    theme(strip.text = element_text(size = 7), legend.position = "bottom",
          panel.spacing = unit(3, "mm"))
}

multi <- names(n_sites)[n_sites >= 2]
cat("Taxa at 2+ sites:", paste(multi, collapse = "; "), "\n")
save_fig("fig_niche_curves", make_plot(multi, ncol = min(3, length(multi))),
         height = 40 + 42 * ceiling(length(multi) / 3))
all_sp <- names(n_sites)
save_fig("fig_niche_curves_all", make_plot(all_sp, ncol = 5),
         height = 30 + 30 * ceiling(length(all_sp) / 5))

pk_out <- peaks_all[order(peaks_all$species, peaks_all$site), c("species", "site", "height", "smooth")]
names(pk_out)[3:4] <- c("peak_height_m", "peak_suitability_smoothed")
write.csv(pk_out, file.path(EXP_DIR, "niche_curve_peaks.csv"), row.names = FALSE)
cat("Saved:", file.path(EXP_DIR, "niche_curve_peaks.csv"), "\n")
print(pk_out, row.names = FALSE)
