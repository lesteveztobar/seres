# fig_decoupling.R -- daytime covariation of shortwave radiation with relative
# humidity and with temperature, all seven sites, Yanayacu set apart as the
# eastern-slope outgroup.
#   fig_decoupling_heatmap.{png,pdf}  rho by relative-height decile (fig:decouplingdecile)
#   fig_decoupling_summary.{png,pdf}  pooled rho per site, styled like
#                                     fig_night_covariation (fig:decouplingsummary)
# Inputs (climate_decoupling.R): output/tables/climate/climate_decoupling_by_height_band.csv,
#   climate_decoupling_overall.csv. Elevation: output/tables/exports/site_elevation.csv.
#   sbatch scripts/04_results_figures/run_fig.sh fig_decoupling.R
source("scripts/04_results_figures/fig_common.R")

elev <- site_elevation()
site_lab <- function(s) sprintf("%s (%s m)", s, format(round(elev[s]), big.mark = ",", trim = TRUE))
PAIR <- c("swdown~relhum" = "Radiation vs relative humidity", "swdown~temp" = "Radiation vs temperature")

# ── heatmap ──────────────────────────────────────────────────────────────────
d <- read.csv("output/tables/climate/climate_decoupling_by_height_band.csv")
d <- d[d$binning == "relative_decile" & d$pair %in% names(PAIR), ]
lower <- as.numeric(sub("^[\\[(]([^,]+),.*$", "\\1", d$height_band))
lev <- unique(d$height_band[order(lower)])
d$band <- factor(d$height_band, levels = lev, labels = sub(",", "–", gsub("[][()]", "", lev)))
ord <- names(sort(elev[unique(d$site)]))
d$site_f <- factor(site_lab(d$site), levels = site_lab(ord))
d$group <- slope_group(d$site)
d$pair_f <- factor(PAIR[d$pair], levels = PAIR)
d$label <- ifelse(is.na(d$rho), "", sub("^-", "−", sprintf("%.2f", d$rho)))
cat(sprintf("Heatmap: %d cells, %d sites\n", nrow(d), length(ord)))

ph <- ggplot(d, aes(x = band, y = site_f, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 1.75 + TXT_UP - (1 / ggplot2::.pt), colour = "#1a1a19") +
  facet_grid(group ~ pair_f, scales = "free_y", space = "free_y", labeller = slope_short) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
                       limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), name = "Spearman ρ") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Relative height in canopy (0 = ground, 1 = canopy top)", y = NULL) +
  theme_seres() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(face = "bold", hjust = 0, size = 10),
        strip.text.y = element_text(face = "bold", angle = 0, hjust = 0, size = 9),
        strip.background.y = element_rect(fill = "grey93", colour = NA),
        panel.spacing.y = unit(4, "mm"), panel.spacing.x = unit(4, "mm"),
        legend.position = "bottom", legend.key.width = unit(14, "mm"), legend.key.height = unit(3, "mm"))
save_fig("fig_decoupling_heatmap", ph, height = 110)

# ── cross-site summary ───────────────────────────────────────────────────────
o <- read.csv("output/tables/climate/climate_decoupling_overall.csv")
o <- o[o$pair %in% names(PAIR), ]
o$site_f <- factor(site_lab(o$site), levels = site_lab(ord))
o$group <- slope_group(o$site)
o$pair_f <- factor(PAIR[o$pair], levels = PAIR)
seg <- reshape(o[, c("site_f", "group", "pair", "rho")], idvar = c("site_f", "group"), timevar = "pair", direction = "wide")
names(seg) <- sub("rho.swdown~relhum", "rh", sub("rho.swdown~temp", "tp", names(seg)))
ps <- ggplot(o, aes(y = site_f)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_segment(data = seg, aes(x = rh, xend = tp, yend = site_f), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(x = rho, shape = pair_f, fill = pair_f), colour = "black", stroke = STROKE_W, size = 2.6) +
  scale_shape_manual(values = c(21, 24), name = NULL, labels = c("with relative humidity", "with temperature")) +
  scale_fill_manual(values = c("grey20", "white"), name = NULL, labels = c("with relative humidity", "with temperature")) +
  scale_x_continuous(limits = c(-1, 0.5), breaks = seq(-1, 0.5, 0.25)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = slope_short) +
  labs(x = "Spearman ρ with shortwave radiation\n(daylight hours, all height tiers pooled)", y = NULL) +
  theme_seres() +
  theme(legend.position = "top", legend.justification = "left",
        strip.text.y = element_text(face = "bold", angle = 0, hjust = 0, size = 9),
        strip.background.y = element_rect(fill = "grey93", colour = NA),
        panel.spacing.y = unit(3, "mm"))
save_fig("fig_decoupling_summary", ps, width = 140, height = 80)
print(o[order(o$pair, o$rho), c("site", "pair", "rho", "n")], row.names = FALSE)
