# plot_decoupling_heatmap.R -- heatmap of Spearman rho between shortwave
# radiation and (a) relative humidity, (b) temperature, by site and
# relative-height band, from climate_decoupling_by_height_band.csv
# (climate_decoupling.R). Sites are ordered by footprint elevation where one
# is available. Diverging scale, neutral at rho = 0, fixed to [-1, 1].
# Usage: Rscript plot_decoupling_heatmap.R [binning]   (default relative_decile)
suppressPackageStartupMessages(library(ggplot2))
source("scripts/02_model/config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
binning <- if (length(args) >= 1) args[1] else "relative_decile"

d <- read.csv(file.path(OUTPUT_DIR, "climate_decoupling_by_height_band.csv"), stringsAsFactors = FALSE)
d <- d[d$binning == binning, ]
if (!nrow(d)) stop("No rows for binning '", binning, "'")

# Height bands in numeric order of their lower bound.
lower <- as.numeric(sub("^[\\[(]([^,]+),.*$", "\\1", d$height_band))
band_levels <- unique(d$height_band[order(lower)])
band_labels <- gsub("[][()]", "", band_levels)
band_labels <- sub(",", "–", band_labels)
d$band <- factor(d$height_band, levels = band_levels, labels = band_labels)

# Sites low -> high elevation (bottom -> top); sites without a footprint elevation go last.
elev_path <- file.path(OUTPUT_DIR, "footprint_elevation.csv")
elev <- if (file.exists(elev_path)) read.csv(elev_path, stringsAsFactors = FALSE) else data.frame(site = character(), elev_mean = numeric())
sites <- unique(d$site)
e <- elev$elev_mean[match(sites, elev$site)]
ord <- order(is.na(e), e)
site_lab <- ifelse(is.na(e), sprintf("%s (elevation n/a)", sites),
                   sprintf("%s (%s m)", sites, trimws(format(round(e), big.mark = ","))))
d$site_f <- factor(d$site, levels = sites[ord], labels = site_lab[ord])

pair_lab <- c("swdown~relhum" = "Radiation vs relative humidity", "swdown~temp" = "Radiation vs temperature")
d$pair_f <- factor(ifelse(d$pair %in% names(pair_lab), pair_lab[d$pair], d$pair), levels = unique(c(pair_lab, d$pair)))
d$pair_f <- droplevels(d$pair_f)
d$label <- ifelse(is.na(d$rho), "", sub("^-", "−", sprintf("%.2f", d$rho)))

p <- ggplot(d, aes(x = band, y = site_f, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.7, colour = "#1a1a19") +
  facet_wrap(~pair_f, nrow = 1) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
                       limits = c(-1, 1), na.value = "white", name = "Spearman ρ",
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = if (binning == "relative_decile") "Relative height in canopy (0 = ground, 1 = canopy top)" else "Height above ground (m)",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "grey30"),
        axis.text.y = element_text(colour = "grey15"),
        strip.text = element_text(face = "bold", hjust = 0, size = 10),
        legend.position = "right", legend.key.height = grid::unit(1.1, "cm"),
        panel.spacing = grid::unit(1.2, "lines"))

out <- file.path(OUTPUT_DIR, sprintf("climate_decoupling_heatmap_%s.png", binning))
ggsave(out, plot = p, width = 10.5, height = 0.5 * length(sites) + 1.6, dpi = 300, bg = "white")
cat("Saved:", out, "\n")
cat(sprintf("Cells: %d with rho, %d empty (no daylight data in that band)\n", sum(!is.na(d$rho)), sum(is.na(d$rho))))
