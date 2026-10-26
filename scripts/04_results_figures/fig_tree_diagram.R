# fig_tree_diagram.R -- 2D side view of the modelled tree geometry, replacing
# the 3D point cloud of plot_tree_diagram() (plot_functions.R), where the
# crown read as a flat disc and the trunk was barely visible. Same geometry:
# zones and crown radius come from the model's own helpers
# (.classify_tree_zone(), .effective_crown_radius(); get_colonization.R).
# Each height slice is drawn as a bar of half-width = trunk radius (zones 1-2)
# or effective crown radius (zones 3-5). No in-figure caption.
#   Rscript scripts/04_results_figures/fig_tree_diagram.R
source("scripts/04_results_figures/fig_common.R")
suppressMessages(source("scripts/02_model/plots/plot_functions.R"))

fp <- list(mean_hgt = 8.4, mean_crown_r = 2.0, trunk_r = 0.114)   # plot_tree_diagram() defaults
RES <- 0.05
h <- seq(RES / 2, fp$mean_hgt, by = RES)
d <- data.frame(h = h, zone = vapply(h / fp$mean_hgt, .classify_tree_zone, numeric(1)))
d$r <- ifelse(d$zone <= 2, fp$trunk_r,
              vapply(h / fp$mean_hgt, function(rh) .effective_crown_radius(fp$mean_crown_r, rh), numeric(1)))
d <- d[d$r > 0, ]
zone_lab <- c(`1` = "1 trunk base", `2` = "2 trunk", `3` = "3 lower crown",
              `4` = "4 mid crown", `5` = "5 upper crown")
zone_col <- c(`1` = "#6b4423", `2` = "#8a6d3b", `3` = "#a8d18d", `4` = "#4f9a4f", `5` = "#2d6a2d")
d$zone <- factor(d$zone, levels = names(zone_lab))

bounds <- do.call(rbind, lapply(split(d, d$zone, drop = TRUE), function(z)
  data.frame(zone = z$zone[1], lo = min(z$h) - RES / 2, hi = max(z$h) + RES / 2)))
cat("Zone height bounds (m):\n"); print(bounds, row.names = FALSE)
r_max <- max(d$r); h_rmax <- d$h[which.max(d$r)]

p <- ggplot(d) +
  geom_rect(aes(xmin = -r, xmax = r, ymin = h - RES / 2, ymax = h + RES / 2, fill = zone)) +
  geom_hline(data = bounds[-1, ], aes(yintercept = lo), colour = "grey55", linewidth = 0.25, linetype = "dashed") +
  annotate("segment", x = 0, xend = r_max, y = h_rmax, yend = h_rmax, colour = "black", linewidth = 0.3,
           arrow = grid::arrow(length = unit(1.5, "mm"), ends = "both")) +
  annotate("text", x = r_max / 2, y = h_rmax + 0.25, size = 2.4 + TXT_UP,
           label = sprintf("crown radius %.1f m", r_max)) +
  scale_fill_manual(values = zone_col, labels = zone_lab, name = "Zone") +
  scale_y_continuous(breaks = sort(unique(round(c(0, bounds$lo, fp$mean_hgt), 1))), expand = c(0, 0.1)) +
  coord_fixed(xlim = c(-fp$mean_crown_r * 1.2, fp$mean_crown_r * 1.2)) +
  labs(x = "Distance from trunk axis (m)", y = "Height above ground (m)") +
  theme_seres() + theme(legend.position = "right", panel.grid.major.x = element_blank())
save_fig("fig_tree_diagram", p, width = HALF_W + 20, height = 110)
