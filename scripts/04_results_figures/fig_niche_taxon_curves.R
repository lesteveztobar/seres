# fig_niche_taxon_curves.R --
# Outputs: fig_niche_by_site_curves.{png,pdf}   canopy niche suitability heatmap by site (weighted richness)
#          fig_niche_multisite_curves.{png,pdf} curves for species identified at 2+ sites
#   sbatch scripts/04_results_figures/run_fig.sh fig_niche_taxon_curves.R
source("scripts/04_results_figures/fig_common.R")
suppressMessages(source("scripts/02_model/plots/plot_functions.R"))
suppressPackageStartupMessages(library(ggrepel))

SMOOTH_M <- 2.0; PEAK_MIN <- 5
rows <- .niche_score_rows_all_sites_cached(PROCESSED_DIR, 0.4)$rows
df <- data.frame(site = as.character(rows$site), species = as.character(rows$species),
                 height = rows$height, score = rows$after)
df <- df[is.finite(df$score), ]
df <- do.call(rbind, lapply(split(df, list(df$site, df$species), drop = TRUE), function(d) {
  d <- d[order(d$height), ]
  d$smooth <- vapply(d$height, function(h) mean(d$score[abs(d$height - h) <= SMOOTH_M / 2 + 1e-9]), numeric(1))
  d
}))
is_morpho <- function(x) grepl("^[A-Z][a-z]+\\s+sp", x)
ident <- df[!is_morpho(df$species), ]
n_sites <- tapply(ident$site, ident$species, function(s) length(unique(s)))
multi <- names(n_sites)[n_sites >= 2]
cat("Identified species:", length(n_sites), "| at 2+ sites:", paste(multi, collapse = "; "), "\n")

# One colour per species, fixed across both figures: the multi-site species
# take the four most distinct Okabe-Ito hues, the rest Paul Tol's muted set.
okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
tol   <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
           "#CC6677", "#882255", "#AA4499", "#6699CC", "#AA4466", "#4477AA")
others <- setdiff(sort(names(n_sites)), multi)
TAXON_COL <- c(setNames(okabe[seq_along(multi)], sort(multi)), setNames(tol[seq_along(others)], others))

make_plot <- function(d, ncol, top = 115) {
  d$site <- factor(d$site, levels = intersect(SITE_LEVELS, d$site))
  pk <- do.call(rbind, lapply(split(d, list(d$site, d$species), drop = TRUE), function(x) x[which.max(x$smooth), ]))
  pk <- pk[pk$smooth >= PEAK_MIN, ]
  pk$lab <- taxon_expr(pk$species, abbreviate = TRUE)
  sp <- sort(unique(d$species))
  ggplot(d, aes(x = height, y = smooth, colour = species)) +
    geom_line(linewidth = 0.6) +
    geom_point(data = pk, aes(fill = species), shape = 21, colour = "black", size = 1.7, stroke = STROKE_W,
               show.legend = FALSE) +
    geom_text_repel(data = pk, aes(label = lab), parse = TRUE, size = 2.3 + TXT_UP + 1 / ggplot2::.pt, show.legend = FALSE,
                    min.segment.length = 0, segment.size = 0.25, box.padding = 0.45, seed = 1,
                    max.overlaps = Inf, direction = "both", ylim = c(109, top + 4), nudge_y = 3, force = 3) +
    scale_colour_manual(values = TAXON_COL[sp], labels = taxon_labels(TRUE), name = "Taxon unit") +
    scale_fill_manual(values = TAXON_COL[sp], guide = "none") +
    scale_y_continuous(limits = c(0, top), breaks = c(0, 50, 100)) +
    facet_wrap(~site, ncol = ncol, scales = "free_x") +
    labs(x = "Height above ground (m)", y = "Climate-niche suitability (0–100)") +
    theme_seres() +
    theme(strip.text = element_text(face = "bold", size = 9.5), legend.position = "bottom",
          legend.text = element_text(size = 8.5), panel.spacing = unit(4, "mm")) +
    guides(colour = guide_legend(nrow = 3, override.aes = list(linewidth = 1)))
}

# fig_niche_by_site_curves reimagined as a heatmap:
# Color represents how many species that height tier is suitable for,
# weighted by the suitability score (sum of after-ceiling-rescale score / 100).
elev <- site_elevation()
sites_ordered <- intersect(names(sort(elev)), unique(df$site))
site_labels <- sprintf("%s\n(%s m)", sites_ordered, trimws(format(round(elev[sites_ordered]), big.mark = ",")))
names(site_labels) <- sites_ordered

agg_all <- aggregate(smooth ~ site + height, data = df, FUN = function(x) sum(x / 100))
names(agg_all)[3] <- "w_species"
agg_all$site_f <- factor(agg_all$site, levels = sites_ordered, labels = site_labels[sites_ordered])

p1 <- ggplot(agg_all, aes(x = site_f, y = height, fill = w_species)) +
  geom_tile(width = 0.92, height = 0.4) +
  scale_fill_viridis_c(option = "mako", direction = -1,
                       name = "Suitability-weighted\ntaxon-unit richness",
                       limits = c(0, ceiling(max(agg_all$w_species)))) +
  scale_y_continuous(breaks = seq(0, 50, by = 10), limits = c(0, 50), expand = c(0.01, 0)) +
  labs(x = "Site (elevation order)", y = "Height above ground (m)") +
  theme_seres(base_size = 10) +
  theme(axis.text.x = element_text(size = 8.5, lineheight = 0.95),
        legend.position = "right",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.width = unit(3.5, "mm"),
        legend.key.height = unit(12, "mm"),
        panel.grid.major.x = element_blank())

save_fig("fig_niche_by_site_curves", p1, width = 190, height = 95)
p2 <- make_plot(ident[ident$species %in% multi, ], ncol = 2, top = 130)
save_fig("fig_niche_multisite_curves", p2 + guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1))),
         height = 105)
cat("Sites in multisite curves figure:", paste(sort(unique(ident$site[ident$species %in% multi])), collapse = ", "), "\n")
