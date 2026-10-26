# fig_competition.R -- Results figure 2.2: height shift (joint minus
# isolated) against joint-run sample size, one point per site-taxon pair.
# x = n_multispecies (individuals in the joint run, pooled over replicates).
# Encoding (2026-10-02 redesign): colour = site (SITE_COL) on every point;
# outline always black; Holm significance = full vs faded fill of the same
# site colour; shape = taxon for the four key taxa, circle for all others
# (only shapes 21-25 take a separate fill and outline, so at most four
# taxa can have their own shape next to the "other taxa" circle).
# No in-figure caption -- that text lives in the LaTeX caption.
#   Rscript scripts/04_results_figures/fig_competition.R
source("scripts/04_results_figures/fig_common.R")
suppressPackageStartupMessages(library(ggrepel))

src <- "output/tables/competition/competition_analysis.csv"
d_all <- read.csv(src)
cat(sprintf("%s: %d rows\n", src, nrow(d_all)))

# The manuscript reports 30 pairs at five sites (sec:res-competition). The
# 31st row is Yanayacu Ida sp.: both runs put every individual on the lowest
# tier (mean 0.5 m in both), D ~ 0, so it carries no comparison.
drop <- d_all$site == "Yanayacu" & d_all$species == "Ida sp."
d <- d_all[!drop, ]
cat(sprintf("Dropped %d row(s): Yanayacu / Ida sp.; plotted %d\n", sum(drop), nrow(d)))

d$filled <- d$sig != "ns"                      # Holm-adjusted
d$site <- factor(d$site, levels = intersect(SITE_LEVELS, d$site))

# Key taxa: Holm-significant taxa first (smallest adjusted p), then taxa
# present at the most sites; four shapes available besides "other taxa".
sp <- unique(d$species)
tx <- data.frame(species = sp,
                 p = vapply(sp, function(x) min(d$ks_p_holm[d$species == x]), numeric(1)),
                 n_sites = vapply(sp, function(x) sum(d$species == x), integer(1)))
tx$any_sig <- tx$species %in% d$species[d$filled]
tx <- tx[order(!tx$any_sig, tx$p, -tx$n_sites), ]
key <- head(tx$species, 4)
other_lab <- "Other taxa"
d$taxon <- factor(ifelse(d$species %in% key, d$species, other_lab), levels = c(key, other_lab))
taxon_shape <- setNames(c(22, 23, 24, 25, 21), c(key, other_lab))
cat("Key taxa:", paste(key, collapse = "; "), "\n")

# Significant pairs whose taxon has no own shape keep a text label.
d$label <- ifelse(d$filled & d$taxon == other_lab, paste0(d$species, " (", d$site, ")"), NA)
cat(sprintf("Holm-significant: %d; text-labelled: %d\n", sum(d$filled), sum(!is.na(d$label))))
print(d[d$filled, c("site", "species", "n_multispecies", "height_shift_m", "ks_D", "sig")], row.names = FALSE)

sig_lab <- c("TRUE" = "Holm-adjusted p < 0.05", "FALSE" = "not significant")
p <- ggplot(d, aes(x = n_multispecies, y = height_shift_m)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
  # fill layer (alpha = significance) under an opaque outline layer, so the
  # faded not-significant points keep a solid black outline
  geom_point(aes(size = ks_D, shape = taxon, fill = site, alpha = filled),
             colour = "transparent", stroke = 0) +
  geom_point(aes(size = ks_D, shape = taxon), fill = NA, colour = "black", stroke = STROKE_W,
             show.legend = FALSE) +
  geom_text_repel(aes(label = label), size = 2.3 + TXT_UP, colour = "grey15", fontface = "italic",
                  min.segment.length = 0, segment.colour = "grey55", segment.size = 0.25,
                  box.padding = 0.4, max.overlaps = Inf, seed = 1, na.rm = TRUE) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_size_area(max_size = 6, name = "KS D", breaks = c(0.01, 0.1, 0.2, 0.4)) +
  scale_shape_manual(values = taxon_shape, name = "Taxon") +
  scale_fill_manual(values = SITE_COL[levels(d$site)], name = "Site") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), labels = sig_lab, name = NULL,
                     breaks = c("TRUE", "FALSE")) +
  guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 2.6, colour = "black", stroke = STROKE_W)),
         shape = guide_legend(order = 2, override.aes = list(size = 2.6, fill = "grey60", colour = "black", stroke = STROKE_W),
                              label.theme = element_text(face = "italic", size = 9)),
         alpha = guide_legend(order = 3, override.aes = list(shape = 21, size = 2.6, fill = "grey30", colour = "black", stroke = STROKE_W)),
         size  = guide_legend(order = 4, override.aes = list(shape = 21, fill = "grey70", colour = "black", stroke = STROKE_W))) +
  labs(x = "Individuals in the joint run (pooled over 5 replicates; log scale)",
       y = "Height shift, joint − isolated (m)") +
  theme_seres() + theme(legend.position = "right")
save_fig("fig_competition", p, width = FULL_W, height = 105)
