# fig_competition_versions.R -- two layouts of the competition/isolation
# result (2026-10-03), to compare against fig_competition.R's single panel.
#   fig_competition_by_site.{png,pdf}:  one panel per site; shape = taxon for
#     the six taxa modelled at 2+ sites, circle + italic label for the
#     single-site taxa.
#   fig_competition_by_taxon.{png,pdf}: one panel per taxon unit; shape = site.
# In both: y = mean height shift (joint - isolated), x = individuals in the
# joint run (log), point area = KS D, filled = Holm-adjusted p < 0.05,
# white = not significant. Yanayacu's single pair (Ida sp., no comparison) is dropped.
#   sbatch scripts/04_results_figures/run_fig.sh fig_competition_versions.R
source("scripts/04_results_figures/fig_common.R")
suppressPackageStartupMessages(library(ggrepel))

d <- read.csv("output/tables/competition/competition_analysis.csv")
d <- d[!(d$site == "Yanayacu" & d$species == "Ida sp."), ]
d$sig <- d$sig != "ns"
d$site <- factor(d$site, levels = intersect(SITE_LEVELS, d$site))
n_at <- table(d$species)
multi <- names(sort(n_at[n_at >= 2], decreasing = TRUE))
cat(sprintf("%d pairs, %d taxa, %d at 2+ sites: %s\n", nrow(d), length(n_at), length(multi), paste(multi, collapse = "; ")))
sig_scale <- function() scale_alpha_identity()
size_scale <- scale_size_area(max_size = 5, name = "KS D", breaks = c(0.01, 0.1, 0.2, 0.3))
xs <- scale_x_log10(labels = scales::label_comma(), expand = expansion(mult = 0.12))
ys <- scale_y_continuous(expand = expansion(mult = 0.12))

# ── version 1: one panel per site ────────────────────────────────────────────
fill_shapes <- c(23, 24, 22, 25)                       # fillable shapes for taxa that can be significant
any_sig <- tapply(d$sig, d$species, any)
ord_multi <- c(multi[any_sig[multi]], multi[!any_sig[multi]])
shape_map <- setNames(c(fill_shapes, 4, 3)[seq_along(ord_multi)], ord_multi)
stopifnot(!any(any_sig[names(shape_map)[!shape_map %in% 21:25]]))  # unfillable shapes only for never-significant taxa
single_lab <- "Single-site taxon (labelled)"
d$taxon_key <- factor(ifelse(d$species %in% multi, d$species, single_lab), levels = c(ord_multi, single_lab))
shape_all <- c(shape_map, setNames(21, single_lab))
d$fillcol <- ifelse(d$sig, SITE_COL[as.character(d$site)], "white")
d$lab <- ifelse(d$species %in% multi, NA, taxon_expr(d$species, abbreviate = TRUE))
leg_lab <- function(x) parse(text = ifelse(x == single_lab, sprintf("'%s'", x), taxon_expr(x)))

p1 <- ggplot(d, aes(x = n_multispecies, y = height_shift_m)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_point(aes(shape = taxon_key, size = ks_D, fill = I(fillcol)), colour = "black", stroke = STROKE_W) +
  geom_text_repel(aes(label = lab), parse = TRUE, size = 2.1 + TXT_UP, colour = "grey25", na.rm = TRUE,
                  min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60",
                  box.padding = 0.3, max.overlaps = Inf, seed = 1) +
  xs + ys + size_scale +
  scale_shape_manual(values = shape_all, labels = leg_lab, name = "Taxon") +
  facet_wrap(~site, ncol = 3) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 2.6, fill = "grey60")),
         size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey80"))) +
  labs(x = "Individuals in the joint run (pooled over 5 replicates; log scale)",
       y = "Height shift, joint − isolated (m)",
       caption = "Filled: Holm-adjusted p < 0.05 (site colour); white: not significant.") +
  theme_seres() +
  theme(strip.text = element_text(face = "bold", size = 9.5), legend.position = c(0.84, 0.25),
        legend.box = "horizontal", legend.text = element_text(size = 8.5))
save_fig("fig_competition_by_site", p1, height = 120)

# ── version 2: one panel per taxon unit ──────────────────────────────────────
site_shape <- setNames(c(21, 22, 23, 24, 25), levels(d$site))
sp_ord <- c(multi, sort(setdiff(names(n_at), multi)))
d$taxon_f <- factor(taxon_expr(d$species, abbreviate = TRUE), levels = taxon_expr(sp_ord, abbreviate = TRUE))
p2 <- ggplot(d, aes(x = n_multispecies, y = height_shift_m)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_point(aes(shape = site, size = ks_D, fill = I(fillcol)), colour = "black", stroke = STROKE_W) +
  xs + ys + size_scale +
  scale_shape_manual(values = site_shape, name = "Site") +
  facet_wrap(~taxon_f, ncol = 6, labeller = label_parsed) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 2.6, fill = "grey60")),
         size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey80"))) +
  labs(x = "Individuals in the joint run (pooled over 5 replicates; log scale)",
       y = "Height shift, joint − isolated (m)",
       caption = "Filled: Holm-adjusted p < 0.05 (site colour); white: not significant.") +
  theme_seres() +
  theme(strip.text = element_text(size = 8.5), legend.position = "bottom",
        axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1), panel.spacing = unit(2, "mm"))
save_fig("fig_competition_by_taxon", p2, height = 150)
print(d[d$sig, c("site", "species", "height_shift_m", "ks_D")], row.names = FALSE)
