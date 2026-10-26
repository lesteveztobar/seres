# fig_heldout.R -- held-out validation figures (sec:res-vertposition,
# sec:res-threeway).
#   fig_heldout_heights.{png,pdf}: modelled height distribution at t = 50 per
#     replicate (abundance-weighted boxplots, all stages) against the heights
#     of the held-out individuals. Maquipucuna and Mashpi meet the 4-record
#     minimum and are the tested sites; the others are shown for reference.
#   fig_heldout_scores.{png,pdf}: each held-out individual's climate-niche
#     suitability under the three predictors, joined by a line; bars are
#     medians; paired Wilcoxon tests recomputed here and printed.
# Inputs: output/tables/exports/heldout_{modelled,observed}_heights.csv,
#   output/tables/validation/held_out_validation_scores_realistic.csv
#   sbatch scripts/04_results_figures/run_fig.sh fig_heldout.R
source("scripts/04_results_figures/fig_common.R")

TESTED <- c("Maquipucuna", "Mashpi")
mod <- read.csv(file.path(EXP_DIR, "heldout_modelled_heights.csv"))
obs <- read.csv(file.path(EXP_DIR, "heldout_observed_heights.csv"))
sites <- intersect(SITE_LEVELS, union(mod$site, obs$site))
nh <- as.integer(table(obs$site)[sites])
strip <- setNames(sprintf("%s%s\n%d held out\n(%s)", sites, ifelse(sites == OUTGROUP, " (outgroup)", ""), nh,
                          ifelse(sites %in% TESTED, "tested", "not tested")), sites)
mod$site_f <- factor(strip[mod$site], levels = strip); obs$site_f <- factor(strip[obs$site], levels = strip)
mod$x <- factor(paste("Rep", mod$rep), levels = c(paste("Rep", 1:5), "Held out"))
obs$x <- factor("Held out", levels = levels(mod$x))
mod$tested <- mod$site %in% TESTED
# one row per individual (weighted boxplots would need quantreg, absent on the cluster)
mod <- mod[rep(seq_len(nrow(mod)), mod$n_individuals), ]

ph <- ggplot() +
  geom_boxplot(data = mod, aes(x = x, y = height_m, fill = site, alpha = tested),
               colour = "grey25", linewidth = 0.3, outlier.shape = NA, width = 0.6) +
  geom_boxplot(data = obs, aes(x = x, y = observed_height), fill = "white", colour = "grey25",
               linewidth = 0.3, outlier.shape = NA, width = 0.6) +
  geom_point(data = obs, aes(x = x, y = observed_height, fill = site), shape = 21, colour = "black",
             stroke = STROKE_W, size = 1.8, position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  scale_fill_manual(values = SITE_COL, guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.35), guide = "none") +
  facet_wrap(~site_f, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = "Height above ground (m)") +
  theme_seres() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        strip.text = element_text(size = 8.5, lineheight = 1.05))
save_fig("fig_heldout_heights", ph, height = 80)

# ── three predictors ─────────────────────────────────────────────────────────
s <- read.csv("output/tables/validation/held_out_validation_scores_realistic.csv")
stopifnot(all(tapply(s$score_own_position, paste(s$site, s$species, s$observed_height), function(v) diff(range(v)) < 1e-9)))
s <- s[s$rep == 1, ]
s$id <- seq_len(nrow(s))
w_site <- wilcox.test(s$score_own_position, s$score_site_mean, paired = TRUE)
w_hgt  <- wilcox.test(s$score_own_position, s$score_height_only, paired = TRUE)
cat(sprintf("n = %d | own > site-shared in %d (V = %.1f, p = %.3f) | own > height-only in %d (V = %.1f, p = %.3f)\n",
            nrow(s), sum(s$score_own_position > s$score_site_mean), w_site$statistic, w_site$p.value,
            sum(s$score_own_position > s$score_height_only), w_hgt$statistic, w_hgt$p.value))
PRED <- c(score_site_mean = "Site-shared\nclimate", score_height_only = "Height only\n(raster mean)",
          score_own_position = "Own voxel")
long <- do.call(rbind, lapply(names(PRED), function(k)
  data.frame(id = s$id, site = s$site, predictor = PRED[[k]], score = s[[k]])))
long$predictor <- factor(long$predictor, levels = PRED)
long$site <- factor(long$site, levels = intersect(SITE_LEVELS, long$site))
long$outgroup <- long$site == OUTGROUP
med <- aggregate(score ~ predictor, long, median)
lab <- sprintf("Own voxel higher than:\n  site-shared in %d of %d individuals\n  height only in %d of %d individuals",
               sum(s$score_own_position > s$score_site_mean), nrow(s),
               sum(s$score_own_position > s$score_height_only), nrow(s))
site_leg <- setNames(ifelse(levels(long$site) == OUTGROUP, paste(OUTGROUP, "(outgroup)"), levels(long$site)), levels(long$site))
ps <- ggplot(long, aes(x = predictor, y = score)) +
  geom_line(aes(group = id, colour = site, linetype = outgroup), linewidth = 0.4, alpha = 0.8) +
  geom_point(aes(fill = site), shape = 21, colour = "black", stroke = STROKE_W, size = 1.9) +
  geom_crossbar(data = med, aes(y = score, ymin = score, ymax = score), width = 0.45, linewidth = 0.5,
                colour = "black", middle.linewidth = 0.6) +
  annotate("label", x = 0.55, y = 98, label = lab, hjust = 0, vjust = 1, size = 2.1 + TXT_UP, colour = "grey20",
           fill = alpha("white", 0.85), linewidth = 0.2, label.padding = unit(1.8, "mm"), lineheight = 1.1) +
  scale_y_sqrt(limits = c(0, 100), breaks = c(0, 1, 5, 10, 25, 50, 100)) +
  scale_colour_manual(values = SITE_COL, labels = site_leg, name = "Site") +
  scale_fill_manual(values = SITE_COL, labels = site_leg, name = "Site") +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"), guide = "none") +
  labs(x = NULL, y = "Climate-niche suitability of held-out individual\n(0–100, square-root scale)") +
  theme_seres() + theme(legend.position = "right")
save_fig("fig_heldout_scores", ps, width = 150, height = 90)
print(med)
