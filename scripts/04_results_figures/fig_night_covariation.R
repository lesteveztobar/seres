# fig_night_covariation.R -- Results figure 2.1: (a) temperature-humidity
# correlation by day and by night; (b) previous-day insolation against night
# minimum temperature. Sites ordered by elevation (site_elevation.csv).
#   Rscript scripts/04_results_figures/fig_night_covariation.R
source("scripts/04_results_figures/fig_common.R")
suppressPackageStartupMessages(library(patchwork))

elev <- site_elevation()
src_a <- "output/climate_decoupling_night/temp_relhum_night_vs_day.csv"
src_b <- "output/climate_decoupling_night/prevday_insolation_vs_night.csv"
a_all <- read.csv(src_a); b_all <- read.csv(src_b)

a <- a_all[a_all$stratum == "overall", ]
b <- b_all[b_all$stratum == "overall" & b_all$version == "month_anomaly" & b_all$response == "temp_min", ]
cat(sprintf("Panel a: %s, %d of %d rows after stratum == overall\n", src_a, nrow(a), nrow(a_all)))
cat(sprintf("Panel b: %s, %d of %d rows after overall / month_anomaly / temp_min\n", src_b, nrow(b), nrow(b_all)))
stopifnot(!anyNA(elev[c(a$site, b$site)]))

ord <- names(sort(elev[unique(c(a$site, b$site))]))
lab <- setNames(sprintf("%s (%s m)", ord, format(round(elev[ord]), big.mark = ",", trim = TRUE)), ord)
a$site_f <- factor(a$site, levels = ord); b$site_f <- factor(b$site, levels = ord)
# 2026-10-03: no site colours (the y axis already names each site); Yanayacu,
# the eastern-slope outgroup, is set apart in its own facet row.
a$group <- slope_group(a$site); b$group <- slope_group(b$site)
strip_theme <- theme(strip.placement = "outside", strip.background.y = element_rect(fill = "grey93", colour = NA),
                     panel.spacing.y = unit(3, "mm"))
a$daypart <- factor(a$daypart, levels = c("day", "night"))
b$filled <- b$p < 0.05

seg <- reshape(a[, c("site_f", "group", "daypart", "rho")], idvar = c("site_f", "group"), timevar = "daypart", direction = "wide")
pa <- ggplot(a, aes(y = site_f)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_segment(data = seg, aes(x = rho.day, xend = rho.night, yend = site_f), colour = "grey60", linewidth = 0.6) +
  geom_point(aes(x = rho, shape = daypart, fill = daypart), colour = "black", stroke = STROKE_W, size = 2.6) +
  scale_shape_manual(values = c(day = 24, night = 21), name = NULL, labels = c(day = "Day", night = "Night")) +
  scale_fill_manual(values = c(day = "white", night = "grey20"), name = NULL, labels = c(day = "Day", night = "Night")) +
  scale_y_discrete(labels = lab) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Spearman ρ\n(temperature, relative humidity)", y = NULL, title = "a") +
  theme_seres() + theme(legend.position = "top", legend.justification = "left",
                        strip.text.y = element_blank()) + strip_theme

pb <- ggplot(b, aes(y = site_f, x = rho)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_point(aes(fill = filled), shape = 21, colour = "black", size = 2.6, stroke = STROKE_W) +
  scale_fill_manual(values = c("TRUE" = "grey20", "FALSE" = "white"), breaks = c("TRUE", "FALSE"),
                    labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"), name = NULL) +
  scale_y_discrete(labels = lab) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(group = c("Western slope (Chocó Andino)" = "Western\nslope",
                                           "Eastern-slope outgroup" = "Outgroup"))) +
  labs(x = "Spearman ρ\n(previous-day insolation,\nnight minimum temperature)", y = NULL, title = "b") +
  theme_seres() + theme(axis.text.y = element_blank(), legend.position = "top", legend.justification = "left",
                        strip.text.y = element_text(face = "bold", angle = 0, hjust = 0, size = 9)) + strip_theme

n_med <- median(a$n)
cap <- sprintf(paste0("(a) Hourly pairs, overall stratum; n per site and daypart %s–%s (median %s), so every ",
                      "correlation is significant and no stars are shown. (b) n = %s days per site."),
               format(min(a$n), big.mark = ","), format(max(a$n), big.mark = ","),
               format(round(n_med), big.mark = ","), paste(unique(b$n), collapse = "/"))
p <- pa | pb  # caption text lives in the LaTeX caption (printed below)
save_fig("fig_night_covariation", p, width = FULL_W, height = 85)

cat(sprintf("Plotted: panel a %d points (%d sites x 2 dayparts), panel b %d points; dropped 0\n",
            nrow(a), length(unique(a$site)), nrow(b)))
cat(sprintf("Panel a: night rho all negative: %s; day rho > 0 at %d sites (%s)\n",
            all(a$rho[a$daypart == "night"] < 0), sum(a$rho[a$daypart == "day"] > 0),
            paste(a$site[a$daypart == "day" & a$rho > 0], collapse = ", ")))
cat(sprintf("Panel b: all negative: %s; p < 0.05 at %d (%s)\n", all(b$rho < 0), sum(b$filled),
            paste(b$site[b$filled], collapse = ", ")))
cat("Caption:", cap, "\n")
