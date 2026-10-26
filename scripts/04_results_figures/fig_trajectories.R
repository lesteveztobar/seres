# fig_trajectories.R -- Results figure (realistictraj): population
# trajectories per taxon unit. Colour = site (SITE_COL), shape + linetype =
# life stage (S/J/A), columns = model configuration, rows = taxon unit. Lines
# are the mean of the replicates; markers every 5 years. y is log10(1 + n) so
# seedling and juvenile counts stay visible next to adults.
# Replaces the six base-R abundance_<site>_realistic_h0.40_rep1.png panels.
# Configurations are drawn only if their colonization_<site>_<tag>_h0.40.rds
# exists; best_case_kcal = best_case.rds rerun after the 2026-09-29 K
# calibration (the older colonization_*_best_case_h0.40.rds predate it).
# Per-taxon totals are cached to output/tables/exports/trajectories_by_taxon_<tag>.csv.
# Outputs: fig_trajectories_expanding / _steady.{png,pdf} (taxa at 2+ sites, split 2026-10-03),
#          fig_trajectories_single.{png,pdf} (taxa at one site, appendix)
#   Rscript scripts/04_results_figures/fig_trajectories.R
source("scripts/04_results_figures/fig_common.R")

PROCESSED_DIR <- "data/processed"
SITES <- c("Mashpi", "Saloya", "Maquipucuna", "MindoTarabita", "MindoMirador", "Yanayacu")
CONFIGS <- c(realistic = "Realistic", best_case_kcal = "Best case")
STAGE_LAB <- c(S = "Seedling", J = "Juvenile", A = "Adult")

extract_tag <- function(tag) {
  cache <- file.path(EXP_DIR, sprintf("trajectories_by_taxon_%s.csv", tag))
  files <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s_h0.40.rds", SITES, tag))
  have <- file.exists(files)
  if (!any(have)) return(NULL)
  if (file.exists(cache) && file.mtime(cache) >= max(file.mtime(files[have]))) {
    cat("Using cache", cache, "\n"); return(read.csv(cache))
  }
  out <- do.call(rbind, lapply(which(have), function(i) {
    res <- readRDS(files[i])
    runs <- Filter(Negate(is.null), if (!is.null(res$runs)) res$runs else list(res))
    rm(res)
    d <- do.call(rbind, lapply(seq_along(runs), function(r) {
      run <- runs[[r]]
      do.call(rbind, lapply(c("S", "J", "A"), function(st) {
        m <- apply(run[[paste0("abundance", st)]], c(4, 5), sum)   # t x species
        data.frame(site = SITES[i], rep = r, stage = st,
                   taxon = rep(run$species_ids, each = nrow(m)),
                   t = rep(seq_len(nrow(m)), ncol(m)), n = as.vector(m))
      }))
    }))
    rm(runs); gc()
    cat(sprintf("%s / %s: %d rows, %d reps\n", SITES[i], tag, nrow(d), length(unique(d$rep))))
    d
  }))
  out$config <- tag
  write.csv(out, cache, row.names = FALSE); cat("Saved:", cache, "\n")
  out
}

d <- do.call(rbind, lapply(names(CONFIGS), extract_tag))
cat("Configurations present:", paste(unique(d$config), collapse = ", "), "\n")
m <- aggregate(n ~ config + site + taxon + stage + t, data = d, FUN = mean)
m$config <- factor(CONFIGS[m$config], levels = CONFIGS[names(CONFIGS) %in% m$config])
m$site <- factor(m$site, levels = intersect(SITE_LEVELS, m$site))
m$stage <- factor(STAGE_LAB[m$stage], levels = STAGE_LAB)

ref <- m[m$config == levels(m$config)[1], ]
n_sites <- tapply(ref$site, ref$taxon, function(s) length(unique(s)))
lab <- function(x) ifelse(grepl("^[A-Z][a-z]+ sp", x), x, sub("^([A-Z])[a-z]+\\s+", "\\1. ", x))

# 2026-10-03: the main-text figures show the realistic configuration only.
# best_case_kcal runs exceed the runaway ceiling during spin-up, so their
# trajectories are flat carried-forward states and would need their own
# explanation; they stay in the per-taxon cache for reference.
make_plot <- function(keep, log_y = TRUE, ncol = 1) {
  dd <- m[m$taxon %in% keep & m$config == "Realistic", ]
  dd$taxon <- factor(taxon_expr(dd$taxon, abbreviate = TRUE), levels = unique(taxon_expr(sort(keep), abbreviate = TRUE)))
  pts <- dd[dd$t %% 5 == 0, ]
  yscale <- if (log_y) {
    scale_y_continuous(trans = scales::log1p_trans(), breaks = c(0, 1, 10, 100, 1000, 10000),
                       labels = scales::label_comma())
  } else scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA))
  facet <- if (nlevels(droplevels(dd$config)) > 1) {
    facet_grid(taxon ~ config, labeller = labeller(taxon = label_parsed))
  } else facet_wrap(~taxon, ncol = ncol, labeller = label_parsed)
  ggplot(dd, aes(x = t, y = n, colour = site, group = interaction(site, stage))) +
    geom_line(aes(linetype = stage), linewidth = 0.4) +
    geom_point(data = pts, aes(shape = stage, fill = site), colour = "black", size = 1.4, stroke = 0.3) +
    yscale +
    scale_colour_manual(values = SITE_COL[levels(droplevels(dd$site))], name = "Site", drop = TRUE) +
    scale_fill_manual(values = SITE_COL[levels(droplevels(dd$site))], guide = "none", drop = TRUE) +
    scale_shape_manual(values = c(Seedling = 21, Juvenile = 24, Adult = 22), name = "Stage") +
    scale_linetype_manual(values = c(Seedling = "dotted", Juvenile = "dashed", Adult = "solid"), name = "Stage") +
    facet +
    labs(x = "Year", y = sprintf("Individuals (mean of replicates%s)", if (log_y) "; log scale" else "")) +
    theme_seres() +
    theme(strip.text = element_text(size = 9.5, hjust = 0), legend.position = "bottom",
          legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))
}

multi <- names(n_sites)[n_sites >= 2]; single <- names(n_sites)[n_sites < 2]
cat("Taxa at 2+ sites:", paste(multi, collapse = "; "), "\n")
# 2026-10-03: the multi-site figure is split by whether the taxon unit
# expands. "Expanding" = mean total abundance (all stages) at the final year
# reaches EXPAND_N at one or more sites; these keep the log axis. The others
# stay near their founder numbers and read better on a linear axis.
EXPAND_N <- 100
fin <- aggregate(n ~ taxon + site, data = ref[ref$t == max(ref$t), ], FUN = sum)
peak_final <- tapply(fin$n, fin$taxon, max)
expanding <- multi[peak_final[multi] >= EXPAND_N]; steady <- setdiff(multi, expanding)
cat("Expanding:", paste(expanding, collapse = "; "), "| steady:", paste(steady, collapse = "; "), "\n")
print(round(peak_final[multi]))
w <- 150
save_fig("fig_trajectories_expanding", make_plot(expanding, log_y = TRUE), width = w, height = 40 + 38 * length(expanding))
save_fig("fig_trajectories_steady", make_plot(steady, log_y = FALSE), width = w, height = 40 + 34 * length(steady))
save_fig("fig_trajectories_single", make_plot(single, log_y = TRUE, ncol = 4), width = FULL_W,
         height = 45 + 34 * ceiling(length(single) / 4))
