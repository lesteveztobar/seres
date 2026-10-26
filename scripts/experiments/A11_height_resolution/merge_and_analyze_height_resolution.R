# merge_and_analyze_height_resolution.R
# Merges the unchanged 0.40 m rows from backup_tier5 with the newly run
# physical-vertical-dispersal (+/- 2.0 m) step files, saves the combined
# height_resolution_<Site>_realistic_75founders.rds files, exports
# output/tables/exports/height_resolution_runs.csv across all 5 sites,
# and prints summary statistics and Kruskal-Wallis / Dunn / Mann-Whitney U tests.

PROCESSED_DIR <- "data/processed"
OUTPUT_DIR    <- "output"

backup_dir <- file.path(PROCESSED_DIR, "backup_tier5")
sites <- c("Maquipucuna", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")
steps_by_site <- list(
  Maquipucuna   = c(0.10, 0.25, 0.40, 0.50, 0.75, 1.00),
  MindoMirador  = c(0.40, 1.00),
  MindoTarabita = c(0.40, 1.00),
  Saloya        = c(0.40, 1.00),
  Yanayacu      = c(0.40, 1.00)
)

all_fin <- list()

for (s in sites) {
  b_file <- file.path(backup_dir, sprintf("height_resolution_%s_realistic_75founders.rds", s))
  b_df <- readRDS(b_file)
  df_040 <- b_df[abs(b_df$height_step - 0.40) < 1e-6, ]
  stopifnot(nrow(df_040) > 0)

  step_dfs <- list(df_040)
  for (st in setdiff(steps_by_site[[s]], 0.40)) {
    st_file <- file.path(PROCESSED_DIR, sprintf("height_resolution_%s_realistic_75founders_step%.2f.rds", s, st))
    if (!file.exists(st_file)) {
      stop("Missing step file: ", st_file)
    }
    st_df <- readRDS(st_file)
    step_dfs[[length(step_dfs) + 1]] <- st_df
  }
  merged <- do.call(rbind, step_dfs)
  merged <- merged[order(merged$height_step, merged$rep, merged$t), ]
  rownames(merged) <- NULL

  out_rds <- file.path(PROCESSED_DIR, sprintf("height_resolution_%s_realistic_75founders.rds", s))
  saveRDS(merged, out_rds)
  cat(sprintf("Saved merged %s (%d rows)\n", out_rds, nrow(merged)))

  fin <- do.call(rbind, lapply(split(merged, list(merged$height_step, merged$rep), drop = TRUE), function(d) {
    d <- d[d$t == max(d$t), ]
    data.frame(
      site = s,
      spacing_m = d$height_step,
      n_tiers = d$n_heights,
      maxDispZ_tiers = max(1L, as.integer(round(2.0 / d$height_step))),
      maxDispZ_m = max(1L, as.integer(round(2.0 / d$height_step))) * d$height_step,
      rep = d$rep,
      t_final = d$t,
      final_total = d$total,
      recruited = d$recruited,
      extinct = d$extinct,
      regime = d$regime,
      stringsAsFactors = FALSE
    )
  }))
  fin <- fin[order(fin$spacing_m, fin$rep), ]
  all_fin[[s]] <- fin
}

hr_runs <- do.call(rbind, all_fin)
rownames(hr_runs) <- NULL
exp_csv <- file.path(OUTPUT_DIR, "tables", "exports", "height_resolution_runs.csv")
write.csv(hr_runs, exp_csv, row.names = FALSE)
cat(sprintf("\nExported %s (%d rows)\n\n", exp_csv, nrow(hr_runs)))

for (s in sites) {
  cat("====================================================================\n")
  cat("SITE:", s, "\n")
  d <- all_fin[[s]]
  agg <- do.call(rbind, lapply(split(d, d$spacing_m), function(x) {
    data.frame(
      spacing_m = x$spacing_m[1],
      n_tiers = x$n_tiers[1],
      maxDispZ_tiers = x$maxDispZ_tiers[1],
      maxDispZ_m = x$maxDispZ_m[1],
      mean_total = round(mean(x$final_total), 1),
      sd_total = round(sd(x$final_total), 1),
      min_total = min(x$final_total),
      max_total = max(x$final_total),
      n_recruited = sum(x$recruited),
      n_extinct = sum(x$extinct),
      n_bounded = sum(x$regime == "bounded")
    )
  }))
  print(agg, row.names = FALSE)

  kw <- kruskal.test(final_total ~ factor(spacing_m), data = d)
  cat(sprintf("Kruskal-Wallis: chi2 = %.4f, df = %d, p = %.4g\n",
              kw$statistic, kw$parameter, kw$p.value))

  if (length(unique(d$spacing_m)) == 2) {
    x04 <- d$final_total[abs(d$spacing_m - 0.40) < 1e-6]
    x10 <- d$final_total[abs(d$spacing_m - 1.00) < 1e-6]
    wt <- wilcox.test(x04, x10, exact = FALSE)
    cat(sprintf("Wilcoxon rank-sum (0.40 vs 1.00 m): W = %.1f, p = %.4g | ratio (1.0/0.4) = %.2f\n",
                wt$statistic, wt$p.value, mean(x10) / mean(x04)))
  } else if (requireNamespace("dunn.test", quietly = TRUE)) {
    cat("Dunn test (Holm-adjusted):\n")
    dt <- dunn.test::dunn.test(d$final_total, g = factor(d$spacing_m), method = "holm", kw = FALSE, table = FALSE)
    res_dt <- data.frame(comparison = dt$comparisons, Z = round(dt$Z, 3), p_unadj = signif(dt$P, 3), p_holm = signif(dt$P.adjusted, 3))
    print(res_dt, row.names = FALSE)
  }
}

