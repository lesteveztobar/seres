# pawn_analysis_v10.R -- full post-processing of the two-sided sensitivity
# campaign (v10: 18,000 main runs across 6 sites x 2 survival arms x 3 seeds,
# plus 540 p_germ descent runs across 9 half-decade levels).
#
# Outputs written to output/sensitivity/:
#   runs_all_v10_main.csv
#   runs_all_v10_descent.csv
#   regime_split_v10.csv
#   pawn_indices_v10.csv
#   src_v10.csv
#   replicate_cv_v10.csv
#   descent_summary_v10.csv
#   arm_comparison_v10.csv
#   sensitivity_report_v10.txt
suppressPackageStartupMessages({
  .libPaths(c("~/R/marvin_libs", .libPaths()))
  library(rpart)
})
source("scripts/02_model/config/paths.R")

PARAM <- c("p_poll", "p_germ", "p_s1", "S", "n_founders")
REGIMES <- c("extinction", "bounded", "runaway")

main_dir <- file.path(OUTPUT_DIR, "sensitivity_runs_v10", "main")
desc_dir <- file.path(OUTPUT_DIR, "sensitivity_runs_v10", "descent")

main_parts <- sort(list.files(main_dir, pattern = "^part_[0-9]+\\.csv$", full.names = TRUE))
desc_parts <- sort(list.files(desc_dir, pattern = "^part_[0-9]+\\.csv$", full.names = TRUE))
stopifnot(length(main_parts) == 600L, length(desc_parts) == 36L)

d_main_all <- do.call(rbind, lapply(main_parts, read.csv, stringsAsFactors = FALSE))
d_desc_all <- do.call(rbind, lapply(desc_parts, read.csv, stringsAsFactors = FALSE))

RES <- file.path(OUTPUT_DIR, "sensitivity")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
outf <- function(stem, ext) file.path(RES, sprintf("%s_v10.%s", stem, ext))

write.csv(d_main_all, outf("runs_all_main", "csv"), row.names = FALSE)
write.csv(d_desc_all, outf("runs_all_descent", "csv"), row.names = FALSE)

# Clean logicals and regimes
clean_df <- function(df) {
  for (b in c("recruited", "extinct", "exceeded_plausible_density")) {
    df[[b]] <- as.integer(as.logical(df[[b]]))
  }
  df$regime <- factor(df$regime, levels = c(REGIMES, "timeout", "error"))
  df
}
d_main <- clean_df(d_main_all)
d_desc <- clean_df(d_desc_all)

pawn_index <- function(x, y, n_int) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2L) return(c(median = NA_real_, n_usable = 0))
  Fu <- ecdf(y)
  br <- quantile(x, probs = seq(0, 1, length.out = n_int + 1), names = FALSE, type = 8)
  br[1] <- -Inf; br[length(br)] <- Inf
  grp <- cut(x, breaks = br, labels = FALSE, include.lowest = TRUE)
  grid <- sort(unique(y))
  ks <- vapply(seq_len(n_int), function(k) {
    yk <- y[grp == k]
    if (length(yk) < 5) return(NA_real_)
    max(abs(ecdf(yk)(grid) - Fu(grid)))
  }, numeric(1))
  c(median = median(ks, na.rm = TRUE), n_usable = sum(!is.na(ks)))
}

sink(outf("sensitivity_report", "txt"), split = TRUE)
cat(sprintf("v10 Main block: %d runs from %d part files (errors/timeouts: %d)\n",
            nrow(d_main), length(main_parts), sum(!d_main$regime %in% REGIMES)))
cat(sprintf("v10 Descent block: %d runs from %d part files (errors/timeouts: %d)\n",
            nrow(d_desc), length(desc_parts), sum(!d_desc$regime %in% REGIMES)))

SITES <- sort(unique(d_main$site))
ARMS  <- c("current", "corrected")

pawn_rows <- list()
src_rows  <- list()
split_rows <- list()
cv_rows   <- list()
tree_rows <- list()

for (arm in ARMS) {
  cat(sprintf("\n\n================================================================================\n"))
  cat(sprintf("                                ARM: %s\n", toupper(arm)))
  cat(sprintf("================================================================================\n"))
  for (site in SITES) {
    df <- d_main[d_main$arm == arm & d_main$site == site & d_main$regime %in% REGIMES, ]
    df$regime <- droplevels(factor(df$regime, levels = REGIMES))
    n_points <- length(unique(df$point_id))
    cat(sprintf("\n### %s (%s arm): %d valid runs, %d parameter points ###\n",
                site, arm, nrow(df), n_points))

    tbl <- table(factor(df$regime, levels = REGIMES))
    print(tbl)
    split_rows[[paste(arm, site, sep = "_")]] <- data.frame(
      arm = arm, site = site, n_runs = nrow(df), n_points = n_points,
      extinction = as.integer(tbl["extinction"]),
      bounded = as.integer(tbl["bounded"]),
      runaway = as.integer(tbl["runaway"]),
      pct_extinction = round(100 * tbl["extinction"] / nrow(df), 2),
      pct_bounded = round(100 * tbl["bounded"] / nrow(df), 2),
      pct_runaway = round(100 * tbl["runaway"] / nrow(df), 2),
      recruited = sum(df$recruited, na.rm = TRUE),
      pct_recruited = round(100 * sum(df$recruited, na.rm = TRUE) / nrow(df), 2),
      median_final_total = median(df$final_total, na.rm = TRUE),
      iqr_final_total = IQR(df$final_total, na.rm = TRUE),
      mean_z_S = mean(df$mean_z_S, na.rm = TRUE),
      mean_z_J = mean(df$mean_z_J, na.rm = TRUE),
      mean_z_A = mean(df$mean_z_A, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    for (r in REGIMES) df[[paste0("regime_", r)]] <- as.integer(df$regime == r)
    outs <- c("final_total", "recruited", "extinct", paste0("regime_", REGIMES))

    for (n_int in c(5, 10, 20)) {
      for (o in outs) {
        if (length(unique(df[[o]][is.finite(df[[o]])])) < 2) next
        raw <- lapply(PARAM, function(p) pawn_index(df[[p]], df[[o]], n_int))
        idx <- setNames(sapply(raw, `[`, "median"), PARAM)
        n_us <- setNames(sapply(raw, `[`, "n_usable"), PARAM)
        rk <- names(sort(idx, decreasing = TRUE))
        if (n_int == 10 && o %in% c("final_total", "recruited", "regime_runaway", "regime_extinction")) {
          cat(sprintf("  PAWN (10 int) %-18s: %s | %s\n", o,
                      paste(sprintf("%s=%.3f", PARAM, idx), collapse = " "),
                      paste(rk, collapse = " > ")))
        }
        pawn_rows[[length(pawn_rows) + 1]] <- data.frame(
          arm = arm, site = site, n_int = n_int, output = o, t(idx),
          rank1 = rk[1], rank2 = rk[2], rank3 = rk[3], min_usable_bins = min(n_us),
          stringsAsFactors = FALSE
        )
      }
    }

    # Standardised regression coefficients (using log10(p_germ) and raw for others, plus raw all)
    for (o in c("final_total", "recruited", "extinct")) {
      dd <- df[is.finite(df[[o]]), ]
      if (length(unique(dd[[o]])) < 2) next
      z <- as.data.frame(scale(dd[, c(o, PARAM)]))
      m <- lm(reformulate(PARAM, response = o), data = z)
      co <- coef(m)[PARAM]
      src_rows[[length(src_rows) + 1]] <- data.frame(
        arm = arm, site = site, output = o, r2 = summary(m)$r.squared, t(co),
        stringsAsFactors = FALSE
      )
      if (o == "final_total") {
        cat(sprintf("  SRC final_total (R2=%.3f): %s\n",
                    summary(m)$r.squared, paste(sprintf("%s=%+.3f", PARAM, co), collapse = " ")))
      }
    }

    # Classification tree
    dt <- droplevels(df[!is.na(df$regime), ])
    if (length(unique(dt$regime)) >= 2) {
      set.seed(1)
      fit <- rpart(reformulate(PARAM, response = "regime"), data = dt, method = "class",
                   control = rpart.control(cp = 0.002, minbucket = 20, xval = 10))
      pruned <- prune(fit, cp = fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"])
      root_acc <- max(table(dt$regime)) / nrow(dt)
      xerr <- pruned$cptable[nrow(pruned$cptable), "xerror"] * (1 - root_acc)
      first_var <- if (nrow(pruned$frame) > 1) as.character(pruned$frame$var[1]) else "none"
      first_split <- if (!is.null(pruned$splits) && nrow(pruned$splits) > 0) pruned$splits[1, "index"] else NA_real_
      cat(sprintf("  Tree: first_var=%s (split=%.6g), baseline=%.3f, CV_acc=%.3f\n",
                  first_var, first_split, root_acc, 1 - xerr))
      tree_rows[[length(tree_rows) + 1]] <- data.frame(
        arm = arm, site = site, first_var = first_var, first_split = first_split,
        baseline_acc = round(100 * root_acc, 2), cv_acc = round(100 * (1 - xerr), 2),
        n_leaves = sum(pruned$frame$var == "<leaf>"),
        stringsAsFactors = FALSE
      )
    }

    # Replicate CV and modal regime per point
    agg <- do.call(rbind, by(df, df$point_id, function(g) {
      if (nrow(g) < 2) return(NULL)
      ft <- g$final_total
      tb <- table(factor(g$regime, levels = REGIMES))
      data.frame(
        arm = arm, site = site, point_id = g$point_id[1],
        idx_num = as.integer(sub(".*_", "", g$point_id[1])),
        n_rep = nrow(g), mean_total = mean(ft),
        cv_total = if (mean(ft) > 0) sd(ft) / mean(ft) else NA_real_,
        regime_mode = names(which.max(tb)),
        regime_unstable = length(unique(g$regime)) > 1,
        stringsAsFactors = FALSE
      )
    }))
    cv_rows[[paste(arm, site, sep = "_")]] <- agg
    cat(sprintf("  Replicate median CV: %.3f (90th pct: %.3f) | Unstable regime points: %d / %d (%.1f%%)\n",
                median(agg$cv_total, na.rm = TRUE),
                quantile(agg$cv_total, 0.9, na.rm = TRUE),
                sum(agg$regime_unstable), nrow(agg), 100 * mean(agg$regime_unstable)))
  }
}

split_df <- do.call(rbind, split_rows)
pawn_df  <- do.call(rbind, pawn_rows)
src_df   <- do.call(rbind, src_rows)
cv_df    <- do.call(rbind, cv_rows)
tree_df  <- do.call(rbind, tree_rows)

write.csv(split_df, outf("regime_split", "csv"), row.names = FALSE)
write.csv(pawn_df,  outf("pawn_indices", "csv"), row.names = FALSE)
write.csv(src_df,   outf("src", "csv"), row.names = FALSE)
write.csv(cv_df,    outf("replicate_cv", "csv"), row.names = FALSE)
write.csv(tree_df,  outf("trees", "csv"), row.names = FALSE)

# Point-by-point nesting check across sites per arm
cat("\n\n=== POINT-BY-POINT REGIME NESTING AND SPEARMAN CORRELATIONS ===\n")
for (arm in ARMS) {
  cat(sprintf("\n--- Arm: %s ---\n", arm))
  sub_cv <- cv_df[cv_df$arm == arm, ]
  wide_reg <- reshape(sub_cv[, c("idx_num", "site", "regime_mode")],
                      idvar = "idx_num", timevar = "site", direction = "wide")
  wide_tot <- reshape(sub_cv[, c("idx_num", "site", "mean_total")],
                      idvar = "idx_num", timevar = "site", direction = "wide")
  names(wide_tot) <- sub("^mean_total\\.", "", names(wide_tot))
  cor_mat <- cor(wide_tot[, SITES], method = "spearman", use = "pairwise.complete.obs")
  cat("Spearman rank correlation of per-point mean final abundance:\n")
  print(round(cor_mat, 3))
}

# Descent block summary
cat("\n\n=== DESCENT BLOCK SUMMARY (p_germ from 1e-2 to 1e-6) ===\n")
desc_summary <- do.call(rbind, by(d_desc, list(d_desc$arm, d_desc$site, d_desc$p_germ), function(g) {
  if (nrow(g) == 0) return(NULL)
  tbl <- table(factor(g$regime, levels = REGIMES))
  data.frame(
    arm = g$arm[1], site = g$site[1], p_germ = g$p_germ[1], n_runs = nrow(g),
    n_extinction = as.integer(tbl["extinction"]),
    n_bounded = as.integer(tbl["bounded"]),
    n_runaway = as.integer(tbl["runaway"]),
    n_recruited = sum(g$recruited, na.rm = TRUE),
    mean_final_total = round(mean(g$final_total, na.rm = TRUE), 1),
    sd_final_total = round(sd(g$final_total, na.rm = TRUE), 1),
    stringsAsFactors = FALSE
  )
}))
desc_summary <- desc_summary[order(desc_summary$arm, desc_summary$site, -desc_summary$p_germ), ]
write.csv(desc_summary, outf("descent_summary", "csv"), row.names = FALSE)
print(desc_summary, row.names = FALSE)

cat("\n\n=== SUMMARY TABLE BY ARM AND SITE ===\n")
print(split_df, row.names = FALSE)
cat("\n=== CLASSIFICATION TREES ===\n")
print(tree_df, row.names = FALSE)
sink()
