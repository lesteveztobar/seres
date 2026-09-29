# pawn_analysis.R -- generic analysis of a sensitivity design run through
# sensitivity_run.R. Works for any design tag, any set of sites and any set
# of varied parameters.
#
# Usage: Rscript pawn_analysis.R <tag> [ranges_csv]
#   tag         design tag; part files are read from every directory
#               output/sensitivity_runs_<tag>*
#   ranges_csv  the ranges file the design was built from (names the varied
#               parameters); default data/params/sensitivity_ranges_<tag>.csv
# Per site: regime split, PAWN indices (5/10/20 conditioning intervals),
# standardised regression coefficients, regime classification tree, and
# within-parameter-set replicate CV (the stochastic noise floor). Runs with
# a recorded t=50 snapshot that differs from the final state (designs run
# past 50 timesteps) also get PAWN at the snapshot, for a time-horizon check.
# Analyses whatever runs have completed and reports n per site.
suppressPackageStartupMessages({
  .libPaths(c("~/R/marvin_libs", .libPaths()))
  library(rpart)
})
source("scripts/02_model/config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: pawn_analysis.R <tag> [ranges_csv]")
tag <- args[1]
ranges_csv <- if (length(args) >= 2) args[2] else file.path(PARAMS_DIR, sprintf("sensitivity_ranges_%s.csv", tag))
PARAM <- read.csv(ranges_csv, stringsAsFactors = FALSE)$param

dirs <- Sys.glob(file.path(OUTPUT_DIR, sprintf("sensitivity_runs_%s*", tag)))
dirs <- dirs[dir.exists(dirs)]
parts <- unlist(lapply(dirs, list.files, pattern = "^part_.*\\.csv$", full.names = TRUE))
if (!length(parts)) stop("No part files found for tag ", tag)
d_all <- do.call(rbind, lapply(parts, read.csv, stringsAsFactors = FALSE))
n_err <- sum(!(is.na(d_all$error) | d_all$error == ""))
d <- d_all[is.na(d_all$error) | d_all$error == "", ]

RES <- file.path(OUTPUT_DIR, "sensitivity")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
outf <- function(stem, ext) file.path(RES, sprintf("%s_%s.%s", stem, tag, ext))
write.csv(d_all, outf("runs_all", "csv"), row.names = FALSE)

OUT_CONT <- c("final_total", "realised_mean_height")
OUT_BIN <- c("recruited", "extinct", "exceeded_plausible_density")
REGIMES <- c("extinction", "bounded", "runaway")
d$regime <- factor(d$regime, levels = REGIMES)
for (b in OUT_BIN) d[[b]] <- as.integer(as.logical(d[[b]]))
has_snapshot <- "t50_total" %in% names(d) && any(is.finite(d$t50_total)) &&
  any(d$t50_total != d$final_total, na.rm = TRUE)

sink(outf("sensitivity_report", "txt"), split = TRUE)
cat(sprintf("Design '%s': %d runs read from %d part files (%d errored runs dropped)\n", tag, nrow(d_all), length(parts), n_err))
cat(sprintf("Varied parameters: %s\n", paste(PARAM, collapse = ", ")))
SITES <- sort(unique(d$site))

pawn_index <- function(x, y, n_int) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
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

pawn_rows <- list(); src_rows <- list(); split_rows <- list(); cv_rows <- list()
for (site in SITES) {
  df <- d[d$site == site, ]
  n_points <- length(unique(df$point_id))
  cat(sprintf("\n\n################################ %s: %d runs, %d parameter sets ################################\n",
              site, nrow(df), n_points))

  cat("Regime split:\n"); print(table(df$regime, useNA = "ifany"))
  split_rows[[site]] <- data.frame(site = site, n_runs = nrow(df), n_points = n_points,
    t(setNames(vapply(REGIMES, function(r) sum(df$regime == r, na.rm = TRUE), integer(1)), REGIMES)),
    recruited = sum(df$recruited, na.rm = TRUE), median_final_total = median(df$final_total, na.rm = TRUE))

  for (r in REGIMES) df[[paste0("regime_", r)]] <- as.integer(df$regime == r)
  outs <- c(OUT_CONT, OUT_BIN, paste0("regime_", REGIMES), if (has_snapshot) c("t50_total", "t50_realised_mean_height"))
  for (n_int in c(5, 10, 20)) {
    cat(sprintf("\n--- PAWN, %d conditioning intervals ---\n", n_int))
    for (o in outs) {
      if (length(unique(df[[o]][is.finite(df[[o]])])) < 2) {
        cat(sprintf("%-27s : constant in this sample -- no index\n", o)); next
      }
      raw <- lapply(PARAM, function(p) pawn_index(df[[p]], df[[o]], n_int))
      idx <- setNames(sapply(raw, `[`, "median"), PARAM)
      n_us <- setNames(sapply(raw, `[`, "n_usable"), PARAM)
      rk <- names(sort(idx, decreasing = TRUE))
      cat(sprintf("%-27s : %s | rank: %s%s\n", o, paste(sprintf("%s=%.3f", PARAM, idx), collapse = "  "),
                  paste(rk, collapse = " > "), if (any(n_us < n_int)) "   [THIN BINS]" else ""))
      pawn_rows[[length(pawn_rows) + 1]] <- data.frame(site = site, n_int = n_int, output = o, t(idx),
        rank1 = rk[1], rank2 = rk[2], min_usable_bins = min(n_us))
    }
  }

  cat("\n--- Standardised regression coefficients ---\n")
  for (o in c(OUT_CONT, OUT_BIN)) {
    dd <- df[is.finite(df[[o]]), ]
    if (length(unique(dd[[o]])) < 2) { cat(sprintf("%-27s : constant -- skipped\n", o)); next }
    z <- as.data.frame(scale(dd[, c(o, PARAM)]))
    m <- lm(reformulate(PARAM, response = o), data = z)
    co <- coef(m)[PARAM]
    cat(sprintf("%-27s (R2=%.3f) : %s\n", o, summary(m)$r.squared, paste(sprintf("%s=%+.3f", PARAM, co), collapse = "  ")))
    src_rows[[length(src_rows) + 1]] <- data.frame(site = site, output = o, r2 = summary(m)$r.squared, t(co))
  }

  cat("\n--- Regime classification tree ---\n")
  dt <- droplevels(df[!is.na(df$regime), ])
  if (length(unique(dt$regime)) < 2) {
    cat("regime is constant at this site -- no tree.\n")
  } else {
    set.seed(1)
    fit <- rpart(reformulate(PARAM, response = "regime"), data = dt, method = "class",
                 control = rpart.control(cp = 0.002, minbucket = 20, xval = 10))
    pruned <- prune(fit, cp = fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"])
    print(pruned)
    root_acc <- max(table(dt$regime)) / nrow(dt)
    xerr <- pruned$cptable[nrow(pruned$cptable), "xerror"] * (1 - root_acc)
    cat(sprintf("majority-class accuracy %.3f | pruned-tree 10-fold CV accuracy %.3f\n", root_acc, 1 - xerr))
    saveRDS(pruned, file.path(RES, sprintf("regime_tree_%s_%s.rds", tag, site)))
    if (nrow(pruned$frame) > 1) {
      png(file.path(RES, sprintf("regime_tree_%s_%s.png", tag, site)), width = 1100, height = 800, res = 110, type = "cairo")
      par(mar = c(1, 1, 2, 1)); plot(pruned, uniform = TRUE, margin = 0.08)
      text(pruned, use.n = TRUE, cex = 0.75); title(sprintf("Regime tree -- %s (%s)", site, tag))
      dev.off()
    }
  }

  agg <- do.call(rbind, by(df, df$point_id, function(g) {
    if (nrow(g) < 2) return(NULL)
    ft <- g$final_total
    data.frame(site = site, point_id = g$point_id[1], n_rep = nrow(g), mean_total = mean(ft),
               cv_total = if (mean(ft) > 0) sd(ft) / mean(ft) else NA_real_,
               regime_mode = names(which.max(table(g$regime))),
               regime_unstable = length(unique(g$regime)) > 1)
  }))
  cv_rows[[site]] <- agg
  cat("\n--- Within-parameter-set replicate CV of final total (noise floor) ---\n")
  for (rg in setdiff(REGIMES, "extinction")) {
    cv <- agg$cv_total[agg$regime_mode == rg]
    if (!length(cv)) next
    cat(sprintf("%-8s (n=%d parameter sets): ", rg, sum(!is.na(cv))))
    print(round(quantile(cv, c(.1, .25, .5, .75, .9, .95, 1), na.rm = TRUE), 3))
  }
  cat(sprintf("parameter sets whose replicates disagree on regime: %d / %d\n", sum(agg$regime_unstable), nrow(agg)))
}

pawn_df <- do.call(rbind, pawn_rows)
write.csv(pawn_df, outf("pawn_indices", "csv"), row.names = FALSE)
write.csv(do.call(rbind, src_rows), outf("src", "csv"), row.names = FALSE)
write.csv(do.call(rbind, split_rows), outf("regime_split", "csv"), row.names = FALSE)
write.csv(do.call(rbind, cv_rows), outf("replicate_cv", "csv"), row.names = FALSE)

cat("\n\n################################ CROSS-SITE SUMMARY ################################\n")
print(do.call(rbind, split_rows), row.names = FALSE)
cat("\nPAWN ranking of final_total, 10 intervals, by site:\n")
sub <- pawn_df[pawn_df$n_int == 10 & pawn_df$output == "final_total", ]
for (i in seq_len(nrow(sub))) {
  idx <- sort(unlist(sub[i, PARAM]), decreasing = TRUE)
  cat(sprintf("%-14s %s\n", sub$site[i], paste(sprintf("%s (%.2f)", names(idx), idx), collapse = " > ")))
}
sink()
