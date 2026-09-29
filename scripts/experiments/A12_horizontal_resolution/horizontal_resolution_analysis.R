# horizontal_resolution_analysis.R -- paired test of pooled vs voxel climate
# mode (A12). Replicates share seeds across modes (run_replicated(): seed =
# rep), so final total abundance is paired by (site, rep).
#   - per site: paired Wilcoxon signed-rank on the n_reps pairs (with 5 pairs
#     the smallest attainable two-sided p is 0.0625, so the effect size and
#     its bootstrap CI are the primary result);
#   - across sites: paired Wilcoxon on the per-site mean of each mode (n = sites).
source("scripts/02_model/config/paths.R")

SITES <- c("Maquipucuna", "Mashpi", "Yanayacu", "MindoMirador", "MindoTarabita", "Saloya")
N_BOOT <- 10000L
set.seed(20260930)

final_by_rep <- function(site, mode) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_horizontal_resolution_%s_h0.40.rds", site, mode))
  if (!file.exists(f)) return(NULL)
  s <- readRDS(f)$summary
  d <- s[s$t == max(s$t), c("rep", "total", "recruited")]
  names(d)[2:3] <- paste0(c("total_", "recruited_"), mode)
  d
}

pairs_all <- list(); rows <- list()
for (site in SITES) {
  p <- final_by_rep(site, "pooled"); v <- final_by_rep(site, "voxel")
  if (is.null(p) || is.null(v)) { message("Skipping ", site, " -- missing a mode"); next }
  d <- merge(p, v, by = "rep")
  d$site <- site
  d$ratio <- d$total_voxel / d$total_pooled
  pairs_all[[site]] <- d
  w <- tryCatch(wilcox.test(d$total_voxel, d$total_pooled, paired = TRUE, exact = TRUE), error = function(e) NULL)
  boot <- vapply(seq_len(N_BOOT), function(i) {
    k <- sample.int(nrow(d), replace = TRUE)
    mean(d$total_voxel[k]) / mean(d$total_pooled[k])
  }, numeric(1))
  rows[[site]] <- data.frame(
    site = site, n_pairs = nrow(d),
    mean_pooled = mean(d$total_pooled), mean_voxel = mean(d$total_voxel),
    ratio_voxel_over_pooled = mean(d$total_voxel) / mean(d$total_pooled),
    ratio_ci_lo = quantile(boot, 0.025, names = FALSE), ratio_ci_hi = quantile(boot, 0.975, names = FALSE),
    wilcoxon_V = if (!is.null(w)) unname(w$statistic) else NA, wilcoxon_p = if (!is.null(w)) w$p.value else NA,
    recruited_pooled = sum(d$recruited_pooled), recruited_voxel = sum(d$recruited_voxel))
}
res <- do.call(rbind, rows)
res$wilcoxon_p_holm <- p.adjust(res$wilcoxon_p, method = "holm")
print(res, row.names = FALSE, digits = 4)

cat("\n== Across sites: paired Wilcoxon on per-site means (n = sites) ==\n")
wa <- wilcox.test(res$mean_voxel, res$mean_pooled, paired = TRUE, exact = TRUE)
print(wa)

write.csv(res, file.path(OUTPUT_DIR, "horizontal_resolution_paired_test.csv"), row.names = FALSE)
write.csv(do.call(rbind, pairs_all), file.path(OUTPUT_DIR, "horizontal_resolution_pairs.csv"), row.names = FALSE)
cat(sprintf("\nAcross-site Wilcoxon: V = %.1f, p = %.4f (n = %d sites)\nDone.\n", wa$statistic, wa$p.value, nrow(res)))
