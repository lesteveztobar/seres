# founder_number_summary.R
# Clean summary of the post-fix founder_number_v6 sweep, 5 sites: recruited
# fraction (final-timestep S+J>0 definition) and final total abundance per
# n_founders level. Avoids .classify_result_shape()'s "recruited" column
# being misread as a second swept parameter (a real, separate bug in that
# function -- not fixed here, worked around by reading result$summary
# directly).
source("scripts/02_model/config/paths.R")

sites <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Yanayacu")
for (s in sites) {
  path <- file.path(PROCESSED_DIR, sprintf("colonization_%s_founder_number_v6_h0.40.rds", s))
  if (!file.exists(path)) { message(s, ": no result"); next }
  d <- readRDS(path)
  final_t <- max(d$t)
  final <- d[d$t == final_t, ]
  agg <- aggregate(cbind(total, recruited) ~ param_value, final, mean)
  agg <- agg[order(as.numeric(as.character(agg$param_value))), ]
  cat(sprintf("\n=== %s (final t=%d, n=%d levels) ===\n", s, final_t, nrow(agg)))
  names(agg) <- c("n_founders", "mean_total_final", "frac_recruited_final")
  print(agg, row.names = FALSE)
  cat(sprintf("Overall recruited fraction across all levels: %.1f%% (%d/%d combos)\n",
              100 * mean(final$recruited), sum(final$recruited), nrow(final)))
}
