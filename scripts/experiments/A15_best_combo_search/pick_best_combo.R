# pick_best_combo.R -- reads the 32-combo targeted-search factorial result
# per site, picks the combo with the highest RECRUITED fraction (the
# locked persistence metric, methods_update_report.md) across its 3
# replicates, tie-broken by mean final total abundance. Writes one
# single-value params RDS per site (data/params/best_combo_<site>.rds) --
# NOT the winning combo's own factorial-summary rows, which is why this
# has to be a re-run: run_factorial_experiment() only returns a summary
# data.frame (totalS/J/A, extinct, recruited per rep x timestep), not the
# per-replicate obs_val/abundance arrays held_out_validation.R needs --
# those only come from a real single-params run_replicated() call.
source("scripts/02_model/config/paths.R")

SITES <- c("Maquipucuna", "Mashpi", "Yanayacu", "MindoMirador", "MindoTarabita", "Saloya")
base <- readRDS("data/params/best_combo_targeted_search.rds")
swept <- c("p_poll", "p_germ", "p_s1", "S", "n_founders")

for (site in SITES) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_best_combo_search_h0.40.rds", site))
  if (!file.exists(f)) { message("Skipping ", site, " -- no ", f, " (job likely timed out/failed -- see checkpoint)"); next }
  df <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    message("Skipping ", site, " -- ", f, " exists but has 0 usable rows (every replicate failed).")
    next
  }
  # One row per (combo, rep, t) -- reduce to per-combo: recruited fraction
  # across reps at the FINAL timestep, mean final total as tiebreak.
  final_t <- max(df$t)
  final_rows <- df[df$t == final_t, ]
  combo_key <- do.call(paste, c(final_rows[swept], sep = "|"))
  agg <- aggregate(cbind(recruited, total) ~ combo_key, data = data.frame(final_rows, combo_key = combo_key),
                   FUN = function(x) x)
  frac_recruited <- sapply(agg$recruited, function(x) mean(as.logical(x)))
  mean_total <- sapply(agg$total, mean)
  best_i <- order(-frac_recruited, -mean_total)[1]
  best_key <- agg$combo_key[best_i]
  best_row <- final_rows[combo_key == best_key, ][1, ]

  params_best <- base
  for (nm in swept) params_best[[nm]] <- best_row[[nm]]
  params_best$n_reps <- 5
  out_path <- file.path("data/params", sprintf("best_combo_%s.rds", site))
  saveRDS(params_best, out_path)
  cat(sprintf("%-15s best combo: p_poll=%.3f p_germ=%.5f p_s1=%.3f S=%.3g n_founders=%.0f | recruited_frac=%.2f mean_total=%.1f -> %s\n",
              site, best_row$p_poll, best_row$p_germ, best_row$p_s1, best_row$S, best_row$n_founders,
              frac_recruited[best_i], mean_total[best_i], out_path))
}
