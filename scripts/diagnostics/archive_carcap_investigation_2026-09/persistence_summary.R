# persistence_summary.R
# Item 2: production-scale persistence validation (best_case.rds /
# realistic.rds), post-fix, 5 sites -- summary of final-timestep abundance,
# persisted, and recruited (final-timestep S+J>0 definition).
source("scripts/02_model/config/paths.R")

sites <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Yanayacu")
configs <- c("best_case_v6", "realistic_v6")

for (cfg in configs) {
  cat(sprintf("\n########## %s ##########\n", cfg))
  for (s in sites) {
    path <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s_h0.40.rds", s, cfg))
    if (!file.exists(path)) { message(s, ": no result"); next }
    res <- readRDS(path)
    d <- res$summary
    final_t <- max(d$t)
    final <- d[d$t == final_t, ]
    n_reps <- nrow(final)
    cat(sprintf("\n=== %s / %s (n_reps=%d, final t=%d) ===\n", s, cfg, n_reps, final_t))
    print(final[, c("rep", "totalS", "totalJ", "totalA", "total", "extinct", "recruited")],
          row.names = FALSE)
    traj <- aggregate(cbind(totalS, totalJ, totalA, total) ~ t, d, mean)
    traj <- traj[traj$t %in% c(1, round(final_t / 2), final_t), ]
    cat("Mean trajectory (t=1, mid, final):\n")
    print(traj, row.names = FALSE)
    cat(sprintf("Persisted: %d/%d | Recruited (final t): %d/%d\n",
                sum(!final$extinct), n_reps, sum(final$recruited), n_reps))
  }
}
