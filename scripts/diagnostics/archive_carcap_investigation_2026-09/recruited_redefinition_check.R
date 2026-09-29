# recruited_redefinition_check.R
# Recomputes "recruited" under the new final-timestep-only definition
# directly from each site's already-saved best_combo_v6 result (per-timestep
# totalS/totalJ are already in $summary -- no rerun of the model needed),
# and reports the mean total-abundance trajectory (t=1, t=25, t=50) per
# replicate-averaged site. Read-only.
source("scripts/02_model/config/paths.R")

sites <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Yanayacu")

for (s in sites) {
  path <- file.path(PROCESSED_DIR, sprintf("colonization_%s_best_combo_v6_h0.40.rds", s))
  if (!file.exists(path)) { message(s, ": no result yet, skipping"); next }
  res <- readRDS(path)
  d <- res$summary
  n_timesteps <- max(d$t)

  by_rep_old <- sapply(split(d, d$rep), function(sub) {
    half <- (n_timesteps %/% 2):n_timesteps
    any(c(sub$totalS[sub$t %in% half], sub$totalJ[sub$t %in% half]) > 0)
  })
  by_rep_new <- sapply(split(d, d$rep), function(sub) {
    (sub$totalS[sub$t == n_timesteps] + sub$totalJ[sub$t == n_timesteps]) > 0
  })

  cat(sprintf("\n=== %s (best_combo_v6, n_timesteps=%d) ===\n", s, n_timesteps))
  cat(sprintf("  Recruited (OLD, any S/J in 2nd half): %d/%d replicates\n", sum(by_rep_old), length(by_rep_old)))
  cat(sprintf("  Recruited (NEW, S/J>0 at final t):    %d/%d replicates\n", sum(by_rep_new), length(by_rep_new)))

  traj <- aggregate(cbind(totalS, totalJ, totalA, total) ~ t, d, mean)
  traj <- traj[traj$t %in% c(1, round(n_timesteps / 2), n_timesteps), ]
  cat("  Mean trajectory (t=1, mid, final):\n")
  print(traj, row.names = FALSE)
}
