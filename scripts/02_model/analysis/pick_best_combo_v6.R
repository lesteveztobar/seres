source("scripts/02_model/config/paths.R")
sites <- c("Maquipucuna","Mashpi","MindoMirador","MindoTarabita","Yanayacu")
base_params <- list(
    beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
    beta1   =  0.10,
    s_S_min =  0.0,  s_S_max =  1.0,
    s_J_min =  1.0,  s_J_max =  7.0,
    s_A_min =  7.0,  s_A_max = 20.0,
    psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_s_base =  0.80,
    cost_repro   =  0.50,
    lambda = 3.23,  Ut = 0.23,
    n_reps = 3
)
for (site in sites) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_reproduction_factorial_v3_v6_h0.40.rds", site))
  if (!file.exists(f)) { message("Skipping ", site, " -- no ", f); next }
  d <- readRDS(f)
  final_t <- max(d$t)
  final <- d[d$t == final_t, ]
  agg <- aggregate(total ~ p_poll + p_germ + p_s1 + n_founders, data = final, FUN = mean)
  best <- agg[which.max(agg$total), ]
  cat(sprintf("%-15s best combo: p_poll=%.2f p_germ=%.4f p_s1=%.2f n_founders=%d (mean final total=%.1f, from %d combos)\n",
    site, best$p_poll, best$p_germ, best$p_s1, best$n_founders, best$total, nrow(agg)))
  params <- base_params
  params$p_poll <- best$p_poll
  params$p_germ <- best$p_germ
  params$p_s1   <- best$p_s1
  params$n_founders <- best$n_founders
  out_path <- file.path(PARAMS_DIR, sprintf("best_combo_v6_%s.rds", site))
  saveRDS(params, out_path)
  cat(sprintf("  Saved to %s\n", out_path))
}
