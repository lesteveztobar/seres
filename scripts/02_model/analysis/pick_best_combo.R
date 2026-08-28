# pick_best_combo.R — for each site's reproduction_factorial_v3 result, find the
# (p_poll, p_germ, p_s1, n_founders) combination with the highest mean final-year
# total abundance, and save a params RDS with that combo substituted in, ready to
# feed straight to run_colonization.R as a single representative run (no swept
# parameter -> hits its `length(swept_params) == 0` branch -> run_replicated()
# saves the full per-species abundanceS/J/A arrays needed for
# plot_species_stage_curves(), unlike the sweep .rds files themselves which only
# keep pooled totals).
#
# Usage: Rscript scripts/02_model/analysis/pick_best_combo.R [site1] [site2] ...
#   Defaults to every site with a reproduction_factorial_v3_h0.40.rds on disk.
source("scripts/02_model/config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  sites <- args
} else {
  files <- list.files(PROCESSED_DIR,
    pattern = "^colonization_.*_reproduction_factorial_v3_h0\\.40\\.rds$")
  sites <- sub("^colonization_(.*)_reproduction_factorial_v3_h0\\.40\\.rds$", "\\1", files)
}

# Base (non-swept) params -- same literature-default values make_params.R uses
# for every OAT/factorial params list (see e.g. its params_reprofactorial_v3
# block); only p_poll/p_germ/p_s1/n_founders get overridden below per-site.
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
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_reproduction_factorial_v3_h0.40.rds", site))
  if (!file.exists(f)) {
    message("Skipping ", site, " -- no ", f)
    next
  }
  d <- readRDS(f)
  final_t <- max(d$t)
  final <- d[d$t == final_t, ]
  agg <- aggregate(total ~ p_poll + p_germ + p_s1 + n_founders, data = final, FUN = mean)
  best <- agg[which.max(agg$total), ]

  cat(sprintf(
    "%-15s best combo: p_poll=%.2f p_germ=%.4f p_s1=%.2f n_founders=%d (mean final total=%.1f, from %d combos)\n",
    site, best$p_poll, best$p_germ, best$p_s1, best$n_founders, best$total, nrow(agg)))

  params <- base_params
  params$p_poll <- best$p_poll
  params$p_germ <- best$p_germ
  params$p_s1   <- best$p_s1
  params$n_founders <- best$n_founders

  out_path <- file.path(PARAMS_DIR, sprintf("best_combo_%s.rds", site))
  saveRDS(params, out_path)
  cat(sprintf("  Saved to %s\n", out_path))
}
