# time_runs.R -- sensitivity step 1: time N serial colonization runs at
# literature parameters, Maquipucuna, against the corrected niche cache, so
# the LHS sample size can be scaled from a measured seconds/run.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("time_runs_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))))

N <- as.integer(Sys.getenv("CANOPY_TIMING_N", unset = "20"))
SITE <- "Maquipucuna"
Sys.setenv(CANOPY_CLIM_MODE = "voxel")

microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", SITE)))
niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
canopy_ceiling <- site_canopy_ceiling(SITE, niches, max(microenv_heights(microenv)))
canopy_grid <- matrix(canopy_ceiling, nrow = 50, ncol = 50)
site <- list(Site = SITE)
forestparams <- site_forestparams(SITE, canopy_ceiling)

params <- list(
  beta0S = -0.24 + 2.889, beta0J = 0.41 + 2.729, beta0A = 1.73 + 2.563, beta1 = 0.10,
  s_S_min = 0.0, s_S_max = 1.0, s_J_min = 1.0, s_J_max = 7.0, s_A_min = 7.0, s_A_max = 20.0,
  psi0S = -3.30 - 2.577, psi0J = -2.70 - 2.619, beta_precip = 3e-4, beta_rh = 0.010,
  sigma = 0.10, delta_s_base = 0.80, cost_repro = 0.50,
  S = 1.76e6, p_poll = 0.30, p_germ = 0.001, p_s1 = 0.45,
  canopy_z = canopy_ceiling, lambda = 1, Ut = 1, n_founders = 30
)

clim_cache <- build_clim_cache(microenv)

times <- numeric(N); totals <- numeric(N); regimes <- character(N)
for (i in seq_len(N)) {
  t0 <- Sys.time()
  r <- suppressMessages(runcolonization(
    site = site, niches = niches, canopy_grid = canopy_grid, microenv = microenv,
    timesteps = 50, resolution = 10, carCap = 1, maxDisp = 5,
    stochastic = FALSE, Visualize = FALSE, spinup = 5,
    parameters = params, forestparams = forestparams, seed = i, clim_cache = clim_cache
  ))
  times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  totals[i] <- tail(r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA, 1)
  regimes[i] <- r$regime
  log_msg(sprintf("run %2d/%d: %.1f s | final total %.0f | regime %s | k_bind_t %s",
                  i, N, times[i], totals[i], regimes[i],
                  if (is.na(r$k_bind_t)) "NA" else r$k_bind_t))
}

cat(sprintf("\n=== TIMING: %d serial runs, %s, literature params, 50 ts, voxel mode ===\n", N, SITE))
cat(sprintf("per-run seconds: mean %.1f | sd %.1f | min %.1f | max %.1f | median %.1f\n",
            mean(times), sd(times), min(times), max(times), median(times)))
cat(sprintf("(first run %.1f s -- includes cache warm-up; runs 2..N mean %.1f s)\n",
            times[1], mean(times[-1])))
cat(sprintf("final totals: mean %.0f (sd %.0f) | regimes: %s\n",
            mean(totals), sd(totals), paste(names(table(regimes)), table(regimes), sep = "x", collapse = " ")))

sec <- mean(times[-1])
cat(sprintf("\nseconds/run for budgeting = %.1f\n", sec))
for (nm in c("Maq_3000", "Mashpi_900", "total_3900", "halved_2400")) {
  n <- switch(nm, Maq_3000 = 3000, Mashpi_900 = 900, total_3900 = 3900, halved_2400 = 2400)
  cat(sprintf("  %-12s %4d runs = %6.1f CPU-h | 40 concurrent %5.1f h | 128 concurrent %5.1f h\n",
              nm, n, n * sec / 3600, n * sec / 3600 / 40, n * sec / 3600 / 128))
}
cat("\nDone.\n")
