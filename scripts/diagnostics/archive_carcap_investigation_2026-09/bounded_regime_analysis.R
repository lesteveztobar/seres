# bounded_regime_analysis.R -- follow-ups on the literature-parameter
# diagnostic once bounded persistence was confirmed:
#   (1) extend ONE Maquipucuna replicate to 150 timesteps and report
#       whether/where the trajectory plateaus (t=50 was ~86% of a fitted
#       asymptote -- not a demonstration);
#   (2) at t=50, the fraction of OCCUPIABLE voxels (K>0) that are occupied,
#       and the per-voxel occupancy-vs-K distribution -- to test whether
#       the population saturates well below structural capacity because
#       suitable microenvironment, not bark area, is limiting.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("bounded_regime_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))))

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
params <- readRDS("data/params/literature_diagnostic.rds")
params$n_reps <- NULL
params$canopy_z <- canopy_ceiling

# ---- (1) 150-timestep single replicate --------------------------------------
log_msg("Running Maquipucuna, literature params, seed=1, 150 timesteps, voxel mode...")
res <- runcolonization(
  site = site, niches = niches, canopy_grid = canopy_grid, microenv = microenv,
  timesteps = 150, resolution = 10, carCap = 1, maxDisp = 5,
  stochastic = FALSE, Visualize = FALSE, spinup = 5,
  parameters = params, forestparams = forestparams, seed = 1
)
tot <- res$totalabundanceS + res$totalabundanceJ + res$totalabundanceA
traj <- data.frame(t = seq_along(tot), totalS = res$totalabundanceS,
                   totalJ = res$totalabundanceJ, totalA = res$totalabundanceA, total = tot)
write.csv(traj, file.path(OUTPUT_DIR, "bounded_regime_maquipucuna_150ts.csv"), row.names = FALSE)

cat("\n=== (1) Maquipucuna 150-timestep trajectory (seed 1) ===\n")
show_t <- c(1, 5, 10, 25, 50, 75, 100, 125, 150)
print(traj[traj$t %in% show_t, ])
incr <- diff(tot)
cat(sprintf("\nyearly increment: t=50 %.1f | t=75 %.1f | t=100 %.1f | t=125 %.1f | t=149 %.1f\n",
            incr[49], incr[74], incr[99], incr[124], incr[149]))
last10 <- mean(incr[141:149])
cat(sprintf("mean increment over final 10 years: %.2f/yr (%.3f%% of t=150 total)\n",
            last10, 100 * last10 / tot[150]))
cat(sprintf("t=50 total = %.0f  |  t=150 total = %.0f  |  t=50 is %.1f%% of t=150\n",
            tot[50], tot[150], 100 * tot[50] / tot[150]))

# ---- (2) per-voxel occupancy vs K at t=50 (from this run) ------------------
cc <- res$state$carCap_voxel
land <- res$state$landscape
occ50 <- apply(res$abundanceS[, , , 50, , drop = FALSE], 1:3, sum) +
         apply(res$abundanceJ[, , , 50, , drop = FALSE], 1:3, sum) +
         apply(res$abundanceA[, , , 50, , drop = FALSE], 1:3, sum)
occ150 <- apply(res$abundanceS[, , , 150, , drop = FALSE], 1:3, sum) +
          apply(res$abundanceJ[, , , 150, , drop = FALSE], 1:3, sum) +
          apply(res$abundanceA[, , , 150, , drop = FALSE], 1:3, sum)

report_occ <- function(occ, label) {
  kpos <- which(cc > 0)
  n_kpos <- length(kpos)
  n_occ_kpos <- sum(occ[kpos] > 0)
  n_land <- sum(land)
  n_occ_land <- sum(occ[land] > 0)
  cat(sprintf("\n=== (2) per-voxel occupancy vs K -- %s ===\n", label))
  cat(sprintf("landscape voxels: %d | occupiable (K>0): %d | total K: %.0f\n", n_land, n_kpos, sum(cc)))
  cat(sprintf("occupiable voxels that are occupied: %d / %d = %.1f%%\n", n_occ_kpos, n_kpos, 100 * n_occ_kpos / n_kpos))
  cat(sprintf("all landscape voxels that are occupied: %d / %d = %.1f%%\n", n_occ_land, n_land, 100 * n_occ_land / n_land))
  cat(sprintf("total N = %.0f | total K = %.0f | N/K = %.3f\n", sum(occ), sum(cc), sum(occ) / sum(cc)))
  # occupancy relative to local K, among occupied occupiable voxels
  ratio <- occ[kpos] / cc[kpos]
  cat(sprintf("among occupiable voxels: mean occ/K = %.2f | median = %.2f | at/over local K: %d (%.1f%%)\n",
              mean(ratio), median(ratio), sum(ratio >= 1), 100 * sum(ratio >= 1) / n_kpos))
  cat("occ/K distribution (occupiable voxels): ")
  print(round(quantile(ratio, c(.5, .75, .9, .95, .99, 1)), 2))
  # occupied voxels that carry NO K (established outside bark-capacity voxels
  # -- shouldn't happen given the establishment gate, a sanity check)
  cat(sprintf("occupied voxels with K==0: %d\n", sum(occ > 0 & cc == 0)))
}
report_occ(occ50, "t=50")
report_occ(occ150, "t=150")

saveRDS(list(cc = cc, land = land, occ50 = occ50, occ150 = occ150, traj = traj),
        file.path(PROCESSED_DIR, "bounded_regime_maquipucuna_voxels.rds"))
cat("\nDone.\n")
