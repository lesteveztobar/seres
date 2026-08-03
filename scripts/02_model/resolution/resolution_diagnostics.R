# resolution_diagnostics.R
# Efficiency diagnostic: times the climate-cache build (driven by height-tier
# spacing) and a short colonization run (driven by horizontal voxel
# resolution) across several resolution combinations, to find the coarsest
# resolution that's still fine enough before committing to full experiment
# runs at the production 10m x 0.1m resolution.
#
# For each height step, requires microenv_<site>.rds (step 0.1, production)
# or microenv_<site>_h<step>.rds (coarser steps, from height_res_array.sh) to
# already exist — this script does NOT generate microclimate data itself.
#
# Usage: Rscript scripts/02_model/resolution/resolution_diagnostics.R <site> [height_steps] [horiz_resolutions]
#   e.g.: Rscript scripts/02_model/resolution/resolution_diagnostics.R Maquipucuna "0.1,0.25,0.5,1.0" "5,10,20,40"
# Output: output/resolution_diagnostics_<site>.csv
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(parallel)
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

args <- commandArgs(trailingOnly = TRUE)
site_name  <- if (length(args) >= 1) args[1] else "Maquipucuna"
height_steps <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(strsplit(args[2], ",")[[1]]) else c(0.1, 0.25, 0.5, 1.0)
horiz_res    <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(strsplit(args[3], ",")[[1]]) else c(5, 10, 20, 40)
DIAG_TIMESTEPS <- 5   # short run — timing only, not meant to be scientifically meaningful
DIAG_SPINUP    <- 2

N_CORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)

log_msg <- function(msg) message("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

# ── Site setup (mirrors run_colonization.R) ─────────────────────────────
niches <- load_observations()
niches <- niches[
  !is.na(niches$lat) & !is.na(niches$lon) &
  !is.na(niches$Height_m) & !is.na(niches$FinalID), ]

mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == site_name], na.rm = TRUE)
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)
site <- list(Site = site_name)

forestparams <- list(
  stems_per_ha = 298, mean_hgt = 8.4, sd_hgt = 3.5,
  mean_crown_r = 2.0, sd_crown_r = 0.8, trunk_r = 0.114,
  branch_density = 3.0, epiphyte_footprint_m2 = 0.02
)

params <- list(
  beta0S = -0.24, beta0J = 0.41, beta0A = 1.73, beta1 = 0.10,
  s_S_min = 0.0, s_S_max = 1.0, s_J_min = 1.0, s_J_max = 7.0,
  s_A_min = 7.0, s_A_max = 20.0,
  psi0S = -3.30, psi0J = -2.70, beta_precip = 3e-4, beta_rh = 0.010,
  sigma = 0.10, delta_s_base = 0.80, cost_repro = 0.50,
  p_poll = 0.30, p_germ = 0.001, p_s1 = 0.45,
  canopy_z = mean_canopy, lambda = 1, Ut = 1
)

# ── Run one height-step's worth of combinations ─────────────────────────────────
run_height_step <- function(step) {
  suffix <- if (step != 0.1) sprintf("_h%.2f", step) else ""
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, suffix))
  if (!file.exists(microenv_path)) {
    log_msg(sprintf("Skipping height step %.2fm -- no %s (run height_res_array.sh first)",
                    step, microenv_path))
    return(NULL)
  }

  microenv <- readRDS(microenv_path)
  n_heights <- length(microenv_heights(microenv))
  log_msg(sprintf("Height step %.2fm: %d tiers. Building climate cache...", step, n_heights))

  t0 <- Sys.time()
  clim_cache <- build_clim_cache(microenv)
  clim_build_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_msg(sprintf("Height step %.2fm: climate cache built in %.1fs", step, clim_build_sec))

  # Horizontal resolution doesn't touch the climate cache, so these can run
  # in parallel, all sharing the one clim_cache we just built.
  rows <- mclapply(horiz_res, function(res) {
    t1 <- Sys.time()
    result <- tryCatch(
      runcolonization(
        site = site, niches = niches, canopy_grid = canopy_grid, microenv = microenv,
        timesteps = DIAG_TIMESTEPS, resolution = res, carCap = 5, maxDisp = 10,
        spinup = DIAG_SPINUP, Visualize = FALSE, parameters = params,
        forestparams = forestparams, clim_cache = clim_cache
      ),
      error = function(e) { message("ERROR [height=", step, " res=", res, "]: ", e$message); NULL }
    )
    sim_sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    data.frame(
      height_step = step, n_heights = n_heights, clim_build_sec = clim_build_sec,
      horiz_resolution = res,
      xDim = if (!is.null(result)) result$xDim else NA,
      yDim = if (!is.null(result)) result$yDim else NA,
      sim_sec = sim_sec, ok = !is.null(result)
    )
  }, mc.cores = N_CORES)

  do.call(rbind, rows)
}

all_rows <- lapply(height_steps, run_height_step)
summary_df <- do.call(rbind, Filter(Negate(is.null), all_rows))

if (is.null(summary_df) || nrow(summary_df) == 0) {
  stop("No height-step variants available — run height_res_array.sh for ", site_name, " first.")
}

summary_df$total_sec <- summary_df$clim_build_sec + summary_df$sim_sec
summary_df <- summary_df[order(summary_df$height_step, summary_df$horiz_resolution), ]

out_path <- file.path(OUTPUT_DIR, sprintf("resolution_diagnostics_%s.csv", site_name))
write.csv(summary_df, out_path, row.names = FALSE)

cat("\n== Resolution diagnostic summary ==\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_path))
