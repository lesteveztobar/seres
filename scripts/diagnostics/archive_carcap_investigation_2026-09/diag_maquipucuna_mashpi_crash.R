# diag_maquipucuna_mashpi_crash.R -- reproduces the "missing value where
# TRUE/FALSE needed" crash directly (single replicate, no mclapply) to get
# a real traceback instead of the swallowed error message.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("diag_crash_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))))

Sys.setenv(CANOPY_CLIM_MODE = "voxel")
site_name <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(site_name)) site_name <- "Mashpi"

microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site_name)))
niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
canopy_ceiling <- site_canopy_ceiling(site_name, niches, max(microenv_heights(microenv)))
canopy_grid <- matrix(canopy_ceiling, nrow = 50, ncol = 50)
site <- list(Site = site_name)
forestparams <- site_forestparams(site_name, canopy_ceiling)
params <- readRDS(Sys.getenv("CANOPY_DIAG_PARAMS", unset = "data/params/best_case.rds"))
params$canopy_z <- canopy_ceiling

options(error = function() { traceback(2); quit(status = 1) })
result <- runcolonization(
  site = site, niches = niches, canopy_grid = canopy_grid, microenv = microenv,
  timesteps = 50, resolution = 10, carCap = 1, maxDisp = 5,
  stochastic = FALSE, Visualize = FALSE, spinup = 5,
  parameters = params, forestparams = forestparams, seed = 1
)
cat("Completed without error.\n")
