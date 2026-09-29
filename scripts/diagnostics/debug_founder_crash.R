# debug_founder_crash.R -- reproduce the "NAs are not allowed in subscripted
# assignments" crash from the founder_number sweep (job 27106698,
# colonization_Maquipucuna_founder_number_h0.40_20260823_113447.log, 60/60
# jobs failed identically regardless of n_founders value) OUTSIDE
# run_one()'s tryCatch, so the error's real call stack survives instead of
# being swallowed. Single core, single value, single rep.
options(error = function() { traceback(2); quit(save = "no", status = 1) })

site_name <- "Maquipucuna"
params_file <- "data/params/n_founders.rds"
height_step <- 0.4

library(parallel)
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

manifest_suffix <- if (height_step != 0.1) sprintf("_h%.2f", height_step) else ""
microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, manifest_suffix))
microenv <- readRDS(microenv_path)
available_heights <- microenv_heights(microenv)
cat(sprintf("Heights available: %d (%.2f-%.2fm)\n", length(available_heights),
            min(available_heights), max(available_heights)))

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == site_name], na.rm = TRUE)
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)
site <- list(Site = site_name)
forestparams <- default_forestparams()

params <- readRDS(params_file)
params$canopy_z <- mean_canopy
params$n_founders <- 30  # pin one value from the swept vector for reproduction

log_file <- tempfile()
log_msg <- function(msg) cat(msg, "\n")

clim_cache <- build_clim_cache(microenv)

cat("\n=== Calling runcolonization() directly (no tryCatch) ===\n")
r <- runcolonization(
  site = site, niches = niches, canopy_grid = canopy_grid, microenv = microenv,
  timesteps = 50, resolution = 10, carCap = 5, maxDisp = 10, spinup = 5,
  Visualize = FALSE, parameters = params, forestparams = forestparams,
  clim_cache = clim_cache, seed = 1
)
cat("SUCCEEDED -- no crash reproduced with this seed/value.\n")
