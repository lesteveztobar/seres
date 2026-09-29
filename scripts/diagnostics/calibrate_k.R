# calibrate_k.R -- calibrates occupiable_bark_fraction so Maquipucuna's
# carrying capacity equals a literature stand density, writes
# data/params/k_calibration.rds (read by default_forestparams()), then
# reports the resulting K table for all 6 modelled sites. Exits non-zero if
# the calibration misses its target by more than TOL, so every colonisation
# job chained afterok on this one is cancelled rather than run on a bad K.
#
# Target (per ha, Maxillariinae adult-equivalents):
#   2,800 orchid stands/ha (Alzate-Q et al. 2019, Flora 260:151463, Veracruz)
#   x 41/1,348   (Maxillariinae share of Mexican orchid species)
#   x 271/41     (Ecuadorian vs Mexican Maxillariinae richness)
#   = 2,800 x 271 / 1,348 = 562.9
# The community share becomes 271/4,355 (POWO); the former x0.68 orchid-share
# factor is dropped because 2,800/ha already counts all orchids.
source("scripts/02_model/config/paths.R")
Sys.setenv(CANOPY_K_CALIBRATION = "")  # never read a stale calibration while calibrating
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("calibrate_k_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))))

TARGET_K_HA <- 2800 * 271 / 1348
MSHARE      <- 271 / 4355
TOL         <- 0.02
REF_SITE    <- "Maquipucuna"
SITES       <- c("Maquipucuna", "Mashpi", "Yanayacu", "MindoMirador", "MindoTarabita", "Saloya")
RESOLUTION  <- 10
CAL_PATH    <- file.path(PARAMS_DIR, "k_calibration.rds")

if (nzchar(Sys.getenv("CANOPY_K_MULT")) && Sys.getenv("CANOPY_K_MULT") != "1") stop("CANOPY_K_MULT must be unset for calibration")

obs <- load_observations()
obs <- obs[!is.na(obs$lat) & !is.na(obs$lon) & !is.na(obs$Height_m) & !is.na(obs$FinalID), ]

site_setup <- function(site) {
  microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site)))
  heights <- microenv_heights(microenv)
  so <- obs[obs$Area_or_Site == site, ]
  bb <- site_landscape_bbox(site)
  if (is.null(bb)) bb <- list(lat = range(so$lat), lon = range(so$lon), n = nrow(so))
  lat_m <- diff(bb$lat) * 111000
  lon_m <- diff(bb$lon) * 111000 * cos(mean(bb$lat) * pi / 180)
  list(heights = heights, so = so, bb = bb,
       xDim = max(round(lon_m / RESOLUTION), 10) + 4,
       yDim = max(round(lat_m / RESOLUTION), 10) + 4,
       area_ha = lat_m * lon_m / 10000,
       fp = site_forestparams(site, site_canopy_ceiling(site, obs, max(heights))))
}

site_K <- function(s, obf) {
  fp <- s$fp
  fp$occupiable_bark_fraction <- obf
  fp$maxillariinae_community_share <- MSHARE
  set.seed(1)
  landscape <- array(FALSE, dim = c(s$xDim, s$yDim, length(s$heights)))
  forest <- build_forest(landscape, s$heights, fp, s$so, RESOLUTION, land_bbox = s$bb)
  list(total_K = sum(forest$carCap_voxel), n_trees = forest$n_trees)
}

ref <- site_setup(REF_SITE)
base <- .default_forestparams_base()
obf <- base$occupiable_bark_fraction * base$maxillariinae_community_share / MSHARE * (TARGET_K_HA / 549)
for (it in 1:8) {
  k_ha <- site_K(ref, obf)$total_K / ref$area_ha
  err <- k_ha / TARGET_K_HA - 1
  cat(sprintf("iter %d: obf=%.6f  K/ha=%.1f  (target %.1f, err %+.2f%%)\n", it, obf, k_ha, TARGET_K_HA, 100 * err))
  if (abs(err) < 0.005) break
  obf <- obf / (1 + err)
}

if (abs(err) > TOL) {
  cat(sprintf("GATE FAILED: Maquipucuna K/ha %.1f is %.2f%% from target -- k_calibration.rds NOT written\n", k_ha, 100 * err))
  quit(status = 1)
}

saveRDS(list(occupiable_bark_fraction = obf, maxillariinae_community_share = MSHARE,
             target_K_per_ha = TARGET_K_HA, achieved_K_per_ha_ref = k_ha, ref_site = REF_SITE,
             source = "Alzate-Q et al. 2019 (2,800 orchid stands/ha) x 271/1,348",
             calibrated_on = format(Sys.time())), CAL_PATH)
cat("Wrote", CAL_PATH, "\n\n")

rows <- lapply(SITES, function(site) {
  s <- if (site == REF_SITE) ref else site_setup(site)
  k <- site_K(s, obf)
  data.frame(site = site, area_ha = round(s$area_ha, 2), n_trees = k$n_trees, total_K = k$total_K,
             K_per_tree = round(k$total_K / k$n_trees, 2), K_per_ha = round(k$total_K / s$area_ha, 1))
})
df <- do.call(rbind, rows)
print(df, row.names = FALSE)
write.csv(df, file.path(OUTPUT_DIR, "carcap_calibrated.csv"), row.names = FALSE)
cat("GATE PASSED\n")
