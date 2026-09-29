# build_site_climate_cache_v2.R -- 2026-09-28. One-off cache build for the
# 7-site microclimate re-analysis (VPD/elevation/vertical-structure task).
# Reuses .build_site_climate_series() (shared_helpers.R) -- the SITE-LEVEL
# (whole-footprint spatially-averaged) hourly series per height tier, the
# only hourly-resolution data that exists in this pipeline (per-pixel
# hourly data was discarded at write time; only monthly/annual per-pixel
# quantiles survive -- see docs note this script's caller writes).
# Output: data/processed/site_climate_cache_v2.rds
#   list(site_pixels = list(site -> data.frame(height, temp, relhum, swdown)),
#        site_elev = named vector, tme = POSIXct[4392] (from Maquipucuna's
#        own height files; sites' own tme are stored per-height in
#        tme_by_site), swdown_negclip = data.frame(site, n_negative, min_raw))
suppressPackageStartupMessages({})
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

SITES <- c("LaElenita", "Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")

Sys.setenv(CANOPY_OBS_CSV = Sys.getenv("CANOPY_OBS_CSV", "data/csv/combinedv6.csv"))
niches_raw <- read.csv(OBSERVATIONS_CSV, stringsAsFactors = FALSE)
niches_raw <- niches_raw[!is.na(niches_raw$lat) & !is.na(niches_raw$lon), ]
niches_elev <- augment_elevation(niches_raw)

# also need FinalID-filtered niches for .build_site_climate_series()'s own
# (separately-handled) site_elev -- but we pass site_elev explicitly below
# from niches_elev (matches shared_helpers.R's own 2026-09-05 LaElenita fix
# rationale: elevation must not depend on identification-confirmed rows).
niches_id <- niches_raw[!is.na(niches_raw$Height_m) & !is.na(niches_raw$FinalID), ]

cat("Building site_pixels (per-height, whole-footprint hourly series) -- this reads every h*.rds under each site's height_dir; heavy I/O, ~220GB total.\n")
t0 <- Sys.time()
sc <- .build_site_climate_series(SITES, niches = niches_id)
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

site_pixels <- sc$site_pixels
# Overwrite site_elev with the minimally-filtered version (lat/lon only,
# not identification-confirmed) -- consistent with shared_helpers.R's own
# fix note; ensures every site (incl. LaElenita) gets a real elevation.
site_elev <- sapply(SITES, function(s) mean(niches_elev$Elevation_final_m[niches_elev$Area_or_Site == s], na.rm = TRUE))

# ---- negative-swdown audit + clip (Step 1.1-adjacent, per user request 'd') --
negclip <- do.call(rbind, lapply(names(site_pixels), function(s) {
  px <- site_pixels[[s]]
  n_neg <- sum(px$swdown < 0, na.rm = TRUE)
  min_raw <- suppressWarnings(min(px$swdown, na.rm = TRUE))
  data.frame(site = s, n_negative = n_neg, min_raw = min_raw, n_total = nrow(px))
}))
for (s in names(site_pixels)) {
  site_pixels[[s]]$swdown <- pmax(site_pixels[[s]]$swdown, 0)
}

# ---- per-height tme (needed since sites don't share a calendar window) ----
tme_by_site <- list()
for (s in SITES) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", s))
  if (!file.exists(microenv_path)) next
  microenv <- readRDS(microenv_path)
  h1 <- sort(microenv_heights(microenv))[1]
  hd <- microenv$.height_dir
  fn <- sprintf("h%.2f.rds", h1)
  x <- readRDS(file.path(hd, fn))
  tme_by_site[[s]] <- x$tme
}

out <- list(site_pixels = site_pixels, site_elev = site_elev,
            tme_by_site = tme_by_site, swdown_negclip = negclip,
            built = Sys.time())
saveRDS(out, file.path(PROCESSED_DIR, "site_climate_cache_v2.rds"))
cat("Saved: ", file.path(PROCESSED_DIR, "site_climate_cache_v2.rds"), "\n")
cat("\nnegative swdown audit (pre-clip):\n"); print(negclip)
cat("\nsite_elev:\n"); print(round(site_elev))
cat("\nn heights per site:\n"); print(sapply(site_pixels, function(px) length(unique(px$height))))
