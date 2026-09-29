# fetch_vhgt_site.R
# Standalone Lang et al. (2023) / ETH GlobalCanopyHeight export for ONE site,
# bypassing the full get_vegetation() pipeline (lcover/lai/refldata not
# needed for this) -- used for Task 6 (site-specific mean tree height).
# Exports via GEE to Drive, downloads, computes mean/median/sd tree height
# within the site's own observation bounding box (buffered), saves both the
# raw raster and a one-row summary CSV. Read-only w.r.t. the model/observations.
#
# Usage: Rscript scripts/diagnostics/fetch_vhgt_site.R <SiteName>
source("scripts/02_model/config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
site_name <- args[1]
if (is.na(site_name)) stop("Usage: Rscript fetch_vhgt_site.R <SiteName>")

niches <- read.csv(OBSERVATIONS_CSV)
obs <- niches[niches$Area_or_Site == site_name & !is.na(niches$lat) & !is.na(niches$lon), ]
if (nrow(obs) == 0) stop("No coordinates for ", site_name)

# 2026-09-02 (v7 rebuild, Phase 1.3): TIGHT buffer around the plot's own
# centroid, not the raw observation bounding box -- the original version of
# this script buffered whatever range min/max(lat/lon) happened to span,
# which for several sites (raw coordinate spread genuinely wide, or --
# Saloya specifically -- corrupted by the 2 mislabeled rows, since relabeled
# to MindoMirador) produced a box on the order of 1,100 km^2, no different
# from sampling a broad regional average (confirmed: Maquipucuna and Mashpi
# both returned ~27m under the old buffering, indistinguishable from each
# other despite being different sites). MEDIAN (not mean, robust to any
# remaining outliers) lat/lon as centroid, fixed 300m half-width -- "a few
# hundred metres around actual plot coordinates", matching the scale of a
# real field survey plot (Task 0's own landscape-extent audit: 7-33 ha at
# every genuine site).
buf_deg <- 300 / 111000
lat0 <- median(obs$lat); lon0 <- median(obs$lon)
xmin <- lon0 - buf_deg; xmax <- lon0 + buf_deg
ymin <- lat0 - buf_deg; ymax <- lat0 + buf_deg
cat(sprintf("%s centroid (median): lat=%.6f lon=%.6f | tight bbox: lon [%.6f, %.6f] lat [%.6f, %.6f]\n",
            site_name, lat0, lon0, xmin, xmax, ymin, ymax))

reticulate::use_python(Sys.getenv("CANOPY_PYTHON",
  unset = "/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python3.12"), required = TRUE)
rgee::ee_Initialize(project = "ee-lizethestevezt")

aoi <- rgee::ee$Geometry$Rectangle(c(xmin, ymin, xmax, ymax))
canopy_height <- rgee::ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")

# Saved under a name distinct from "vhgt.tif" -- that filename is what
# get_vegetation()/the production microclimate pipeline reads (lib.R), and
# per Phase 1.3 this tight-buffer extraction is for forest STRUCTURE
# parameters only; it must never be picked up there. Drive prefix is
# likewise distinct from the original wide-bbox export's
# "canopy_height_<site>" name (still sitting in the rgee_backup Drive
# folder from the pre-v7 run) -- reusing that prefix would have `grepl()`
# match the stale wide-area file below and silently skip the new export.
out_dir <- file.path(RAW_DIR, site_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
vhgt_file <- file.path(out_dir, "vhgt_tight_forestparams.tif")

drive_prefix <- paste0("canopy_height_tight_", site_name)
googledrive::drive_auth(email = "lizethestevezt@gmail.com", cache = "~/.secrets")
folder <- googledrive::drive_find(pattern = "rgee_backup", type = "folder", n_max = 1)
if (nrow(folder) == 0) {
  folder <- googledrive::drive_mkdir("rgee_backup")
}
drive_files <- googledrive::drive_ls(folder)
drive_file <- drive_files[grepl(drive_prefix, drive_files$name), ]

if (nrow(drive_file) == 0) {
  message("Exporting from GEE for ", site_name, "...")
  task <- rgee::ee$batch$Export$image$toDrive(
    image = canopy_height,
    description = paste0("canopy_height_tight_", site_name),
    folder = "rgee_backup",
    fileNamePrefix = drive_prefix,
    region = aoi,
    scale = 10,
    crs = "EPSG:4326"
  )
  task$start()
  # Poll task status directly rather than relying on microclimdata's internal
  # monitor (avoids pulling in the full get_vegetation() dependency chain).
  repeat {
    Sys.sleep(15)
    status <- task$status()
    st <- status$state
    message("  GEE task state: ", st)
    if (st %in% c("COMPLETED", "FAILED", "CANCELLED")) break
  }
  if (st != "COMPLETED") stop("GEE export did not complete: ", st)
  drive_files <- googledrive::drive_ls(folder)
  drive_file <- drive_files[grepl(drive_prefix, drive_files$name), ]
  if (nrow(drive_file) == 0) stop("Export completed but file not found on Drive")
} else {
  message("Found existing Drive export for ", site_name, ", downloading...")
}

googledrive::drive_download(file = drive_file[1, ], path = vhgt_file, overwrite = TRUE)
message("Downloaded to ", vhgt_file)

vhgt <- terra::rast(vhgt_file)
vals <- terra::values(vhgt, na.rm = TRUE)
cat(sprintf("\n%s: n_pixels=%d  mean=%.2f  median=%.2f  sd=%.2f  p99=%.2f  max=%.2f\n",
            site_name, length(vals), mean(vals), median(vals), sd(vals),
            quantile(vals, 0.99), max(vals)))

summary_row <- data.frame(site = site_name, n_pixels = length(vals),
                           mean_hgt = mean(vals), median_hgt = median(vals),
                           sd_hgt = sd(vals), p99_hgt = quantile(vals, 0.99),
                           max_hgt = max(vals))
dir.create("output", showWarnings = FALSE)
out_csv <- "output/vhgt_site_summary_tight.csv"
if (file.exists(out_csv)) {
  prior <- read.csv(out_csv)
  prior <- prior[prior$site != site_name, ]
  summary_row <- rbind(prior, summary_row)
}
write.csv(summary_row, out_csv, row.names = FALSE)
message("Saved ", out_csv)
