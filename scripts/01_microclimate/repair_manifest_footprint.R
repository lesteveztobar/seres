# repair_manifest_footprint.R (2026-09-30)
# Adds the .footprint_idx / .site_name fields to a solo-site microenv manifest
# whose run died after saveRDS(manifest) but before the final block of
# run_microclimate_site.R (job 27989950: the script was edited in place while
# Rscript was still reading it -> "Error: unexpected end of input").
# Recomputes the footprint exactly as run_microclimate_site.R does (unpadded
# obs bbox -> .site_footprint_pixels() on the manifest's own .spatial).
# Usage: Rscript repair_manifest_footprint.R <Site> <res_suffix> <expected_px>
args <- commandArgs(trailingOnly = TRUE)
SITE <- args[1]; SUFFIX <- args[2]; EXPECTED <- as.integer(args[3])
suppressPackageStartupMessages({ library(readr); library(dplyr) })
source("scripts/02_model/config/paths.R")
source("scripts/data_prep/helper_functions.R")
source("scripts/02_model/engine/get_colonization.R")

path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", SITE, SUFFIX))
m <- readRDS(path)
if (!is.null(m$.footprint_idx)) { cat("already has .footprint_idx -- nothing to do\n"); quit(status = 0) }

df1 <- readr::read_csv(OBSERVATIONS_CSV, na = c("", "NA", "N/A"),
                       col_types = COMBINED_COL_TYPES, show_col_types = FALSE) |>
  dplyr::filter(!is.na(Source), Area_or_Site == SITE,
                !is.na(lat), !is.na(lon), !is.na(datetime), !is.na(Height_m)) |>
  dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
fp <- .site_footprint_pixels(range(df1$lat), range(df1$lon), m$.spatial)
cat(sprintf("%s: recomputed footprint = %d px (run logged %d)\n", SITE, length(fp), EXPECTED))
stopifnot(length(fp) == EXPECTED)

file.copy(path, paste0(path, ".bak_pre_footprint"), overwrite = FALSE)
m$.footprint_idx <- fp
m$.site_name <- SITE
saveRDS(m, path)
cat(sprintf("patched %s (backup: %s.bak_pre_footprint)\n", path, path))
