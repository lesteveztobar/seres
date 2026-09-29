# pooled_fallback_audit.R -- per site, per variable: what fraction of the
# model's ACTUAL 3D voxel grid (every landscape (x,y,z) triple, mapped
# through clim_pixel_row/clim_pixel_col exactly as .clim_voxel_slice() does)
# carries real per-pixel climate vs. falls back to the pooled site mean via
# .fill_na(). Builds state via init_colonization() only -- no spinup, no
# simulation, no replicates -- a query, not a run.
suppressMessages({
  library(plotly); library(ggplot2); library(patchwork); library(parallel)
})
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("pooled_fallback_audit_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))))

sites <- c("Maquipucuna","Mashpi","Yanayacu","MindoMirador","MindoTarabita","Saloya","LaElenita")
# Real consumer keys: temp/relhum/swdown are read "day" (run_pass2_establish()/
# run_pass3_survive_grow() both key day-only for these); windspeed is "both"
# (run_pass1_disperse()'s wind_by_height). Annual bucket used here as the
# representative summary (matches Pass1/Pass2's own annual-resolution reads;
# Pass3's monthly reads are a finer partition of the same underlying data).
VAR_DAYPART <- c(temp = "day", relhum = "day", swdown = "day", windspeed = "both")

niches_all <- load_observations()
niches_all <- niches_all[!is.na(niches_all$lat) & !is.na(niches_all$lon) &
                         !is.na(niches_all$Height_m) & !is.na(niches_all$FinalID), ]

rows <- list()
pixel_status <- list()

for (s in sites) {
  # 2026-09-08: pointed at the archived (pre-grid-snap) manifests, not the
  # live path -- the live path is mid-repopulation by the grid-snapped
  # Phase A re-run right now, and racing that (as happened on the first
  # attempt at this job, which crashed when the underlying scratch height
  # files were moved out from under it) is not safe. Maquipucuna/Mashpi/
  # Yanayacu's domains are UNCHANGED by the grid-snap (shift=0), so their
  # numbers here are already final; MindoMirador/MindoTarabita/Saloya/
  # LaElenita's domains shifted, so THEIR numbers here reflect the OLD
  # (pre-snap) domain and should be re-checked once the new run lands.
  microenv_path <- file.path(PROCESSED_DIR, "archive_pre_gridsnap", sprintf("microenv_%s_h0.40.rds", s))
  if (!file.exists(microenv_path)) { message("Skipping ", s, " -- no microenv"); next }
  microenv <- readRDS(microenv_path)
  # .height_dir inside the manifest still names the LIVE scratch path (now
  # mid-repopulation) -- repoint it at the archived scratch copy, same fix
  # as applied to the B4/negative-control archived manifests earlier.
  microenv$.height_dir <- sub(
    "s38leste_hpc-seres/microenv_", "s38leste_hpc-seres/archive_pre_gridsnap/microenv_",
    microenv$.height_dir, fixed = TRUE)
  if (is.null(microenv$.weather)) {
    pm <- readRDS(file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", s)))
    microenv$.weather <- pm[[1]]$weather
  }
  available_heights <- microclimate_heights <- microenv_heights(microenv)
  canopy_ceiling <- site_canopy_ceiling(s, niches_all, max(available_heights))
  canopy_grid <- matrix(canopy_ceiling, nrow = 50, ncol = 50)
  site_obj <- list(Site = s)
  forestparams <- site_forestparams(s, canopy_ceiling)
  held_out <- get_held_out_split(niches_all[niches_all$Area_or_Site == s, ], train_frac = 0.70)
  site_obs <- niches_all[niches_all$Area_or_Site == s, ]
  niches_train <- site_obs[!held_out, ]

  message("Building state for ", s, "...")
  state <- tryCatch(
    init_colonization(site_obj, niches_train, canopy_grid, microenv,
      resolution = 10, carCap = 5, maxDisp = 10,
      params = list(canopy_z = canopy_ceiling, n_founders = 30),
      forestparams = forestparams),
    error = function(e) { message("Skipping ", s, " -- ", conditionMessage(e)); NULL })
  if (is.null(state)) next

  zDim <- state$zDim
  xDim <- state$xDim
  yDim <- state$yDim
  cache <- state$clim_cache_voxel
  pixel_idx_grid <- (as.vector(state$clim_pixel_col) - 1L) * state$clim_cache_voxel$clim_voxel_by_height[[1]]$nr +
                     as.vector(state$clim_pixel_row)

  # 2026-09-09 FIX: the original version only checked the "annual" bucket.
  # An annual quantile pools all 12 months' worth of hours together, so a
  # pixel with genuinely NA data in one specific month can still show a
  # finite annual value (the other 11 months carry it) -- annual-only
  # checking cannot see a per-month gap. Now checks all 12 months
  # individually AND annual, matching exactly what the two real consumers
  # read: run_pass2_establish()/dispersal use "annual" (Pass1/2),
  # run_pass3_survive_grow() uses the SPECIFIC calendar month (monthly
  # survival/growth) -- both need auditing, not just one.
  for (v in names(VAR_DAYPART)) {
    daypart <- VAR_DAYPART[[v]]
    for (month_label in c("annual", as.character(1:12))) {
      finite_count <- 0L
      total_count <- 0L
      na_pixel_rowcol <- list()
      for (zi in seq_len(zDim)) {
        height_key <- names(cache$clim_voxel_by_height)[zi]
        entry <- cache$clim_voxel_by_height[[height_key]]
        if (is.null(entry)) next
        key <- sprintf("%s_%s_%s", month_label, daypart, v)
        q_mat <- entry$quantiles[[key]]
        if (is.null(q_mat)) next
        pos <- match(pixel_idx_grid, entry$footprint_full_idx)
        # A voxel's value is "real" iff its pixel's quantile row has NO NA
        # among the 5 quantiles (.sample_quantile()'s own anyNA() gate).
        pixel_ok <- !is.na(pos) & apply(q_mat[pos, , drop = FALSE], 1, function(r) !anyNA(r))
        finite_count <- finite_count + sum(pixel_ok)
        total_count  <- total_count + length(pixel_ok)
        if (zi == 1) {
          bad <- which(!pixel_ok)
          if (length(bad) > 0) {
            na_pixel_rowcol[[v]] <- data.frame(
              row = as.vector(state$clim_pixel_row)[bad], col = as.vector(state$clim_pixel_col)[bad])
          }
        }
      }
      frac_real <- if (total_count > 0) finite_count / total_count else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        site = s, variable = v, month = month_label, xDim = xDim, yDim = yDim, zDim = zDim,
        n_voxel_var_combinations = total_count,
        frac_real_data = frac_real, frac_pooled_fallback = 1 - frac_real
      )
    }
  }

  # Spatial pattern of NA raster pixels at height tier 1, for temp (day):
  entry1 <- cache$clim_voxel_by_height[[names(cache$clim_voxel_by_height)[1]]]
  q1 <- entry1$quantiles[["annual_day_temp"]]
  all_pixel_idx <- entry1$footprint_full_idx
  pixel_ok1 <- apply(q1, 1, function(r) !anyNA(r))
  nr <- entry1$nr; nc <- entry1$nc
  rows_all <- ((all_pixel_idx - 1L) %% nr) + 1L
  cols_all <- ((all_pixel_idx - 1L) %/% nr) + 1L
  touched <- unique(data.frame(row = as.vector(state$clim_pixel_row), col = as.vector(state$clim_pixel_col)))
  touched$pos <- match((touched$col - 1L) * nr + touched$row, all_pixel_idx)
  touched$valid <- pixel_ok1[touched$pos]
  pixel_status[[s]] <- touched
  message(sprintf("  %s: %d distinct raster pixels touched by the landscape, %d valid, %d NA",
                  s, nrow(touched), sum(touched$valid), sum(!touched$valid)))
  print(touched)
}

out <- do.call(rbind, rows)
write.csv(out, "output/pooled_fallback_audit.csv", row.names = FALSE)
cat("\n=== POOLED FALLBACK FRACTION, PER SITE PER VARIABLE (whole 3D voxel grid) ===\n")
print(out, row.names = FALSE)

saveRDS(pixel_status, "output/pooled_fallback_pixel_status.rds")
