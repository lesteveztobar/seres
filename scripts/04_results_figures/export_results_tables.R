# export_results_tables.R -- long-format CSV exports behind the Results
# figures (brief 2026-10-02, Phase 1). Reads v9-era outputs only; writes only
# to output/tables/exports/. Never touches output/sensitivity_runs_v10/.
#
#   vertical_profiles.csv               1.1  footprint basis (same series as A09)
#   vertical_profiles_hourly_rasterpool.csv  1.1 companion, true hourly basis
#   heldout_modelled_heights.csv        1.2  site, rep, height_m, stage, n_individuals
#   heldout_observed_heights.csv        1.2  site, species, observed_height
#   persistence_trajectories.csv        1.3
#   height_resolution_runs.csv          1.4
#   niche_curves.csv / niche_breadth.csv / niche_records.csv   1.5
#   site_elevation.csv                  (support for the figures)
#
# Section 1.1 provenance: per-pixel HOURLY series do not exist (discarded at
# write time; only per-pixel 5/25/50/75/95 quantiles survive). The footprint
# file therefore aggregates each footprint pixel's annual per-height MEDIAN
# (the .build_site_climate_series_footprint() series A09's slopes use) across
# pixels: p10/p90 are SPATIAL, n_hours is NA, n_pixels is given. The hourly
# companion uses site_climate_cache_v2.rds (whole-raster spatial mean per
# hour): p10/p90 are TEMPORAL and n_hours is real. Shortwave is day-only
# (swdown > 0) in both, matching every per-pixel radiation consumer.
#
# Usage: sbatch scripts/04_results_figures/run_export_results_tables.sh

source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

EXP_DIR <- file.path(OUTPUT_DIR, "tables", "exports")
dir.create(EXP_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(!grepl("sensitivity_runs_v10", EXP_DIR))
put <- function(df, name) {
  write.csv(df, file.path(EXP_DIR, name), row.names = FALSE)
  cat(sprintf(">> %s: %d rows\n", name, nrow(df)))
}
SITES7 <- c("Mashpi", "Saloya", "Maquipucuna", "MindoTarabita", "MindoMirador", "Yanayacu", "LaElenita")
VAR_LABEL <- c(temp = "temperature", relhum = "relhum", swdown = "shortwave")
RUN_SITES <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")

cat("Observations:", OBSERVATIONS_CSV, "\nNiche cache: ", NICHE_CACHE_PATH, "\n\n")

# ── site elevation (one source for all seven) ────────────────────────────────
sc <- readRDS(file.path(PROCESSED_DIR, "site_climate_cache_v2.rds"))
elev <- data.frame(site = names(sc$site_elev), elevation_m = unname(sc$site_elev),
                   source = "site_climate_cache_v2.rds$site_elev (mean observation elevation)")
fe <- read.csv(file.path(OUTPUT_DIR, "tables", "climate", "footprint_elevation.csv"))
elev$footprint_elev_mean_m <- fe$elev_mean[match(elev$site, fe$site)]
put(elev[order(elev$elevation_m), ], "site_elevation.csv")

# ── 1.1 vertical profiles ────────────────────────────────────────────────────
cat("\n== 1.1 vertical profiles ==\n")
niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
niches <- augment_elevation(niches)                      # identical to A09
fc <- .build_site_climate_series_footprint(SITES7, niches)
cat("Footprint sites built:", paste(names(fc$pixels), collapse = ", "), "\n")
prof <- do.call(rbind, lapply(names(fc$pixels), function(s) {
  px <- fc$pixels[[s]]
  do.call(rbind, lapply(names(VAR_LABEL), function(v) {
    do.call(rbind, lapply(split(px, px$height), function(d) {
      y <- d[[v]]; y <- y[is.finite(y)]
      data.frame(site = s, height_m = d$height[1], relative_height = d$rel_height[1],
                 variable = VAR_LABEL[[v]], mean = mean(y), median = median(y),
                 p10 = unname(quantile(y, 0.10)), p90 = unname(quantile(y, 0.90)),
                 n_hours = NA_integer_, n_pixels = length(y), n_pixels_na = nrow(d) - length(y),
                 basis = "footprint pixels, per-pixel annual median; p10/p90 across pixels")
    }))
  }))
}))
put(prof, "vertical_profiles.csv")

ceil_hourly <- sapply(names(sc$site_pixels), function(s)
  site_canopy_ceiling(s, niches, max(sc$site_pixels[[s]]$height)))
prof_h <- do.call(rbind, lapply(names(sc$site_pixels), function(s) {
  px <- sc$site_pixels[[s]]
  do.call(rbind, lapply(names(VAR_LABEL), function(v) {
    do.call(rbind, lapply(split(px, px$height), function(d) {
      y <- d[[v]]
      if (v == "swdown") y <- y[y > 0]
      y <- y[is.finite(y)]
      data.frame(site = s, height_m = d$height[1], relative_height = d$height[1] / ceil_hourly[[s]],
                 variable = VAR_LABEL[[v]], mean = mean(y), median = median(y),
                 p10 = unname(quantile(y, 0.10)), p90 = unname(quantile(y, 0.90)),
                 n_hours = length(y), n_pixels = NA_integer_, n_pixels_na = NA_integer_,
                 basis = "whole-raster spatial mean per hour; p10/p90 across hours; shortwave swdown>0 hours")
    }))
  }))
}))
put(prof_h, "vertical_profiles_hourly_rasterpool.csv")

# ── 1.5 niche (before the big RDS loop so a crash there still leaves these) ──
cat("\n== 1.5 niche ==\n")
obs <- load_observations()
held_out <- get_held_out_split(obs)
n_train <- table(obs$FinalID[!held_out])
n_total <- table(obs$FinalID)
cache <- readRDS(NICHE_CACHE_PATH)
curves <- do.call(rbind, lapply(names(cache), function(sp) do.call(rbind, lapply(names(cache[[sp]]$axes), function(ax) {
  a <- cache[[sp]]$axes[[ax]]
  data.frame(taxon_unit = sp, axis = ax, x = a$x, score = a$score)
}))))
curves$n_train <- as.integer(n_train[curves$taxon_unit])
put(curves, "niche_curves.csv")
old <- read.csv(file.path(OUTPUT_DIR, "tables", "niche", "niche_curves.csv"))
cat(sprintf("Check vs tables/niche/niche_curves.csv: same rows %s, max |score diff| %.3g\n",
            nrow(old) == nrow(curves), max(abs(old$score - curves$score), na.rm = TRUE)))
grids_identical <- sapply(unique(curves$axis), function(ax) {
  xs <- lapply(split(curves$x[curves$axis == ax], curves$taxon_unit[curves$axis == ax]), range)
  length(unique(xs)) == 1
})
cat("x grid range identical across taxa per axis:", paste(names(grids_identical), grids_identical, collapse = ", "), "\n")
breadth <- do.call(rbind, lapply(split(curves, list(curves$taxon_unit, curves$axis), drop = TRUE), function(d) {
  d <- d[order(d$x), ]; x <- d$x; n <- length(x)
  # cell width per grid point (midpoint rule), so unequal spacing is handled
  w <- diff(c(x[1], (x[-1] + x[-n]) / 2, x[n]))
  data.frame(taxon_unit = d$taxon_unit[1], axis = d$axis[1], n_train = d$n_train[1],
             n_records_total = as.integer(n_total[d$taxon_unit[1]]),
             breadth_frac = sum(w[is.finite(d$score) & d$score > 50]) / (x[n] - x[1]),
             x_min = x[1], x_max = x[n])
}))
put(breadth, "niche_breadth.csv")
rec <- as.data.frame(table(taxon_unit = obs$FinalID, site = obs$Area_or_Site), stringsAsFactors = FALSE)
names(rec)[3] <- "n_records"
put(rec[rec$n_records > 0, ], "niche_records.csv")
cat("Taxa in cache without records:", paste(setdiff(names(cache), names(n_total)), collapse = ", "), "\n")

# ── 1.4 height resolution ────────────────────────────────────────────────────
cat("\n== 1.4 height resolution ==\n")
hr_sites <- c("Maquipucuna", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")
hr_fin <- do.call(rbind, lapply(hr_sites, function(s) {
  hr_file <- file.path(PROCESSED_DIR, sprintf("height_resolution_%s_realistic_75founders.rds", s))
  if (!file.exists(hr_file)) return(NULL)
  hr <- readRDS(hr_file)
  fin <- do.call(rbind, lapply(split(hr, list(hr$height_step, hr$rep), drop = TRUE), function(d) {
    d <- d[d$t == max(d$t), ]
    data.frame(site = s, spacing_m = d$height_step, n_tiers = d$n_heights, rep = d$rep,
               t_final = d$t, final_total = d$total, recruited = d$recruited, regime = d$regime)
  }))
  fin[order(fin$spacing_m, fin$rep), ]
}))
put(hr_fin, "height_resolution_runs.csv")
print(aggregate(final_total ~ site + spacing_m + n_tiers, hr_fin, function(x) c(mean = round(mean(x)), sd = round(sd(x)), n = length(x))))

# ERA5 forcing cells are NOT derived here: matching each site's .weather
# series to its ERA5 grid by correlation gave r = 0.09-0.35 and no shared
# series (.weather is site-adjusted), so it identifies nothing. The map reads
# output/tables/audits/home_cell_validity_lookup.csv instead.

# ── 1.2 + 1.3 from the realistic runs (one file in memory at a time) ─────────
cat("\n== 1.2 held-out heights / 1.3 persistence ==\n")
scores <- read.csv(file.path(OUTPUT_DIR, "tables", "validation", "held_out_validation_scores_realistic.csv"))
first_rep <- tapply(scores$rep, scores$site, min)
obs_h <- scores[scores$rep == first_rep[scores$site], c("site", "species", "observed_height")]
put(obs_h, "heldout_observed_heights.csv")
realized <- read.csv(file.path(OUTPUT_DIR, "tables", "validation", "held_out_realized_heights_realistic.csv"))

traj <- list(); mod <- list()
for (s in RUN_SITES) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_realistic_h0.40.rds", s))
  csv <- read.csv(file.path(OUTPUT_DIR, "tables", "run_summaries", sprintf("colonization_%s_realistic_h0.40_summary.csv", s)))
  res <- readRDS(f)
  cols <- names(csv)
  same <- isTRUE(all.equal(as.data.frame(res$summary[, cols]), csv, check.attributes = FALSE))
  cat(sprintf("  %s: run_summaries CSV identical to $summary: %s\n", s, same))
  if (!same) warning(s, ": run_summaries CSV differs from RDS $summary -- exporting RDS")
  traj[[s]] <- data.frame(site = s, experiment = "realistic", res$summary[, c("rep", "t", "totalS", "totalJ", "totalA",
                                                                            "total", "recruited", "extinct", "regime")])
  if (s %in% scores$site) {
    heights <- microenv_heights(readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", s))))
    for (rep_i in seq_along(res$runs)) {
      run <- res$runs[[rep_i]]
      if (is.null(run) || is.null(run$obs_val) || nrow(run$obs_val) == 0) next
      final_t <- dim(run$abundanceA)[4]
      for (st in c("S", "J", "A")) {
        sl <- run[[paste0("abundance", st)]][, , , final_t, , drop = FALSE]
        occ <- which(sl > 0, arr.ind = TRUE)
        if (!nrow(occ)) next
        agg <- tapply(sl[occ], heights[occ[, 3]], sum)
        mod[[length(mod) + 1]] <- data.frame(site = s, rep = rep_i, height_m = as.numeric(names(agg)),
                                             stage = st, n_individuals = as.integer(agg), t_final = final_t)
      }
    }
  }
  rm(res); invisible(gc())
}
traj <- do.call(rbind, traj)
put(traj, "persistence_trajectories.csv")
mod <- do.call(rbind, mod)
put(mod, "heldout_modelled_heights.csv")

# check: summed over stages, must reproduce held_out_realized_heights_realistic.csv exactly
a <- aggregate(n_individuals ~ site + rep + height_m, mod, sum)
b <- aggregate(n ~ site + rep + modelled_height, transform(realized, n = 1L), sum)
names(b) <- c("site", "rep", "height_m", "n_ref")
ab <- merge(a, b, all = TRUE)
cat(sprintf("Stage-summed vs held_out_realized_heights_realistic.csv: %d keys, %d mismatches, totals %d vs %d\n",
            nrow(ab), sum(is.na(ab$n_individuals) | is.na(ab$n_ref) | abs(ab$n_individuals - ab$n_ref) > 1e-9),
            sum(mod$n_individuals), nrow(realized)))
cat("\nDone.\n")
