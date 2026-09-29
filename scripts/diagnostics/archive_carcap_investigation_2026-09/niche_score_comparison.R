# niche_score_comparison.R
# Check 2 (continued) -- item 6: does excluding Saloya's coordinate-outlier
# rows from the background pool change actual niche SCORES, not just the
# key/pixel counts already reported? Scores a common set of points (each
# shared species' own observation climate, pooled across sites, exactly as
# characterize_niches.R already computes it) under both the with-Saloya
# production cache (species_niches_v6.rds/niche_background_density_v6.rds)
# and the without-Saloya comparison cache (*_noSaloya.rds, built earlier via
# EXCLUDE_SALOYA=1), and reports the distribution of per-point score
# differences. Read-only -- writes only to output/, no cache is touched.
#
# Usage: Rscript scripts/diagnostics/archive_carcap_investigation_2026-09/niche_score_comparison.R
# Lizeth Estévez Tobar -- University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

VARS <- NICHE_VARS
HEIGHT_STEP <- 0.4
manifest_suffix <- "_h0.40"

niche_with    <- readRDS(file.path(PROCESSED_DIR, "species_niches_v6.rds"))
niche_without <- readRDS(file.path(PROCESSED_DIR, "species_niches_v6_noSaloya.rds"))

shared_species <- intersect(names(niche_with), names(niche_without))
message(sprintf("Shared species (present in both caches): %d of %d (with-Saloya) / %d (without)",
                length(shared_species), length(niche_with), length(niche_without)))

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
niches <- niches[niches$Area_or_Site != "Saloya", ]  # score only at points both caches could plausibly place a species -- Saloya's own points aren't in the noSaloya run's landscape set
sites <- sort(unique(niches$Area_or_Site))

obs_rows <- list()
for (s in sites) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", s, manifest_suffix))
  if (!file.exists(microenv_path)) { message("Skipping ", s, " -- no microenv"); next }
  microenv <- readRDS(microenv_path)
  obs_site <- niches[niches$Area_or_Site == s, ]
  cc_voxel <- build_clim_cache_voxel(microenv, footprint = NULL)
  obs <- data.frame(lon = obs_site$lon, lat = obs_site$lat, height = obs_site$Height_m)
  presence <- voxel_climate_table(cc_voxel, microenv, obs, VARS)
  obs_rows[[s]] <- cbind(species = rep(obs_site$FinalID, each = 12), as.data.frame(presence))
}
obs_df <- do.call(rbind, obs_rows)
obs_df$species <- as.character(obs_df$species)

results <- do.call(rbind, lapply(shared_species, function(sp) {
  rows <- obs_df[obs_df$species == sp, VARS, drop = FALSE]
  if (nrow(rows) == 0) return(NULL)
  score_with    <- apply(rows, 1, function(r) niche_raw_score(as.list(r), niche_with[[sp]]))
  score_without <- apply(rows, 1, function(r) niche_raw_score(as.list(r), niche_without[[sp]]))
  data.frame(species = sp, score_with = score_with, score_without = score_without,
             diff = score_without - score_with)
}))

message(sprintf("\nScored %d points across %d shared species.", nrow(results), length(shared_species)))
message(sprintf("Score difference (without - with Saloya): median=%.2f  IQR=[%.2f, %.2f]  range=[%.2f, %.2f]",
                median(results$diff), quantile(results$diff, 0.25), quantile(results$diff, 0.75),
                min(results$diff), max(results$diff)))
message(sprintf("Mean |difference| = %.2f", mean(abs(results$diff))))

per_species <- aggregate(diff ~ species, results, function(x) c(median = median(x), mean_abs = mean(abs(x))))
print(per_species)

dir.create("output", showWarnings = FALSE)
write.csv(results, "output/niche_score_comparison_with_without_saloya.csv", row.names = FALSE)
message("\nSaved output/niche_score_comparison_with_without_saloya.csv")
