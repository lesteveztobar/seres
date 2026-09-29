# make_isolation_species_files.R
# Builds one params$species_subset RDS per (site, species) pair -- see
# run_colonization.R's species_file arg / init_colonization()'s
# params$species_subset (get_colonization.R) -- so each species can be run
# alone (no competitors) and compared against the existing multi-species run
# for the same site, for the competition/microhabitat-preference isolation
# test.
#
# One file per (site, species) pair, not per species -- a species observed
# at two sites (e.g. Maxillaria bradei at both Mashpi and MindoTarabita)
# needs a separate isolation run per site, since each site is a different
# physical landscape.
#
# Also writes isolation_manifest.csv (site, species, species_file, exp_tag)
# so the orchestration script and the downstream competition-analysis script
# both iterate over exactly the same list, rather than re-deriving it
# independently and risking drift.
#
# Usage: Rscript scripts/experiments/A10_competition_isolation/make_isolation_species_files.R
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")  # .filter_maxillariinae() -- 2026-08-29

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]

# Filesystem-safe tag from a site/species string: non-alphanumeric runs
# collapse to one underscore, trimmed of leading/trailing underscores.
safe_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

pairs <- unique(niches[, c("Area_or_Site", "FinalID")])
pairs <- pairs[order(pairs$Area_or_Site, pairs$FinalID), ]

manifest <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  site <- pairs$Area_or_Site[i]
  sp   <- pairs$FinalID[i]
  sp_tag  <- safe_tag(sp)
  fname   <- sprintf("isolation_%s_%s.rds", safe_tag(site), sp_tag)
  fpath   <- file.path(PARAMS_DIR, fname)
  saveRDS(sp, fpath)  # params$species_subset expects a character vector
  data.frame(site = site, species = sp, species_file = fpath,
            exp_tag = sprintf("isolation_%s", sp_tag), stringsAsFactors = FALSE)
}))

# 2026-09-29: a second run of this script (different baseline/species set)
# used to silently overwrite the manifest that a prior competition_analysis.R
# run's results depended on, with no record of what changed. Fixed two ways,
# without changing the fixed read path competition_analysis.R relies on:
# (1) an existing manifest is archived (not lost) before being replaced;
# (2) every row of the new manifest is stamped with manifest_generation, so
# any consumer -- or a human diffing two manifests later -- can tell which
# generation of isolation-file assignments a given competition_analysis.R
# run actually used.
manifest_path <- file.path(PARAMS_DIR, "isolation_manifest.csv")
manifest_generation <- format(Sys.time(), "%Y%m%d_%H%M%S")
if (file.exists(manifest_path)) {
  archive_dir <- file.path(PARAMS_DIR, "isolation_manifest_archive")
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  old_stamp <- format(file.info(manifest_path)$mtime, "%Y%m%d_%H%M%S")
  archive_path <- file.path(archive_dir, sprintf("isolation_manifest_%s.csv", old_stamp))
  file.copy(manifest_path, archive_path, overwrite = FALSE)
  warning(sprintf(
    "isolation_manifest.csv already existed (from %s) -- archived to %s before overwriting. Any competition_analysis.R run against the old manifest is no longer reproducible from the live path, but its exact species/site assignments are preserved there.",
    old_stamp, archive_path
  ))
}
manifest$manifest_generation <- manifest_generation
write.csv(manifest, manifest_path, row.names = FALSE)

cat(sprintf("Wrote %d isolation species_subset files (manifest_generation %s):\n",
            nrow(manifest), manifest_generation))
print(manifest[, c("site", "species", "exp_tag")], row.names = FALSE)
cat(sprintf("\nManifest saved: %s\n", manifest_path))
