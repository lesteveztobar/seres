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
# Usage: Rscript scripts/02_model/setup/make_isolation_species_files.R
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")

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

manifest_path <- file.path(PARAMS_DIR, "isolation_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)

cat(sprintf("Wrote %d isolation species_subset files:\n", nrow(manifest)))
print(manifest[, c("site", "species", "exp_tag")], row.names = FALSE)
cat(sprintf("\nManifest saved: %s\n", manifest_path))
