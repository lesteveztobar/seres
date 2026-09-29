source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

cat("=========================================================\n")
cat("PART A: empirical test of n_species=0 array-indexing behavior\n")
cat("=========================================================\n")
# Mirrors run_pass1_disperse()'s `for (sp in 1:state$n_species)` against a
# zero-length species dimension, isolated from the rest of the model.
n_species <- 0L
abundanceA <- array(0L, dim = c(3, 3, 3, 5, n_species))
cat("dim(abundanceA):", paste(dim(abundanceA), collapse=","), "\n")
result <- tryCatch({
  for (sp in 1:n_species) {
    x <- abundanceA[, , , 1, sp]
    cat("  sp=", sp, " -> extracted length ", length(x), "\n", sep="")
  }
  "NO ERROR -- loop completed silently"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("1:n_species with n_species=0 result:", result, "\n\n")

# Compare against seq_len(), used correctly elsewhere (run_spinup()'s founder loop)
result2 <- tryCatch({
  n <- 0
  for (sp in seq_len(n)) { stop("should never reach here") }
  "seq_len(0) loop correctly skipped (0 iterations)"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("seq_len(n_species) with n_species=0 result:", result2, "\n\n")

cat("=========================================================\n")
cat("PART B: per-site species-level record audit (current, post-filter data)\n")
cat("=========================================================\n")
niches <- load_observations()
all_sites <- c("Maquipucuna","Mashpi","MindoMirador","MindoTarabita","LaElenita","Saloya","Yanayacu")

is_morphospecies <- function(id) grepl("\\bsp\\.?\\s*[0-9A-Za-z]*$|\\bsp\\d*$", trimws(id), ignore.case = TRUE)

rows <- lapply(all_sites, function(s) {
  d <- niches[niches$Area_or_Site == s, ]
  n_records <- nrow(d)
  n_distinct_final_id <- length(unique(d$FinalID))
  sp_level <- d[!is_morphospecies(d$FinalID), ]
  n_sp_level_records <- nrow(sp_level)
  n_sp_level_species <- length(unique(sp_level$FinalID))
  sp_counts <- table(sp_level$FinalID)
  n_sp_ge3 <- sum(sp_counts >= 3)
  data.frame(site = s, n_records_total = n_records, n_distinct_FinalID = n_distinct_final_id,
             n_species_level_records = n_sp_level_records, n_species_level_species = n_sp_level_species,
             n_species_level_species_ge3records = n_sp_ge3)
})
audit_df <- do.call(rbind, rows)
print(audit_df, row.names = FALSE)

cat("\n--- FinalID values per site (post-filter) ---\n")
for (s in all_sites) {
  d <- niches[niches$Area_or_Site == s, ]
  cat(sprintf("%-15s (n=%d): %s\n", s, nrow(d), paste(sort(table(d$FinalID)), names(sort(table(d$FinalID))), sep="x", collapse=", ")))
}

cat("\n=========================================================\n")
cat("PART C: does init_colonization() actually crash for a zero-record site?\n")
cat("=========================================================\n")
for (s in c("MindoTarabita", "LaElenita", "Yanayacu")) {
  site_obs <- niches[niches$Area_or_Site == s, ]
  cat(sprintf("%-15s n_rows=%d\n", s, nrow(site_obs)))
  if (nrow(site_obs) == 0) {
    lat_range_m <- suppressWarnings((max(site_obs$lat) - min(site_obs$lat)) * 111000)
    lon_range_m <- suppressWarnings((max(site_obs$lon) - min(site_obs$lon)) * 111000)
    xDim <- max(round(lon_range_m / 10), 10) + 4
    yDim <- max(round(lat_range_m / 10), 10) + 4
    cat(sprintf("  lat_range_m=%s lon_range_m=%s -> xDim=%s yDim=%s (vs a real site's ~100+ cells)\n",
                lat_range_m, lon_range_m, xDim, yDim))
  }
}
