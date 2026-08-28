# competition_analysis.R
# Isolation-run competition test: for each (site, species) pair, compares
# that species' realized height distribution when run ALONE (no
# competitors -- data/params/isolation_<site>_<species>.rds, see
# make_isolation_species_files.R) against the SAME species' realized height
# distribution in the existing multi-species run for that site. Both
# conditions use the same params (only species_subset differs) -- so any
# shift in the height distribution isolates the effect of the other species
# being present (competition/facilitation), separate from the species'
# underlying niche preference (which the isolation run alone reveals).
#
# Per replicate, "realized height distribution" = every individual's height
# tier at the final timestep (S+J+A stages pooled, weighted by abundance --
# i.e. a height tier with 5 individuals contributes 5 copies of that height
# to the sample), pooled across all n_reps replicates. Compared via a
# two-sample Kolmogorov-Smirnov test (distribution shape/location, not just
# the mean) -- a significant result with the isolation sample's mean height
# on one side means the species is displaced toward the other side of its
# preferred range when competitors are present.
#
# Requires: make_isolation_species_files.R already run (isolation_manifest.csv
# + the isolation params files), the isolation colonization runs themselves
# already completed with the SAME baseline_tag's params, and each site's
# multi-species run at that same tag already completed (skips any
# site/species without both).
#
# Usage: Rscript scripts/02_model/analysis/competition_analysis.R [baseline_tag]
#   baseline_tag: exp_tag of the multi-species run to compare against
#   (colonization_<site>_<baseline_tag>_h0.40.rds). Defaults to
#   "realistic_273founders" -- realistic.rds's literature-default 30
#   founders produces populations too small/noisy (final abundance ~4-8,
#   DECLINING from the 30 founders) for this comparison to have any power;
#   see run_full_analysis_pipeline.sh's header note on why 273 founders
#   specifically. Must match whatever exp_tag the isolation runs (step 5,
#   run_full_analysis_pipeline.sh) and the multi-species run (step 4) were
#   actually submitted with.
# Output: output/competition_analysis.csv (one row per site x species)
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

args         <- commandArgs(trailingOnly = TRUE)
baseline_tag <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "realistic_273founders"

manifest_path <- file.path(PARAMS_DIR, "isolation_manifest.csv")
if (!file.exists(manifest_path))
  stop("No isolation manifest at ", manifest_path, " -- run make_isolation_species_files.R first")
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

# Every individual's height (S+J+A pooled) at the final timestep, one
# replicate's worth, repeated by abundance -- e.g. 3 adults at 8.5m
# contributes c(8.5, 8.5, 8.5) to the returned vector. sp_idx: which
# species_ids index this species is in THIS run (differs between an
# isolation run, where it's always 1, and a multi-species run).
.realized_heights <- function(run, sp_idx) {
  if (is.null(run)) return(numeric(0))
  T <- dim(run$abundanceA)[4]
  slice <- run$abundanceS[,,,T,sp_idx] + run$abundanceJ[,,,T,sp_idx] + run$abundanceA[,,,T,sp_idx]
  # slice: [xDim, yDim, zDim] counts -- collapse to total count per height (z)
  per_height <- apply(slice, 3, sum)
  rep(run$heights, times = per_height)
}

# Pools realized heights across every successful replicate in a
# run_replicated()-style result (list(runs=..., summary=...)).
.pooled_heights <- function(result, species_name) {
  if (is.null(result) || is.null(result$runs)) return(NULL)
  runs <- Filter(Negate(is.null), result$runs)
  if (length(runs) == 0) return(NULL)
  unlist(lapply(runs, function(r) {
    sp_idx <- match(species_name, r$species_ids)
    if (is.na(sp_idx)) return(numeric(0))
    .realized_heights(r, sp_idx)
  }))
}

rows <- list()
for (i in seq_len(nrow(manifest))) {
  site <- manifest$site[i]; sp <- manifest$species[i]; exp_tag <- manifest$exp_tag[i]

  iso_path   <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s_h0.40.rds", site, exp_tag))
  multi_path <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s_h0.40.rds", site, baseline_tag))
  if (!file.exists(iso_path) || !file.exists(multi_path)) {
    message("Skipping ", site, " / ", sp, " -- missing ",
            if (!file.exists(iso_path)) iso_path else multi_path)
    next
  }

  iso_heights   <- .pooled_heights(readRDS(iso_path), sp)
  multi_heights <- .pooled_heights(readRDS(multi_path), sp)

  if (is.null(iso_heights) || is.null(multi_heights) ||
      length(iso_heights) == 0 || length(multi_heights) == 0) {
    rows[[length(rows) + 1]] <- data.frame(
      site = site, species = sp,
      n_isolation = length(iso_heights), n_multispecies = length(multi_heights),
      mean_height_isolation = NA, mean_height_multispecies = NA,
      ks_D = NA, ks_p = NA,
      note = "no individuals present at final timestep in one or both runs"
    )
    next
  }

  kt <- tryCatch(suppressWarnings(ks.test(iso_heights, multi_heights)), error = function(e) NULL)
  rows[[length(rows) + 1]] <- data.frame(
    site = site, species = sp,
    n_isolation = length(iso_heights), n_multispecies = length(multi_heights),
    mean_height_isolation = mean(iso_heights), mean_height_multispecies = mean(multi_heights),
    ks_D = if (!is.null(kt)) unname(kt$statistic) else NA,
    ks_p = if (!is.null(kt)) kt$p.value else NA,
    note = ""
  )
}

result_df <- do.call(rbind, rows)
if (is.null(result_df) || nrow(result_df) == 0) stop("No site/species pairs had both runs available.")

result_df$height_shift_m <- result_df$mean_height_multispecies - result_df$mean_height_isolation
result_df$sig <- .sig_stars(result_df$ks_p, na_str = "")  # shared_helpers.R

cat("\n== Competition test: realized height distribution, isolation vs. multi-species ==\n")
cat("height_shift_m = mean height WITH competitors minus mean height ALONE\n")
cat("  (positive = species sits higher with competitors present; negative = lower)\n\n")
print(result_df[, c("site", "species", "n_isolation", "n_multispecies",
                    "mean_height_isolation", "mean_height_multispecies",
                    "height_shift_m", "ks_D", "ks_p", "sig")], row.names = FALSE)

out_path <- file.path(OUTPUT_DIR, "competition_analysis.csv")
write.csv(result_df, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_path))
