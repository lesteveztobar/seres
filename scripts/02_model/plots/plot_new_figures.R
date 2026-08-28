# plot_new_figures.R — orchestrates the three new figure types added
# 2026-08-24 (see plots/plot_functions.R): per-species/per-stage abundance
# curves, the niche-suitability hotbox, and the factorial survival-rate
# heatmap. Separate from plot_all.R (rather than folding in) so it can be
# rerun independently as more sites finish their sweep, without re-running
# plot_all.R's whole existing battery of figures.
#
# Usage: Rscript scripts/02_model/plots/plot_new_figures.R [site1] [site2] ...
#   Defaults to every site with a best_combo colonization result on disk
#   (i.e. every site pick_best_combo.R + the representative run have already
#   completed for -- see scripts/02_model/analysis/pick_best_combo.R).
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  sites <- args
} else {
  files <- list.files(PROCESSED_DIR, pattern = "^colonization_.*_best_combo_h0\\.40\\.rds$")
  sites <- sub("^colonization_(.*)_best_combo_h0\\.40\\.rds$", "\\1", files)
}
cat("Sites:", paste(sites, collapse = ", "), "\n")

safe_plot <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    message("SKIPPED (error) -- ", label, ": ", conditionMessage(e))
  })
}

for (site in sites) {
  cat(sprintf("\n=== %s ===\n", site))

  safe_plot(paste("species/stage curves", site),
    plot_species_stage_curves(site, exp_tag = "best_combo"))

  # 2026-08-28: swapped the per-site suitability hotbox for the combined
  # all-sites grid + tree figures (below, once per run rather than once per
  # site) -- plot_niche_suitability_heatmap() stays available/callable on
  # its own for per-site debugging, just not part of this default batch.

  # 2026-08-25: swapped for the pairwise parameter-effect heatmap (answers
  # "which parameter(s) actually move abundance", per-site) -- the old
  # extinction/survival tile+facet view (plot_factorial_experiment()) stays
  # available/callable on its own, just not part of this default batch.
  safe_plot(paste("factorial pairwise heatmap", site),
    plot_factorial_pairwise_heatmap(site, exp_tag = "reproduction_factorial_v3_h0.40"))
}

# 2026-08-28: combined, all-sites niche-suitability figures -- one grid
# (site x height per species) and one tree-bar version, replacing the old
# per-site suitability hotbox call inside the loop above. Both build their
# own site list from whatever microenv_<site>_h0.40.rds files exist, so no
# `sites` arg needed here. Build the shared (site, species, height,
# before/after) table ONCE and pass it to both -- each figure previously
# rebuilt every site's climate cache independently (14 builds across 7
# sites instead of 7), slow enough on the login node to be killed partway
# through.
niche_built <- .niche_score_rows_all_sites()
safe_plot("niche suitability grid (all sites)", plot_niche_suitability_grid(built = niche_built))
safe_plot("niche height trees (all sites)", plot_niche_height_trees(built = niche_built))

# 2026-08-25: cross-site species comparison -- every species genuinely
# identified (see .confirmed_species_sites(), shared_helpers.R) at 2+
# sites, derived from the observations table rather than hardcoded so a
# newly added site/observation is picked up automatically.
niches <- load_observations()
all_species <- sort(unique(niches$Identification[!is.na(niches$Identification) &
                                                  nzchar(trimws(niches$Identification))]))
multi_site_species <- Filter(function(sp) length(.confirmed_species_sites(sp, niches)) >= 2, all_species)
cat("\nSpecies confirmed at 2+ sites:", paste(multi_site_species, collapse = ", "), "\n")

for (sp in multi_site_species) {
  safe_plot(paste("cross-site comparison", sp), plot_species_across_sites(sp))
}

cat("\nAll done. Figures (where inputs existed) are in", OUTPUT_DIR, "\n")
