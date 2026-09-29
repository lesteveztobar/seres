# plot_niche_all_sites.R — regenerates ONLY the combined all-sites
# niche-suitability figure family (see plot_functions.R's "Combined
# all-sites niche-suitability figures" block):
#
#   niche_suitability_by_site.png               Fig A  (main, quantitative)
#   niche_forest_all_sites.png                  Fig B  (visual abstract)
#   niche_suitability_by_site_rescale.png       supp   (ceiling-rescale effect)
#   niche_suitability_by_site_morphospecies.png supp   (unidentified sp.)
#   niche_multisite_species.png                 supp   (species at 2+ sites)
#
# Split out of plot_new_figures.R (which also sources this, so the two never
# drift) so the niche family can be rebuilt on its own — its .niche_score_
# rows_all_sites() step is a full per-site climate-cache pass and is the
# slow part of that script.
#
# Honours CANOPY_OBS_CSV / CANOPY_NICHE_CACHE / CANOPY_NICHE_BACKGROUND
# (paths.R). To build against the v6 dataset, export those before running —
# see run_niche_figures.sh.
#
# Usage: sbatch scripts/experiments/A17_niche_suitability_plots/run_niche_figures.sh
#        (or, on an interactive node: Rscript scripts/experiments/A17_niche_suitability_plots/plot_niche_all_sites.R)

source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")

.plot_niche_all_sites <- function() {
  safe_plot <- function(label, expr)
    tryCatch(expr, error = function(e)
      message("SKIPPED (error) -- ", label, ": ", conditionMessage(e)))

  cat("Observations:", OBSERVATIONS_CSV, "\n")
  cat("Niche cache: ", NICHE_CACHE_PATH, "\n\n")

  # Build the shared (site, species, height, before/after) table ONCE
  # (cached to data/processed/, keyed by the active CSV + niche cache).
  niche_built <- .niche_score_rows_all_sites_cached()

  safe_plot("niche suitability by site (Fig A)",
    plot_niche_suitability_by_site(built = niche_built))
  safe_plot("niche forest (Fig B)",
    plot_niche_forest(built = niche_built))
  safe_plot("niche rescale comparison (supp)",
    plot_niche_suitability_by_site(built = niche_built, stages = c("before", "after")))
  safe_plot("niche morphospecies (supp)",
    plot_niche_suitability_by_site(built = niche_built, confirmed_only = FALSE))
  safe_plot("niche multi-site species (supp)",
    plot_niche_multisite_species(built = niche_built))

  invisible(niche_built)
}

# Run when invoked as a script (Rscript / sbatch), not when sourced.
if (sys.nframe() == 0L) {
  .plot_niche_all_sites()
  cat("\nDone. Figures (where inputs existed) are in", OUTPUT_DIR, "\n")
}
