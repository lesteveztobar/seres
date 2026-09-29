# plot_niche_per_site.R — regenerates the per-site niche diagnostics only:
#
#   niche_suitability_<site>.png          plot_niche_suitability()
#   niche_profile_curves_<site>.png       plot_niche_profile_curves()
#   niche_suitability_heatmap_<site>.png  plot_niche_suitability_heatmap()
#
# Lifted out of plot_all.R's niche loop so they can be rebuilt against a new
# observations dataset / niche cache without re-running the whole plot_all
# battery. The expensive part is one .niche_plot_context() (fresh climate
# cache) per site, shared across the three figures.
#
# Honours CANOPY_OBS_CSV / CANOPY_NICHE_CACHE / CANOPY_NICHE_BACKGROUND
# (paths.R) -- see run_niche_per_site.sh for the v6 invocation.
#
# Usage: sbatch scripts/experiments/A17_niche_suitability_plots/run_niche_per_site.sh [site ...]
#        (no args -> every site in the active observations CSV)

source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")

.plot_niche_per_site <- function(sites = NULL) {
  safe_plot <- function(label, expr)
    tryCatch(expr, error = function(e)
      message("SKIPPED (error) -- ", label, ": ", conditionMessage(e)))

  if (is.null(sites) || length(sites) == 0) {
    obs <- load_observations()
    sites <- sort(unique(obs$Area_or_Site[!is.na(obs$Area_or_Site) &
                                          nzchar(obs$Area_or_Site)]))
  }
  cat("Observations:", OBSERVATIONS_CSV, "\n")
  cat("Niche cache: ", NICHE_CACHE_PATH, "\n")
  cat("Sites:", paste(sites, collapse = ", "), "\n\n")

  for (site in sites) {
    cat(sprintf("=== %s ===\n", site))
    ctx <- tryCatch(.niche_plot_context(site), error = function(e) {
      message("SKIPPED (error) -- niche context ", site, ": ", conditionMessage(e)); NULL
    })
    if (is.null(ctx)) next
    safe_plot(paste("niche suitability", site),
      plot_niche_suitability(site, context = ctx))
    safe_plot(paste("niche profile curves", site),
      plot_niche_profile_curves(site, context = ctx))
    safe_plot(paste("niche suitability heatmap", site),
      plot_niche_suitability_heatmap(site, context = ctx))
  }
}

if (sys.nframe() == 0L) {
  .plot_niche_per_site(commandArgs(trailingOnly = TRUE))
  cat("\nDone. Figures (where inputs existed) are in", OUTPUT_DIR, "\n")
}
