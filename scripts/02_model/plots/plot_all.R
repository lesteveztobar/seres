# plot_all.R
# Generates every project plot from whatever results currently exist under
# data/processed/ and geojson_to_csv/. Safe to re-run at any point in the
# pipeline — each function in plot_functions.R skips (with a message) if its
# inputs aren't there yet, instead of erroring.
#
# Usage: sh scripts/02_model/plots/run_plots.sh   (wraps: Rscript scripts/02_model/plots/plot_all.R)
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/plots/plot_functions.R")

# Sites derived from OBSERVATIONS_CSV itself, not hardcoded -- a newly added
# site's results get plotted automatically (every function below already
# skips gracefully with a message if that site's inputs aren't there yet).
.plot_all_niches <- load_observations()
SITES <- sort(unique(.plot_all_niches$Area_or_Site[
  !is.na(.plot_all_niches$Area_or_Site) & nzchar(.plot_all_niches$Area_or_Site)]))
EXP   <- names(EXP_PARAM_MAP)

# Wraps a single plotting call so one function's runtime error doesn't halt
# the rest of the script. plot_functions.R's own functions already skip
# gracefully when an input file is simply missing (see file header) -- this
# catches the separate case of a file that exists but errors while being
# read, e.g. plot_temperature_profile()'s stale terra/Rcpp external-pointer
# error on some saved microenv objects (2026-07-16) -- without this, that one
# site's failure previously took down every plot after it in the script,
# including sites/experiments that had nothing to do with the failure.
safe_plot <- function(label, expr) {
  tryCatch(
    expr,
    error = function(e) message("SKIPPED (error) ", label, ": ", conditionMessage(e))
  )
}

cat("== Field site map ==\n")
safe_plot("site map", plot_site_map())

cat("\n== Modelled tree diagram ==\n")
safe_plot("tree diagram", plot_tree_diagram())

cat("\n== Temperature profiles ==\n")
for (site in SITES) safe_plot(paste("temperature profile", site), plot_temperature_profile(site))

cat("\n== Niche suitability (per-axis + before/after normalization) ==\n")
# One .niche_plot_context() build per site (the expensive part -- a fresh
# climate cache), shared between plot_niche_suitability() and
# plot_niche_profile_curves() instead of each rebuilding its own.
for (site in SITES) {
  ctx <- safe_plot(paste("niche context", site), .niche_plot_context(site))
  if (is.null(ctx)) next
  safe_plot(paste("niche suitability", site), plot_niche_suitability(site, context = ctx))
  safe_plot(paste("niche profile curves", site), plot_niche_profile_curves(site, context = ctx))
}

cat("\n== Best-fit 3D comparison (simple model) ==\n")
safe_plot("bestfit 3D comparison", plot_bestfit_3d_comparison())

cat("\n== Colonization sensitivity experiments ==\n")
for (site in SITES) {
  for (tag in EXP) {
    safe_plot(paste("experiment", site, tag), plot_colonization_experiment(site, tag))
  }
}

cat("\n== Reproduction factorial (p_poll x p_germ x p_s1 x n_founders) ==\n")
# reproduction_factorial_v3: levels bracketed realistic.rds -> best_case.rds
# (see report/methods.tex, Full factorial experiment). Loops every site now
# that more than just Maquipucuna has a result -- plot_factorial_experiment()
# already skips gracefully for any site that doesn't have one yet.
for (site in SITES) {
  safe_plot(paste("factorial", site), plot_factorial_experiment(site, "reproduction_factorial_v3_h0.40"))
}

cat("\n== Default (unswept) colonization runs ==\n")
for (site in SITES) {
  safe_plot(paste("default run", site), plot_default_colonization_run(site))
}

cat("\n== Persistence validation runs (best-case / realistic / realistic_273founders) ==\n")
# Loops every site -- plot_default_colonization_run() already skips
# gracefully for any site/tag combination without a result yet (e.g.
# MiradorMindo/Yanayacu's best_case_h0.25, which crashed on 0/5 replicates --
# see project memory on the 2026-07-16 NA-propagation investigation, still
# open). realistic_273founders_h0.25 is the run_full_analysis_pipeline.R
# rerun (2026-07-18) -- see its header note on why 273 founders instead of
# realistic.rds's literature-default 30.
for (site in SITES) {
  safe_plot(paste("best_case", site), plot_default_colonization_run(site, exp_tag = "best_case_h0.40"))
  safe_plot(paste("realistic", site), plot_default_colonization_run(site, exp_tag = "realistic_h0.40"))
  safe_plot(paste("realistic_273founders", site),
            plot_default_colonization_run(site, exp_tag = "realistic_273founders_h0.40"))
}

cat("\nAll done. Plots (where inputs existed) are in", OUTPUT_DIR, "\n")
