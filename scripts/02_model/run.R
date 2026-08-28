# run.R — single entry point for the 02_model R scripts (setup,
# single-run driver, diagnostics, plots).
#
# The scripts dispatched to below remain separate files on purpose:
# get_colonization.R is the ~2000-line model engine itself (not a driver, not
# dispatched here — every driver below sources it directly), and each driver/
# diagnostic script has its own top-level argument parsing and side effects
# that would collide if concatenated into one file. This dispatcher just
# gives you one place to remember instead of fifteen filenames.
#
# For the SLURM-submitted .sh scripts (which carry their own #SBATCH resource
# headers and can't be dispatched from here), see run.sh in this directory.
#
# Usage: Rscript scripts/02_model/run.R <subcommand> [args...]
# Run from: /home/s38leste_hpc/canopymicroenv/
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
HERE <- "scripts/02_model"
all_args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: Rscript scripts/02_model/run.R <subcommand> [args...]\n\n",
    "Setup:\n",
    "  params                        Build every parameter-sweep RDS (make_params.R). No args.\n",
    "  isolation-files                Build per-(site,species) isolation params + manifest\n",
    "                                  (make_isolation_species_files.R). No args.\n",
    "  niches [HEIGHT_STEP=0.4]       Pool cross-site species niches (characterize_niches.R).\n\n",
    "Colonization runs:\n",
    "  onesite [SITE] [PARAMS_FILE] [EXP_TAG] [HEIGHT_STEP] [SPECIES_FILE]\n",
    "                                  Single-site colonization run (run_colonization.R).\n",
    "                                  Defaults: Maquipucuna, literature params, \"default\", 0.4, all species.\n\n",
    "Resolution justification:\n",
    "  resolution-diagnostics [SITE] [HEIGHT_STEPS] [HORIZ_RES]\n",
    "                                  Timing only (resolution_diagnostics.R). Comma-separated lists.\n",
    "  resolution-experiment [SITE] [HEIGHT_STEPS] [PARAMS_FILE]\n",
    "                                  Full outcome comparison (height_resolution_experiment.R).\n\n",
    "Analysis / diagnostics (read-only unless noted):\n",
    "  climate-variation               Vertical/elevational climate variation tests\n",
    "                                  (climate_variation_test.R). No args.\n",
    "  competition [BASELINE_TAG=realistic_273founders]\n",
    "                                  Isolation-vs-multi-species comparison (competition_analysis.R).\n",
    "  summarize                       Summarize all results so far (summarize_all_results.R). No args.\n",
    "  check-niche [SITE] [EXTRA_SPECIES] [HEIGHT_STEP]\n",
    "                                  Inspect per-species niche suitability (check_niche_suitability.R).\n",
    "  check-transitions [PARAMS_PATH]  Inspect stage-transition rates (check_transition_rates.R).\n",
    "  check-run [SITE] [EXP_TAG]       Inspect one completed colonization run (check_colonization_run.R).\n\n",
    "Plots:\n",
    "  plots                            Generate every project figure from whatever results exist\n",
    "                                    (plot_all.R). No args.\n",
    sep = ""
  )
}

if (length(all_args) < 1) { usage(); quit(status = 1) }
cmd <- all_args[1]
rest <- if (length(all_args) >= 2) all_args[-1] else character(0)

# Sub-scripts read commandArgs(trailingOnly = TRUE) directly, so re-expose
# just the remaining args to whichever one we dispatch to.
commandArgs <- function(...) rest  # nolint: shadows base::commandArgs on purpose

script <- switch(cmd,
  params               = "setup/make_params.R",
  "isolation-files"     = "setup/make_isolation_species_files.R",
  niches                = "setup/characterize_niches.R",
  onesite                = "run/run_colonization.R",
  "resolution-diagnostics" = "resolution/resolution_diagnostics.R",
  "resolution-experiment"  = "resolution/height_resolution_experiment.R",
  "climate-variation"    = "analysis/climate_variation_test.R",
  competition            = "analysis/competition_analysis.R",
  summarize               = "analysis/summarize_all_results.R",
  "check-niche"          = "diagnostics/check_niche_suitability.R",
  "check-transitions"    = "diagnostics/check_transition_rates.R",
  "check-run"             = "diagnostics/check_colonization_run.R",
  plots                   = "plots/plot_all.R",
  NULL
)

if (is.null(script)) { usage(); quit(status = 1) }
source(file.path(HERE, script))
