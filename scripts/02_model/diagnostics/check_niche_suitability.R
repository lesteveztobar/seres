# check_niche_suitability.R
# Diagnostic: for each species observed at a site (plus any extra species
# named on the command line — see below), reports the distribution of
# per-axis (temp/relhum/swdown) 0-100 suitability scores across every height
# tier, plus the combined score both BEFORE and AFTER the two-stage
# normalization (see niche_raw_score()/niche_ceiling()/niche_overall_score()
# in get_colonization.R):
#   BEFORE = geometric mean of the three axis scores — "multiply, then
#     normalize back to 100" (prof feedback 2026-07-15), scale-consistent
#     unlike dividing by a fixed constant.
#   AFTER  = BEFORE rescaled against a per-site ceiling (prof feedback
#     2026-07-15 follow-up) — corrects for the three axes' individual optima
#     almost never coinciding at one real height tier, which would otherwise
#     cap every species well below 100 everywhere at this site. Two modes:
#       - species observed at this site: ceiling = best score at its own
#         observed presence heights here (the realized-niche framing).
#       - species with NO observations at this site (pass extra species as
#         a 2nd CLI arg to ask "how would species X, characterized
#         elsewhere, do in a landscape it's never been recorded in?"):
#         ceiling = best score anywhere in this landscape's own vertical
#         profile instead. A genuinely weaker claim than the observed case —
#         see the "(ceiling mode: ...)" line each species prints.
# Checking both together shows exactly how much of a lift the ceiling
# rescale is giving — a large gap between BEFORE's max and 100 means this
# site's climate profile doesn't offer a height where all three axes are
# simultaneously ideal for that species, so AFTER is doing real work.
#
# Read-only, no cluster job needed — only needs paths.R + get_colonization.R
# (no sf/GDAL/PROJ/GEOS/UDUNITS), so no module-loading wrapper is required.
# For a plotted version (histograms), see plot_niche_suitability() in
# plot_functions.R / run_plots.sh.
#
# Usage: Rscript scripts/02_model/diagnostics/check_niche_suitability.R [site] [extra_species] [height_step]
#   extra_species: comma-separated FinalID values to check even if never
#   observed at this site (e.g. to preview a species_subset transplant —
#   see run_colonization.R's species_file arg).
#   height_step: defaults to 0.25 (production resolution).
#   e.g.: Rscript scripts/02_model/diagnostics/check_niche_suitability.R Maquipucuna SomeSpeciesFromMashpi
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

args          <- commandArgs(trailingOnly = TRUE)
site_name     <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Maquipucuna"
extra_species <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else character(0)
height_step   <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 0.25

niche_cache_path <- file.path(PROCESSED_DIR, "species_niches.rds")
if (!file.exists(niche_cache_path))
  stop("No niche cache at ", niche_cache_path, " -- run characterize_niches.R first")
niche_cache <- readRDS(niche_cache_path)

# height_step defaults to 0.25 (production resolution, matches
# run_colonization.R's manifest_suffix convention) -- NOT the
# unsuffixed 0.1m file, which most sites don't even have.
manifest_suffix <- if (height_step != 0.1) sprintf("_h%.2f", height_step) else ""
microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, manifest_suffix))
if (!file.exists(microenv_path)) stop("No microenv for ", site_name, " at ", microenv_path)
microenv <- readRDS(microenv_path)
heights  <- microenv_heights(microenv)

message("Building climate cache (reads all ", length(heights), " height files once)...")
cc <- build_clim_cache(microenv)
message("Done.")

# Per-height climate scalars for the three niche variables (list-of-lists,
# one per height tier, in the shape niche_axis_scores()/niche_raw_score()
# expect) and the equivalent matrix form (for the ceiling lookup below).
# landscape_clim_vals: complete-case rows only -- this landscape's own
# vertical profile, used as the ceiling fallback for species with no local
# observations at this site (see niche_ceiling()).
height_scalars      <- height_clim_scalars(cc$clim_by_height)
landscape_clim_vals <- height_scalars[stats::complete.cases(height_scalars), , drop = FALSE]
clim_by_height <- lapply(seq_len(nrow(height_scalars)), function(i) as.list(height_scalars[i, ]))
clim_by_height <- Filter(function(cl) !anyNA(unlist(cl)), clim_by_height)

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
site_obs     <- niches[niches$Area_or_Site == site_name, ]
site_species <- sort(unique(c(site_obs$FinalID, extra_species)))

.qsummary <- function(x) {
  q <- quantile(x, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  sprintf("min=%.1f  Q1=%.1f  median=%.1f  Q3=%.1f  max=%.1f", q[1], q[2], q[3], q[4], q[5])
}

cat(sprintf("\n=== %s: niche suitability across %d height tiers ===\n", site_name, length(clim_by_height)))

for (sp in site_species) {
  niche_sp <- niche_cache[[sp]]
  if (is.null(niche_sp)) {
    cat(sprintf("\n%s: no cached niche (no usable observations anywhere)\n", sp))
    next
  }

  # Per-site ceiling: this species' best niche_raw_score() at its own
  # observed presence heights on this site, or (if it has none here) the
  # best score anywhere in this landscape's own vertical profile instead —
  # same logic init_colonization() applies at runtime, see niche_ceiling().
  obs_sp <- site_obs[site_obs$FinalID == sp & !is.na(site_obs$Height_m), ]
  obs_clim_vals <- if (nrow(obs_sp) > 0) {
    do.call(rbind, lapply(obs_sp$Height_m, function(h)
      height_scalars[which.min(abs(heights - h)), , drop = TRUE]))
  } else NULL
  ceiling_mode <- if (is.null(obs_clim_vals)) "landscape-best (no local observations)" else "local observed presence"
  ceiling <- niche_ceiling(niche_sp, obs_clim_vals, landscape_clim_vals)

  axis_scores <- sapply(clim_by_height, function(clim) niche_axis_scores(clim, niche_sp))
  # axis_scores: matrix, rows = variables, cols = height tiers
  before <- apply(axis_scores, 2, function(s) .geomean(s))       # geometric mean, no ceiling
  after  <- pmin(100, 100 * before / ceiling)                    # rescaled to this site's ceiling

  n_obs <- if (is.null(obs_clim_vals)) 0L else nrow(obs_clim_vals)
  cat(sprintf("\n%s  (ceiling = %.1f, mode: %s, %d observed height(s) at this site)\n",
              sp, ceiling, ceiling_mode, n_obs))
  for (v in rownames(axis_scores)) {
    cat(sprintf("  %-8s %s\n", v, .qsummary(axis_scores[v, ])))
  }
  cat(sprintf("  %-8s %s  (geometric mean, before ceiling rescale)\n", "before", .qsummary(before)))
  cat(sprintf("  %-8s %s  (after ceiling rescale, 0-100 scale)\n", "after", .qsummary(after)))
  cat(sprintf("  %d/%d height tiers with nonzero overall suitability\n",
              sum(after > 0), length(after)))
}
