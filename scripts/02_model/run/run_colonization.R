# run_colonization.R
# Canopy colonization model — vertical niche partitioning of epiphytic Maxillariinae
# Runs one site × one parameter set. Designed to be called interactively or
# from a SLURM job array:
#   Rscript run_colonization.R <site> <params_file> <experiment_tag> <height_step> <species_file>
#
# Arguments (all optional — fall back to defaults if omitted):
#   site           Site name matching Area_or_Site in combinedv3.csv
#                  Default: "Maquipucuna"
#   params_file    Path to an RDS file containing the params list
#                  Default: uses the literature params defined below
#   experiment_tag Label appended to output file names for identification
#                  Default: "default"
#   height_step    Microenv height-tier spacing to run at (must already exist —
#                  see height_res_array.sh). Default: 0.25 (production
#                  resolution — validated against 0.1m by
#                  height_resolution_experiment.R, no significant outcome
#                  difference, ~2x cheaper). Appended to the output filename
#                  and log name whenever it isn't 0.1m, so runs at different
#                  resolutions never silently overwrite each other.
#   species_file   Path to an RDS file containing a character vector of
#                  FinalID values (see combinedv3.csv) -- REPLACES the
#                  modeled species list with exactly this set, instead of
#                  every species observed at the site (see
#                  init_colonization()'s params$species_subset,
#                  get_colonization.R). A listed species doesn't need to
#                  have been observed at this site -- e.g. to ask "how would
#                  species X, characterized from other sites, do in a
#                  landscape it's never been recorded in?" A species with no
#                  local observations here gets its niche-match ceiling from
#                  this landscape's own best-available height instead of a
#                  local realized-presence height (see niche_ceiling()) --
#                  a genuinely weaker claim, since it says nothing about
#                  whether the species could actually establish here.
#                  Separate from params_file so the same subset can be
#                  reused across any sensitivity experiment without
#                  duplicating it into every params list. Build one with e.g.
#                  saveRDS(c("SpeciesA", "SpeciesB"), "data/params/species_subset.rds")
#                  Default: NULL (model every species observed at the site).
#
# Output: data/processed/colonization_<site>_<tag>[_h<step>].rds
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(plotly)
library(ggplot2)
library(patchwork)
library(parallel)
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

# ── Command-line arguments (for cluster / batch runs) ────────────────────────
args         <- commandArgs(trailingOnly = TRUE)
site_name    <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Maquipucuna"
params_file  <- if (length(args) >= 2 && nzchar(args[2])) args[2] else NULL
exp_tag      <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "default"
height_step  <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 0.25
species_file <- if (length(args) >= 5 && nzchar(args[5])) args[5] else NULL
manifest_suffix <- if (height_step != 0.1) sprintf("_h%.2f", height_step) else ""

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOGS_DIR,
  sprintf("colonization_%s_%s%s_%s.log", site_name, exp_tag, manifest_suffix,
          format(Sys.time(), "%Y%m%d_%H%M%S")))
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("run_colonization.R started")

# ── 1. Load microenvironment ──────────────────────────────────────────────────
microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, manifest_suffix))
if (!file.exists(microenv_path)) {
  stop(sprintf("No microenv at %s (height_step=%.2f) — run height_res_array.sh %s first.",
              microenv_path, height_step, site_name))
}
microenv      <- readRDS(microenv_path)

# If this RDS was saved before .weather was added, patch it once here.
if (is.null(microenv$.weather)) {
  log_msg("Patching microenv with ERA5 weather from point model...")
  pm <- readRDS(file.path(PROCESSED_DIR, sprintf("pointmodel_%s.rds", site_name)))
  microenv$.weather <- pm[[1]]$weather
  saveRDS(microenv, microenv_path)
  rm(pm)
  log_msg("Patch saved.")
}
available_heights <- microenv_heights(microenv)
log_msg(sprintf("Heights available: %d (%.1f-%.1fm)",
  length(available_heights), min(available_heights), max(available_heights)))

# ── 2. Field observations ─────────────────────────────────────────────────────
niches <- load_observations()
niches <- niches[
  !is.na(niches$lat) & !is.na(niches$lon) &
  !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
log_msg(sprintf("Observations for %s: %d across %d species", site_name,
  sum(niches$Area_or_Site == site_name),
  length(unique(niches$FinalID[niches$Area_or_Site == site_name]))))

# ── 3. Canopy grid (fallback) ─────────────────────────────────────────────────
mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == site_name],
                    na.rm = TRUE)
log_msg(sprintf("Mean canopy height: %.1f m", mean_canopy))
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)  # used only if forestparams=NULL

# ── 4. Site object ────────────────────────────────────────────────────────────
site <- list(Site = site_name)

# ── 4b. Forest structure parameters ──────────────────────────────────────────
# Structural values from Myster (2017) primary Maquipucuna cloud forest plots:
#   mean dsh = 22.7 cm → trunk_r = 0.114 m
#   stem density = 272–324 trees/ha (≥10 cm dsh) → stems_per_ha = 298
# Crown geometry from pantropical allometry (Williams et al. 2019 review).
# Epiphyte footprint: ~4×5 cm pseudobulb cluster = 0.02 m² (field estimate).
forestparams <- list(
  stems_per_ha          = 298,    # Myster (2017) Table 5, mean of 4 primary MR plots
  mean_hgt              = 8.4,    # m — mean canopy height at this elevation
  sd_hgt                = 3.5,
  mean_crown_r          = 2.0,    # m — crown radius
  sd_crown_r            = 0.8,
  trunk_r               = 0.114,  # m — mean dsh 22.7 cm / 2 (Myster 2017)
  branch_density        = 3.0,    # m² branch surface per m² projected crown area
  epiphyte_footprint_m2 = 0.02    # m² bark area per Maxillariinae individual
)

# ── 5. Parameters ─────────────────────────────────────────────────────────────
# If a params file was passed, load it — otherwise use literature defaults.
# Sensitivity-experiment RDS files (built by make_params.R) store the swept
# field as a vector of candidate values (e.g. p_poll = c(0.05, ..., 0.70));
# every other field stays scalar. Step 7 detects that vector and sweeps over
# it via run_experiment() instead of doing a single runcolonization() call.
#   Rscript run_colonization.R Maquipucuna data/params/p_poll.rds pollination_success
if (!is.null(params_file) && file.exists(params_file)) {
  params <- readRDS(params_file)
  log_msg(sprintf("Loaded params from %s", params_file))
} else {
  log_msg("Using default literature params.")
  params <- list(
    # ── s(z, e): survival (monthly-compounded — see survival_logit() in ────
    # get_colonization.R for the p_month = p_annual^(1/12) derivation) ─────
    beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
    beta1   =  0.10,
    s_S_min =  0.0,  s_S_max =  1.0,
    s_J_min =  1.0,  s_J_max =  7.0,
    s_A_min =  7.0,  s_A_max = 20.0,
    # ── g(z'|z, e): growth / stage transitions (monthly-compounded — ───────
    # q_month = 1-(1-p_annual)^(1/12), see growth_prob() in get_colonization.R)
    psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_s_base =  0.80,
    cost_repro   =  0.50,
    # ── p_r(s) × f_s(s): fecundity ──────────────────────────────────────
    p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
    # ── d(x'|x): dispersal ──────────────────────────────────────────────
    canopy_z = mean_canopy,  lambda = 1,  Ut = 1,
    # ── spin-up ────────────────────────────────────────────────────────────
    # Founders per species, decoupled from however many field observations
    # that species has (see get_colonization.R::run_spinup()). Note: even at
    # the most generous swept p_poll/p_germ/p_s1 combination, expected seed
    # output is ~0.001/adult/year at founder size (z=s_A_min) — with only
    # ~30 founders that's still under 1 expected seed/year, so this may need
    # to go higher once real run results come back thin.
    n_founders = 30
  )
}
# canopy_z is site-specific (mean canopy height); always set it from this
# site's observations, overriding whatever a shared sensitivity-experiment
# params file may have carried.
params$canopy_z <- mean_canopy
if (is.null(params$n_founders)) params$n_founders <- 30
if (!is.null(species_file)) {
  if (!file.exists(species_file)) stop("No species_file at ", species_file)
  params$species_subset <- readRDS(species_file)
  log_msg(sprintf("species_file: restricting to %d species from %s",
                  length(params$species_subset), species_file))
}

# ── 6. Sanity check ───────────────────────────────────────────────────────────
# Shape/column checks alone pass on a structurally valid but entirely-NaN
# data frame -- that's exactly what happened against a stale height-file
# reader (get_colonization.R, lookup_climate_by_height()) before it was
# updated for the new h$tme format: every check below the `all(...)` line
# was satisfied while every climate value was NaN. Assert finiteness on
# every climate column explicitly, not just presence, so a similar reader/
# format mismatch fails loudly here instead of silently propagating into
# run_replicated().
clim_test <- lookup_climate_by_height(available_heights[1], microenv)
stopifnot(
  is.data.frame(clim_test),
  nrow(clim_test) > 0,
  nrow(clim_test) %% 24 == 0,
  all(c("month", "temp", "relhum", "windspeed", "swdown", "precip", "winddir") %in% names(clim_test))
)
clim_cols <- c("temp", "relhum", "windspeed", "swdown", "precip", "winddir")
finite_frac <- sapply(clim_cols, function(col) mean(is.finite(clim_test[[col]])))
if (any(finite_frac == 0)) {
  stop(sprintf(
    "lookup_climate_by_height() returned an all-non-finite column (%s) for height %.2f -- climate data is not being read correctly (stale reader vs. height-file format?). Column finite fractions: %s",
    paste(names(finite_frac)[finite_frac == 0], collapse = ", "),
    available_heights[1],
    paste(sprintf("%s=%.2f", names(finite_frac), finite_frac), collapse = ", ")
  ))
}
log_msg(sprintf("Climate check passed: %.1f°C mean temp, %.0f mm/yr precip",
  mean(clim_test$temp, na.rm = TRUE),
  mean(clim_test$precip, na.rm = TRUE) * 8760))

# ── 7. Run colonization model ─────────────────────────────────────────────────
# A one-at-a-time sensitivity params file has exactly one vector-valued
# field — sweep it with run_experiment(). A factorial params file (e.g.
# make_params.R's reproduction factorial) has several — cross them with
# run_factorial_experiment(). A plain params file (or literature defaults)
# has none, so run the model once via runcolonization().
# RUN_TIMESTEPS: shared run length (years) for every branch below -- raised
# from 30 to 50 (2026-07-24) while microenv generation was still running, so
# it applies to every colonization run launched after this point (including
# persistence validation) rather than only new experiment types. Existing
# colonization_*.rds result files (best_case/realistic, reproduction_
# factorial_v3, OAT sweeps) were all generated at timesteps=30 AND under
# pre-size-tracking/pre-niche-rewrite code, so they need regenerating anyway
# to be comparable to anything run after this change -- see methods.tex,
# One-at-a-time experiments / Full factorial experiment, for the description
# that needs to be updated to 50 once fresh results land under this setting.
RUN_TIMESTEPS <- 50
swept_params <- names(params)[vapply(params, length, integer(1)) > 1]

if (length(swept_params) > 1) {
  N_CORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
  if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)
  N_REPS <- 1  # replicates per combo — factorials get large fast (see paper's 729-combo x1 design)
  param_values <- setNames(lapply(swept_params, function(nm) params[[nm]]), swept_params)
  n_combos <- prod(vapply(param_values, length, integer(1)))
  # Checkpoint by the last-listed swept parameter (expand.grid's slowest-
  # varying dimension, so the natural "outer loop" grouping) so a walltime
  # kill loses at most one block's worth of jobs instead of the whole sweep.
  checkpoint_var  <- swept_params[length(swept_params)]
  checkpoint_path <- file.path(PROCESSED_DIR,
    sprintf("colonization_%s_%s%s_checkpoint.rds", site_name, exp_tag, manifest_suffix))
  log_msg(sprintf("Factorial sweep %s: %d combos x %d reps on %d cores (checkpointing by %s to %s)...",
                  paste(swept_params, collapse = " x "), n_combos, N_REPS, N_CORES,
                  checkpoint_var, checkpoint_path))
  result <- run_factorial_experiment(param_values, base = params,
                                     n_reps = N_REPS, timesteps = RUN_TIMESTEPS, spinup = 5,
                                     checkpoint_var = checkpoint_var,
                                     checkpoint_path = checkpoint_path)
} else if (length(swept_params) == 1) {
  swept_param <- swept_params
  N_CORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
  if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)
  N_REPS <- 3  # replicates per swept value
  log_msg(sprintf("Sweeping %s across %d values [%s] with %d reps on %d cores...",
                  swept_param, length(params[[swept_param]]),
                  paste(params[[swept_param]], collapse = ", "), N_REPS, N_CORES))
  result <- run_experiment(swept_param, params[[swept_param]], base = params,
                            n_reps = N_REPS, timesteps = RUN_TIMESTEPS, spinup = 5)
} else {
  # No swept parameter. params$n_reps (default 1) controls how many
  # independent replicates to run — more than 1 is useful to tell "genuinely
  # blocked" apart from "one unlucky stochastic draw" (see
  # run_replicated()). n_reps=1 behaves like a single runcolonization() call,
  # just wrapped in the same list(runs=, summary=) structure for consistency.
  N_CORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
  if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)
  N_REPS <- if (!is.null(params$n_reps)) params$n_reps else 1
  log_msg(sprintf("Starting colonization run [%s | %s | height_step=%.2f]: %d replicate(s) on %d cores...",
                  site_name, exp_tag, height_step, N_REPS, N_CORES))
  result <- run_replicated(params, n_reps = N_REPS, timesteps = RUN_TIMESTEPS, spinup = 5)
}

out_path <- file.path(PROCESSED_DIR,
  sprintf("colonization_%s_%s%s.rds", site_name, exp_tag, manifest_suffix))
saveRDS(result, out_path)
log_msg(sprintf("Done. Results saved to %s", out_path))

# ── 8. Plots (interactive only) ───────────────────────────────────────────────
if (interactive()) {
  if (is.data.frame(result) && length(swept_params) == 1) {
    print(plot_experiment(result, swept_params))
  } else if (is.data.frame(result)) {
    message("Factorial result — no dedicated plot yet; inspect the data frame directly ",
            "(columns: ", paste(swept_params, collapse = ", "), ", t, totalS/J/A, extinct).")
  } else if (is.list(result) && !is.null(result$runs)) {
    # No swept parameter: result$runs is one runcolonization() output per
    # replicate. Plot the first successful one.
    first_ok <- Filter(Negate(is.null), result$runs)
    if (length(first_ok) > 0) {
      plot_abundance(first_ok[[1]])
      plot_3d_abundance(first_ok[[1]])
    } else {
      message("No replicate succeeded — nothing to plot.")
    }
  }
}
