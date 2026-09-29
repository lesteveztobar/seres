# height_resolution_experiment.R
# Outcome diagnostic (as opposed to resolution_diagnostics.R's timing-only
# diagnostic): runs the actual colonization model at each microenv height
# resolution (0.1 / 0.25 / 0.5 / 1.0m) with real timesteps/spinup and several
# replicates each, to check whether height resolution changes RESULTS
# (population trajectories, extinction rate) and not just runtime.
#
# Defaults to realistic.rds (2026-07-24, was best_case.rds) -- literature-
# calibrated values throughout, rather than the most generous corner tested,
# so the resolution comparison reflects the parameter regime the model is
# actually meant to represent. CAVEAT: as of the last check under the OLD
# code/timesteps, realistic.rds collapsed toward extinction rather than
# persisting (see report/methods.tex, Persistence validation) -- that result
# predates the niche density-ratio rewrite, the monthly precip fix, and the
# S/J size-tracking overhaul, so it needs re-verifying under current code
# before trusting this comparison. If realistic.rds still goes extinct
# across most/all height steps, every trajectory converges to zero
# regardless of resolution and the comparison becomes uninformative --
# check data/processed/height_resolution_<site>_realistic.rds's extinction
# rate before reading anything into the height-step effect size. Any params
# file can still be passed explicitly -- e.g. best_case.rds, to fall back to
# the original non-zero-trajectories guarantee if that turns out necessary.
#
# Horizontal resolution is held fixed at 10 (matches production run_one()) --
# only height resolution varies here.
#
# Requires microenv_<site>.rds (0.1m) and microenv_<site>_h<step>.rds for
# coarser steps to already exist (height_res_array.sh).
#
# Usage: Rscript scripts/experiments/A11_height_resolution/height_resolution_experiment.R <site> [height_steps] [params_file]
#   e.g.: Rscript scripts/experiments/A11_height_resolution/height_resolution_experiment.R Maquipucuna "0.1,0.25,0.5,1.0" data/params/best_case.rds
# Output: data/processed/height_resolution_<site>_<params_tag>.rds (tidy:
#   height_step, n_heights, rep, t, totalS/J/A, total, extinct) -- params_tag
#   is the params_file's basename (e.g. "best_case", "realistic"), so runs
#   with different params files don't silently overwrite each other.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(parallel)
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

args         <- commandArgs(trailingOnly = TRUE)
site_name    <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Maquipucuna"
height_steps <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(strsplit(args[2], ",")[[1]]) else c(0.1, 0.25, 0.5, 1.0)
params_file  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else file.path(PARAMS_DIR, "realistic.rds")
params_tag   <- tools::file_path_sans_ext(basename(params_file))

TIMESTEPS  <- 50  # matches production run_colonization.R (RUN_TIMESTEPS, raised from 30 2026-07-24)
SPINUP     <- 5
RESOLUTION <- 10  # fixed -- isolates height-resolution effect only

N_CORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)

source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg()

# ── Site setup (mirrors run_colonization.R) ─────────────────────────
niches <- load_observations()
niches <- niches[
  !is.na(niches$lat) & !is.na(niches$lon) &
  !is.na(niches$Height_m) & !is.na(niches$FinalID), ]

# v7 (Phase 1.6): shared ceiling formula, site_canopy_ceiling() (shared_helpers.R).
# No microenv is loaded yet at this point (loaded per-resolution-step,
# below) to supply a fallback -- stop() loudly rather than silently
# propagate NaN if this site ever lacks measured CanopyHeight_m (in
# practice this script only targets Maquipucuna, which has it).
# 2026-09-30: site_canopy_ceiling() always evaluates its fallback since the
# 2026-09-10 change, so the old stop() fallback aborted every run. Use the
# production (0.4 m) manifest's top height tier, exactly as run_colonization.R does.
.ref_heights <- microenv_heights(readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site_name))))
mean_canopy <- site_canopy_ceiling(site_name, niches, max(.ref_heights))
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)
site <- list(Site = site_name)

forestparams <- site_forestparams(site_name, mean_canopy)  # same per-site rule as run_colonization.R (was default_forestparams(): Maquipucuna geometry at every site)

params <- readRDS(params_file)
params$canopy_z <- mean_canopy
N_REPS <- if (!is.null(params$n_reps)) params$n_reps else 5
log_msg(sprintf("Loaded params from %s (n_reps=%d)", params_file, N_REPS))

# ── Run one height step's worth of replicates ───────────────────────────────
run_height_step <- function(step) {
  suffix <- if (step != 0.1) sprintf("_h%.2f", step) else ""
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, suffix))
  if (!file.exists(microenv_path)) {
    log_msg(sprintf("Skipping height step %.2fm -- no %s (run height_res_array.sh first)",
                    step, microenv_path))
    return(NULL)
  }

  # run_one()/run_replicated() (get_colonization.R) read `microenv` as a
  # global via lexical scoping -- update it here so they pick up this step's.
  microenv <<- readRDS(microenv_path)
  n_heights <- length(microenv_heights(microenv))
  log_msg(sprintf("Height step %.2fm: %d tiers. Building climate cache...", step, n_heights))
  clim_cache <- build_clim_cache(microenv)

  log_msg(sprintf("Height step %.2fm: running %d replicate(s) x %d timesteps (spinup %d)...",
                  step, N_REPS, TIMESTEPS, SPINUP))
  out <- run_replicated(params, n_reps = N_REPS, timesteps = TIMESTEPS,
                        spinup = SPINUP, clim_cache = clim_cache)

  df <- out$summary
  if (is.null(df) || nrow(df) == 0) {
    log_msg(sprintf("Height step %.2fm: no replicate succeeded.", step))
    return(NULL)
  }
  df$height_step <- step
  df$n_heights   <- n_heights
  df
}

all_rows   <- lapply(height_steps, run_height_step)
summary_df <- do.call(rbind, Filter(Negate(is.null), all_rows))

if (is.null(summary_df) || nrow(summary_df) == 0) {
  stop("No height-step variants produced usable results for ", site_name)
}

out_path <- file.path(PROCESSED_DIR, sprintf("height_resolution_%s_%s.rds", site_name, params_tag))
saveRDS(summary_df, out_path)

# ── Descriptive comparison ──────────────────────────────────────────────────
final_df <- summary_df[summary_df$t == TIMESTEPS, ]

cat("\n== Final total abundance (t=", TIMESTEPS, ") by height step ==\n", sep = "")
print(aggregate(total ~ height_step, data = final_df,
                FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))))

# `extinct` is constant within a (height_step, rep) -- one row per replicate.
extinct_per_rep <- aggregate(extinct ~ height_step + rep, data = summary_df,
                             FUN = function(x) x[1])
cat("\n== Extinction rate by height step ==\n")
print(aggregate(extinct ~ height_step, data = extinct_per_rep, FUN = mean))

# Supplementary omnibus + pairwise checks -- indicative only given small
# n_reps. The omnibus Kruskal-Wallis only says the groups differ somewhere;
# the Dunn post-hoc (Holm-corrected) pairwise comparisons answer the actual
# production question -- which coarser resolution(s), if any, are still
# statistically indistinguishable from the finest (0.1m, the reference/most-
# accurate resolution)? Same method as plot_temperature_profile()'s height
# comparison (plot_functions.R) for consistency across the project.
if (length(unique(final_df$height_step)) > 1 && nrow(final_df) >= 4) {
  kt <- tryCatch(kruskal.test(total ~ factor(height_step), data = final_df),
                error = function(e) NULL)
  if (!is.null(kt)) {
    cat("\n== Kruskal-Wallis test: final total abundance ~ height_step ==\n")
    print(kt)
  }

  if (!requireNamespace("dunn.test", quietly = TRUE)) install.packages("dunn.test")
  # label=TRUE (not FALSE): comparisons need the actual height_step values
  # ("0.1 - 0.25") for the "vs. finest" filter below to match against --
  # label=FALSE gives bare group indices ("1 - 2") instead, which silently
  # produced an empty "vs. finest" sub-table on the first real run
  # (2026-07-18, Maquipucuna/realistic_273founders -- the full pairwise
  # table was still correct, just unfiltered).
  dunn_res <- tryCatch(
    dunn.test::dunn.test(final_df$total, factor(final_df$height_step),
                         method = "holm", kw = FALSE, label = TRUE),
    error = function(e) NULL
  )
  if (!is.null(dunn_res)) {
    dunn_df <- data.frame(
      comparison = dunn_res$comparisons,
      p_adj      = dunn_res$P.adjusted,
      sig        = .sig_stars(dunn_res$P.adjusted)  # shared_helpers.R
    )
    cat("\n== Pairwise Dunn tests (Holm-adjusted): every height_step pair ==\n")
    print(dunn_df[order(dunn_df$p_adj), ], row.names = FALSE)

    # Sub-table: only the pairs that include the finest resolution actually
    # tested -- "ns" here means that coarser step is a safe stand-in for the
    # finest one; "*"/"**"/"***" means it isn't.
    finest_label     <- as.character(sort(unique(final_df$height_step))[1])
    comp_pairs       <- strsplit(dunn_df$comparison, " - ")
    involves_finest  <- vapply(comp_pairs, function(p) finest_label %in% trimws(p), logical(1))
    vs_finest        <- dunn_df[involves_finest, ]
    cat(sprintf("\n-- Coarser resolutions vs. the finest tested (%sm) --\n", finest_label))
    print(vs_finest[order(vs_finest$p_adj), ], row.names = FALSE)
    cat("(ns = statistically indistinguishable from the finest resolution -- safe to use)\n")
  }
}

cat(sprintf("\nSaved: %s\n", out_path))
