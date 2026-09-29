# sensitivity_run.R -- sensitivity step 3: evaluate one array task's slice
# of the LHS design and write one tidy CSV row per run.
#
# env: CANOPY_ARRAY_IDX (1-based), CANOPY_ARRAY_N (total tasks),
#      SLURM_CPUS_PER_TASK (mclapply width), CANOPY_CLIM_MODE=voxel,
#      CANOPY_SENS_DESIGN (design rds, default sensitivity_design.rds),
#      CANOPY_SENS_OUTDIR (subdir of OUTPUT_DIR, default sensitivity_runs),
#      CANOPY_SENS_TIMESTEPS (default 50)
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")
source("scripts/02_model/lib_logging.R")

# 2026-09-12: CANOPY_ARRAY_IDX/CANOPY_ARRAY_N used to default to "1" when
# unset -- silently turning a forgotten env var into "task 1 of 1, run the
# ENTIRE design serially" instead of an error. Every launcher sets both;
# require them explicitly rather than guessing.
.require_env <- function(name) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v)) stop(sprintf("%s is not set -- required (no default; a missing array index/count must not silently mean \"run everything as task 1 of 1\").", name))
  v
}
IDX <- as.integer(.require_env("CANOPY_ARRAY_IDX"))
NTASK <- as.integer(.require_env("CANOPY_ARRAY_N"))
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
TIMESTEPS <- as.integer(Sys.getenv("CANOPY_SENS_TIMESTEPS", unset = "50"))
DESIGN_PATH <- Sys.getenv("CANOPY_SENS_DESIGN", unset = "data/params/sensitivity_design.rds")
OUT_SUBDIR <- Sys.getenv("CANOPY_SENS_OUTDIR", unset = "sensitivity_runs")
Sys.setenv(CANOPY_CLIM_MODE = "voxel")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("sens_run_%s_%03d_%s.log", OUT_SUBDIR, IDX, format(Sys.time(), "%Y%m%d_%H%M%S"))))

design <- readRDS(DESIGN_PATH)
SKIP <- strsplit(Sys.getenv("CANOPY_SKIP_SITES", unset = ""), "[ ,]+")[[1]]
if (length(SKIP) && any(nzchar(SKIP))) {
  message(sprintf("CANOPY_SKIP_SITES=%s: dropping %d design rows", paste(SKIP, collapse = ","), sum(design$site %in% SKIP)))
  design <- design[!design$site %in% SKIP, ]
}
design <- design[order(design$site, design$run_id), ]
slice_idx <- which((seq_len(nrow(design)) - 1L) %% NTASK + 1L == IDX)
mine <- design[slice_idx, ]
log_msg(sprintf("task %d/%d: %d of %d runs (design=%s, timesteps=%d)", IDX, NTASK, nrow(mine), nrow(design), DESIGN_PATH, TIMESTEPS))

OUT_DIR <- file.path(OUTPUT_DIR, OUT_SUBDIR)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(OUT_DIR, sprintf("part_%03d.csv", IDX))

# 2026-09-24: PARAM_BASE used to be a hard-coded literal list with
# lambda = 1, Ut = 1 (run_colonization.R's fallback defaults) -- but every
# project params file (realistic.rds, best_case.rds, the factorials) carries
# the calibrated dispersal values lambda = 3.23, Ut = 0.23. The whole
# 2026-09-11 sensitivity analysis therefore ran with a different dispersal
# kernel from the rest of the project, silently. Derive the base from the
# project's own realistic.rds instead, so it cannot drift again; only the
# five swept parameters are overwritten per run.
PARAM_BASE <- readRDS("data/params/realistic.rds")
for (nm in c("p_poll", "p_germ", "p_s1", "S", "n_founders", "n_reps")) PARAM_BASE[[nm]] <- NULL
stopifnot(!is.null(PARAM_BASE$lambda), !is.null(PARAM_BASE$Ut), abs(PARAM_BASE$lambda - 3.23) < 1e-9)

# per-site fixtures (microenv, niches, ceiling, forestparams, clim_cache) built once
.fix <- new.env()
site_fixture <- function(site) {
  if (!is.null(.fix[[site]])) return(.fix[[site]])
  microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site)))
  niches <- load_observations()
  niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                   !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
  ceil <- site_canopy_ceiling(site, niches, max(microenv_heights(microenv)))
  # Reuse ONLY the pooled hourly clim_cache across this site's runs (cheap,
  # identical every run -- same as run_replicated()). The per-voxel caches
  # are rebuilt inside each run: they are footprint-restricted (~50 pixels)
  # so the rebuild is ~2 min but every .clim_voxel_slice() is then a
  # trivial match(); a reused full-raster voxel cache made those match()es
  # scan ~140k pixels and ran ~5x slower overall.
  f <- list(microenv = microenv, niches = niches, ceiling = ceil,
            canopy_grid = matrix(ceil, 50, 50),
            forestparams = site_forestparams(site, ceil),
            heights = microenv_heights(microenv),
            clim_cache = build_clim_cache(microenv))
  .fix[[site]] <- f
  f
}

realised_mean_height <- function(run, heights, at_t = NULL) {
  ft <- if (is.null(at_t)) dim(run$abundanceA)[4] else at_t
  h <- numeric(0)
  for (arr in list(run$abundanceA, run$abundanceJ, run$abundanceS)) {
    sl <- arr[, , , ft, , drop = FALSE]
    occ <- which(sl > 0, arr.ind = TRUE)
    if (nrow(occ) == 0) next
    h <- c(h, rep(heights[occ[, 3]], times = sl[occ]))
  }
  if (length(h) == 0) return(c(mean = NA_real_, n = 0))
  c(mean = mean(h), n = length(h))
}

one_run <- function(row) {
  # 2026-09-12: `stage` is now written by the harness itself, once, here --
  # not patched on after the fact by whichever combine script happens to
  # run later (that pattern is exactly what let 3,000 stage-2 rows go out
  # with no stage tag at all on 2026-09-11, silently pooled into "stage 1"
  # by a downstream fallback). Falls back to 1L only if the design itself
  # never had a `stage` column (stage 1's own design predates the column).
  stage_val <- if (!is.null(row$stage) && !is.na(row$stage)) as.integer(row$stage) else 1L
  fx <- site_fixture(row$site)
  p <- PARAM_BASE
  p$p_poll <- row$p_poll; p$p_germ <- row$p_germ; p$p_s1 <- row$p_s1
  p$S <- row$S; p$n_founders <- as.integer(row$n_founders)
  p$canopy_z <- fx$ceiling
  # optional per-run carrying-capacity multiplier (K sweep): only the PRODUCT
  # occupiable_bark_fraction x maxillariinae_community_share enters K, so the
  # sweep is on that product's multiplier (see build_forest()).
  Sys.setenv(CANOPY_K_MULT = if (!is.null(row$k_mult) && !is.na(row$k_mult)) format(row$k_mult, digits = 10) else "1")
  t0 <- Sys.time()
  r <- tryCatch(suppressMessages(runcolonization(
    site = list(Site = row$site), niches = fx$niches, canopy_grid = fx$canopy_grid,
    microenv = fx$microenv, timesteps = TIMESTEPS, resolution = 10, carCap = 1, maxDisp = 5,
    stochastic = FALSE, Visualize = FALSE, spinup = 5,
    parameters = p, forestparams = fx$forestparams, seed = row$seed,
    clim_cache = fx$clim_cache
  )), error = function(e) structure(list(err = conditionMessage(e)), class = "runfail"))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "runfail")) {
    return(data.frame(row[, c("run_id", "point_id", "site", "rep", "seed",
                              "p_poll", "p_germ", "p_s1", "S", "n_founders")],
                      final_total = NA, final_S = NA, final_J = NA, final_A = NA,
                      recruited = NA, extinct = NA, regime = "error", k_bind_t = NA,
                      exceeded_plausible_density = NA, total_K = NA,
                      realised_mean_height = NA, n_individuals = NA,
                      t50_total = NA, t50_realised_mean_height = NA,
                      stage = stage_val, k_mult = if (!is.null(row$k_mult)) row$k_mult else 1,
                      run_seconds = secs, error = r$err, stringsAsFactors = FALSE))
  }
  if (is.null(r$k_bind_t)) {
    log_msg(sprintf("WARN [%s]: runcolonization() returned no k_bind_t -- engine/harness version mismatch? Recording NA, not silently assuming a value.", row$run_id))
  }
  n_t <- length(r$totalabundanceS)
  fS <- r$totalabundanceS[n_t]; fJ <- r$totalabundanceJ[n_t]; fA <- r$totalabundanceA[n_t]
  rh <- realised_mean_height(r, fx$heights)
  # transient check (2026-09-11): when TIMESTEPS > 50, also snapshot t=50 so
  # the SAME 300-run subset yields a controlled t=50-vs-t=150 PAWN
  # comparison -- not two different samples.
  t50 <- if (TIMESTEPS > 50 && n_t >= 50) {
    list(total = r$totalabundanceS[50] + r$totalabundanceJ[50] + r$totalabundanceA[50],
         h = realised_mean_height(r, fx$heights, at_t = 50)["mean"])
  } else list(total = NA_real_, h = NA_real_)
  data.frame(row[, c("run_id", "point_id", "site", "rep", "seed",
                     "p_poll", "p_germ", "p_s1", "S", "n_founders")],
             final_total = fS + fJ + fA, final_S = fS, final_J = fJ, final_A = fA,
             recruited = (fS + fJ) > 0,
             extinct = all(r$totalabundanceA[(n_t %/% 2):n_t] == 0),
             regime = r$regime,
             k_bind_t = if (is.null(r$k_bind_t)) NA_integer_ else r$k_bind_t,
             exceeded_plausible_density = isTRUE(r$exceeded_plausible_density),
             total_K = r$total_K,
             realised_mean_height = unname(rh["mean"]),
             n_individuals = unname(rh["n"]),
             t50_total = unname(t50$total), t50_realised_mean_height = unname(t50$h),
             stage = stage_val, k_mult = if (!is.null(row$k_mult)) row$k_mult else 1,
             run_seconds = secs, error = NA_character_, stringsAsFactors = FALSE)
}

rows <- split(mine, seq_len(nrow(mine)))
res <- parallel::mclapply(rows, function(rr) {
  out <- tryCatch(one_run(rr), error = function(e) NULL)
  log_msg(sprintf("  %s done", rr$run_id))
  out
}, mc.cores = NCORES)
res <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(res, out_path, row.names = FALSE)
log_msg(sprintf("task %d: wrote %d rows to %s", IDX, nrow(res), out_path))
