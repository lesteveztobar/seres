# sensitivity_run_v10.R -- 2026-10-02 campaign harness (two calibration arms).
# Derived from sensitivity_run.R; differences:
#   * a task is selected by the design's own `task` column (interleaved by
#     arm x site in build_v10_design.R), not by round-robin over a site sort;
#   * per-row survival intercepts beta0S/J/A (the calibration arm);
#   * per-run elapsed-time cap (CANOPY_RUN_TIMEOUT_S, default 10800 s) --
#     a capped run is recorded with regime = "timeout", not lost;
#   * INCREMENTAL output: every finished run writes its own row file and
#     per-timestep summary file under <outdir>/task_NNN/ (atomic rename),
#     so a task killed at walltime keeps every run it finished; the parent
#     then consolidates them into part_NNN.csv and summary_NNN.csv;
#   * per-timestep summary (rep, t, totalS/J/A, total, extinct, recruited,
#     regime) plus recruits_t / n_adults_t / flowering_adults_t, recorded
#     by wrapping run_pass1_disperse()/run_pass2_establish() here -- the
#     engine itself is not modified and spatial arrays are never written.
#
# env: CANOPY_ARRAY_IDX (task), CANOPY_SENS_DESIGN, CANOPY_SENS_OUTDIR,
#      SLURM_CPUS_PER_TASK, CANOPY_SENS_TIMESTEPS (default 50),
#      CANOPY_RUN_TIMEOUT_S (default 10800)
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/plots/plot_functions.R")
source("scripts/02_model/lib_logging.R")

.require_env <- function(name) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v)) stop(sprintf("%s is not set -- required.", name))
  v
}
IDX <- as.integer(.require_env("CANOPY_ARRAY_IDX"))
DESIGN_PATH <- .require_env("CANOPY_SENS_DESIGN")
OUT_SUBDIR <- .require_env("CANOPY_SENS_OUTDIR")
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
TIMESTEPS <- as.integer(Sys.getenv("CANOPY_SENS_TIMESTEPS", unset = "50"))
TIMEOUT_S <- as.numeric(Sys.getenv("CANOPY_RUN_TIMEOUT_S", unset = "10800"))
Sys.setenv(CANOPY_CLIM_MODE = "voxel")
log_msg <- make_log_msg(log_file = file.path(LOGS_DIR, sprintf("sens_run_%s_%03d_%s.log",
  gsub("/", "_", OUT_SUBDIR), IDX, format(Sys.time(), "%Y%m%d_%H%M%S"))))

design <- readRDS(DESIGN_PATH)
mine <- design[design$task == IDX, ]
if (!nrow(mine)) stop(sprintf("task %d: no rows in %s", IDX, DESIGN_PATH))
stopifnot(length(unique(mine$site)) == 1L)
log_msg(sprintf("task %d: %d runs, site=%s arm=%s block=%s (design=%s, timesteps=%d, timeout=%gs, cores=%d)",
                IDX, nrow(mine), mine$site[1], paste(unique(mine$arm), collapse = ","),
                mine$block[1], DESIGN_PATH, TIMESTEPS, TIMEOUT_S, NCORES))

OUT_DIR <- file.path(OUTPUT_DIR, OUT_SUBDIR)
TASK_DIR <- file.path(OUT_DIR, sprintf("task_%03d", IDX))
dir.create(TASK_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(OUT_DIR, sprintf("part_%03d.csv", IDX))
sum_path <- file.path(OUT_DIR, sprintf("summary_%03d.csv", IDX))

PARAM_BASE <- readRDS("data/params/realistic.rds")
for (nm in c("p_poll", "p_germ", "p_s1", "S", "n_founders", "n_reps", "beta0S", "beta0J", "beta0A")) PARAM_BASE[[nm]] <- NULL
stopifnot(!is.null(PARAM_BASE$lambda), abs(PARAM_BASE$lambda - 3.23) < 1e-9, PARAM_BASE$delta_s_base == 0.80)

site <- mine$site[1]
microenv <- readRDS(file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site)))
niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) & !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
ceil <- site_canopy_ceiling(site, niches, max(microenv_heights(microenv)))
fx <- list(microenv = microenv, niches = niches, ceiling = ceil, canopy_grid = matrix(ceil, 50, 50),
           forestparams = site_forestparams(site, ceil), heights = microenv_heights(microenv),
           clim_cache = build_clim_cache(microenv))
rm(microenv)

# ---- recorders: wrap the pass functions, values passed through unchanged ----
.rec <- new.env()
.orig_pass1 <- run_pass1_disperse
.orig_pass2 <- run_pass2_establish
run_pass1_disperse <- function(state, abundanceA, size_A, t, stochastic = FALSE) {
  if (dim(abundanceA)[4] == TIMESTEPS) {
    a <- abundanceA[, , , t, ]; z <- size_A[, , , t, ]
    occ <- which(a > 0)
    .rec$adults[t] <- sum(a[occ])
    .rec$flowering[t] <- sum(a[occ] * flowering_prob(z[occ]))
  }
  .orig_pass1(state, abundanceA, size_A, t, stochastic = stochastic)
}
run_pass2_establish <- function(state, abundanceS, abundanceJ, abundanceA, size_S, Disp, t, stochastic = FALSE) {
  out <- .orig_pass2(state, abundanceS, abundanceJ, abundanceA, size_S, Disp, t, stochastic = stochastic)
  if (dim(abundanceS)[4] == TIMESTEPS) {
    .rec$recruits[t + 1] <- sum(out$S[, , , t + 1, ]) - sum(abundanceS[, , , t + 1, ])
  }
  out
}

stage_mean_z <- function(ab, sz) {
  occ <- which(ab > 0)
  if (!length(occ)) return(NA_real_)
  sum(ab[occ] * sz[occ]) / sum(ab[occ])
}

ID_COLS <- c("run_id", "block", "arm", "task", "point_id", "site", "rep", "seed",
             "p_poll", "p_germ", "p_s1", "S", "n_founders", "beta0S", "beta0J", "beta0A")

one_run <- function(row) {
  p <- PARAM_BASE
  p$p_poll <- row$p_poll; p$p_germ <- row$p_germ; p$p_s1 <- row$p_s1
  p$S <- row$S; p$n_founders <- as.integer(row$n_founders)
  p$beta0S <- row$beta0S; p$beta0J <- row$beta0J; p$beta0A <- row$beta0A
  p$canopy_z <- fx$ceiling
  Sys.setenv(CANOPY_K_MULT = "1")
  .rec$adults <- rep(NA_real_, TIMESTEPS); .rec$flowering <- rep(NA_real_, TIMESTEPS)
  .rec$recruits <- rep(NA_real_, TIMESTEPS)
  t0 <- Sys.time()
  setTimeLimit(elapsed = TIMEOUT_S, transient = TRUE)
  r <- tryCatch(suppressMessages(runcolonization(
    site = list(Site = row$site), niches = fx$niches, canopy_grid = fx$canopy_grid,
    microenv = fx$microenv, timesteps = TIMESTEPS, resolution = 10, carCap = 1, maxDisp = 5,
    stochastic = FALSE, Visualize = FALSE, spinup = 5,
    parameters = p, forestparams = fx$forestparams, seed = row$seed,
    clim_cache = fx$clim_cache
  )), error = function(e) structure(list(err = conditionMessage(e)), class = "runfail"))
  setTimeLimit(elapsed = Inf)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ids <- row[, ID_COLS]
  if (inherits(r, "runfail")) {
    timed_out <- grepl("elapsed time limit", r$err, fixed = TRUE)
    return(list(row = data.frame(ids, final_total = NA, final_S = NA, final_J = NA, final_A = NA,
      recruited = NA, extinct = NA, regime = if (timed_out) "timeout" else "error",
      k_bind_t = NA, exceeded_plausible_density = NA, total_K = NA,
      mean_z_S = NA, mean_z_J = NA, mean_z_A = NA,
      run_seconds = secs, error = r$err, stringsAsFactors = FALSE), summary = NULL))
  }
  n_t <- length(r$totalabundanceS)
  fS <- r$totalabundanceS[n_t]; fJ <- r$totalabundanceJ[n_t]; fA <- r$totalabundanceA[n_t]
  extinct <- all(r$totalabundanceA[(n_t %/% 2):n_t] == 0)
  recruited <- (fS + fJ) > 0
  row_out <- data.frame(ids, final_total = fS + fJ + fA, final_S = fS, final_J = fJ, final_A = fA,
    recruited = recruited, extinct = extinct, regime = r$regime,
    k_bind_t = if (is.null(r$k_bind_t)) NA_integer_ else r$k_bind_t,
    exceeded_plausible_density = isTRUE(r$exceeded_plausible_density), total_K = r$total_K,
    mean_z_S = stage_mean_z(r$abundanceS, r$size_S), mean_z_J = stage_mean_z(r$abundanceJ, r$size_J),
    mean_z_A = stage_mean_z(r$abundanceA, r$size_A),
    run_seconds = secs, error = NA_character_, stringsAsFactors = FALSE)
  summ <- data.frame(run_id = row$run_id, block = row$block, arm = row$arm, site = row$site,
    point_id = row$point_id, rep = row$rep, seed = row$seed, t = seq_len(n_t),
    totalS = r$totalabundanceS, totalJ = r$totalabundanceJ, totalA = r$totalabundanceA,
    total = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
    extinct = extinct, recruited = recruited, regime = r$regime,
    recruits_t = .rec$recruits[seq_len(n_t)], n_adults_t = .rec$adults[seq_len(n_t)],
    flowering_adults_t = .rec$flowering[seq_len(n_t)], stringsAsFactors = FALSE)
  list(row = row_out, summary = summ)
}

.write_atomic <- function(df, path) {
  tmp <- paste0(path, ".tmp")
  write.csv(df, tmp, row.names = FALSE)
  file.rename(tmp, path)
}

rows <- split(mine, seq_len(nrow(mine)))
done <- parallel::mclapply(rows, function(rr) {
  f_row <- file.path(TASK_DIR, paste0(rr$run_id, ".csv"))
  if (file.exists(f_row)) { log_msg(sprintf("  %s already done -- skipped", rr$run_id)); return(TRUE) }
  out <- tryCatch(one_run(rr), error = function(e) NULL)
  if (is.null(out)) { log_msg(sprintf("  %s harness error", rr$run_id)); return(FALSE) }
  if (!is.null(out$summary)) .write_atomic(out$summary, file.path(TASK_DIR, paste0(rr$run_id, "_summary.csv")))
  .write_atomic(out$row, f_row)
  log_msg(sprintf("  %s done: regime=%s %.0fs", rr$run_id, out$row$regime, out$row$run_seconds))
  TRUE
}, mc.cores = NCORES, mc.preschedule = FALSE)

row_files <- file.path(TASK_DIR, paste0(mine$run_id, ".csv"))
have <- file.exists(row_files)
res <- do.call(rbind, lapply(row_files[have], read.csv, stringsAsFactors = FALSE))
sum_files <- file.path(TASK_DIR, paste0(mine$run_id, "_summary.csv"))
summ <- do.call(rbind, lapply(sum_files[file.exists(sum_files)], read.csv, stringsAsFactors = FALSE))
if (!all(have)) {
  if (!is.null(res)) .write_atomic(res, paste0(out_path, ".incomplete"))
  log_msg(sprintf("task %d: INCOMPLETE -- %d of %d runs have no result (%s)", IDX, sum(!have), nrow(mine),
                  paste(mine$run_id[!have], collapse = ",")))
  quit(save = "no", status = 1)
}
.write_atomic(res, out_path)
if (!is.null(summ)) .write_atomic(summ, sum_path)
log_msg(sprintf("task %d: wrote %d rows to %s and %d summary rows to %s", IDX, nrow(res), out_path,
                if (is.null(summ)) 0L else nrow(summ), sum_path))
