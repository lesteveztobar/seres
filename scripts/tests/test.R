opts <- options(stringsAsFactors = FALSE)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg(timestamp_fmt = "long", flush_output = TRUE)

bench <- function(label, expr) {
    log_msg("Starting ", label, " run")
    t0 <- Sys.time()
    res <- tryCatch(
        eval(expr, envir = parent.frame()),
        error = function(e) {
            log_msg("ERROR in ", label, ": ", conditionMessage(e))
            stop(e)
        }
    )
    log_msg(
        "Finished ", label, " in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        " seconds"
    )
    res
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
if (length(script_file) > 0) {
    script_path <- sub("^--file=", "", script_file[1])
} else {
    stop("Cannot determine script path. Run this script with Rscript.")
}
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."))

# CORRECTED: config/paths.R now lives under scripts/02_model/config/, not
# directly under BASE_DIR -- config/ is a sibling of tests/ within 02_model,
# so it's one level up from script_dir (tests/), not three.
paths_file <- file.path(dirname(script_dir), "config", "paths.R")
if (!file.exists(paths_file)) {
    stop("config/paths.R not found at ", paths_file)
}
source(paths_file)

root_dir <- if (exists("paths") && is.list(paths) && !is.null(paths$root)) {
    normalizePath(paths$root)
} else if (exists("PROJECT_ROOT")) {
    normalizePath(PROJECT_ROOT)
} else {
    project_root
}

scripts_dir <- if (exists("paths") && is.list(paths) && !is.null(paths$scripts)) {
    file.path(root_dir, paths$scripts)
} else {
    file.path(root_dir, "scripts")
}

source_engine_files <- function() {
    engine_dir <- file.path(scripts_dir, "02_model", "engine")
    if (!dir.exists(engine_dir)) stop("Engine directory not found: ", engine_dir)
    files <- sort(list.files(engine_dir, pattern = "\\.R$", full.names = TRUE))
    for (f in files) source(f)

    extras <- c(
        file.path(scripts_dir, "02_model", "setup", "make_params.R")
    )
    for (path in extras) {
        if (!file.exists(path)) {
            stop("Required helper file not found: ", path)
        }
        source(path)
    }
}

find_microenv <- function() {
    pat <- "(la|La).*elenita.*h0\\.5.*\\.rds$"
    found <- list.files(root_dir,
        pattern = pat, recursive = TRUE,
        full.names = TRUE, ignore.case = TRUE
    )
    if (length(found) > 0) {
        return(found[[1]])
    }

    pat2 <- "elenita.*\\.rds$"
    found2 <- list.files(root_dir,
        pattern = pat2, recursive = TRUE,
        full.names = TRUE, ignore.case = TRUE
    )
    if (length(found2) > 0) {
        return(found2[[1]])
    }
    NULL
}

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    opts <- list(
        compare = FALSE,
        microenv = NULL,
        timesteps = 20L,
        spinup = 5L,
        seed = 42L,
        n_founders = 100L,
        site = "LaElenita",
        visual = FALSE,
        output = file.path(project_root, "scripts", "02_model", "tests", "test_results")
    )
    for (arg in args) {
        if (arg == "--compare") {
            opts$compare <- TRUE
        } else if (grepl("^--microenv=", arg)) {
            opts$microenv <- sub("^--microenv=", "", arg)
        } else if (grepl("^--timesteps=", arg)) {
            opts$timesteps <- as.integer(sub("^--timesteps=", "", arg))
        } else if (grepl("^--spinup=", arg)) {
            opts$spinup <- as.integer(sub("^--spinup=", "", arg))
        } else if (grepl("^--seed=", arg)) {
            opts$seed <- as.integer(sub("^--seed=", "", arg))
        } else if (grepl("^--n_founders=", arg)) {
            opts$n_founders <- as.integer(sub("^--n_founders=", "", arg))
        } else if (grepl("^--output=", arg)) opts$output <- sub("^--output=", "", arg)
    }
    opts
}

summary_xy <- function(xy) {
    values <- as.numeric(xy)
    values <- values[is.finite(values)]
    if (length(values) == 0) {
        return(NULL)
    }
    c(
        min = min(values),
        max = max(values),
        range = max(values) - min(values),
        mean = mean(values),
        sd = sd(values),
        median = median(values)
    )
}

# Reports spatial spread across pixels of each variable's MEDIAN quantile
# value, per (month, daypart) combo -- the reduced-format equivalent of what
# this used to report per raw HOUR (raw hourly per-pixel arrays no longer
# exist on disk as of the 2026-08-05 write-time quantile-reduction move; see
# get_colonization.R's .compute_voxel_quantiles()/get_clim_voxel()). Still
# answers the same underlying question ("is there real per-pixel spatial
# variance in this height's climate?"), just at the (month, daypart)
# granularity the new format actually stores rather than per raw hour.
inspect_voxel_quantiles <- function(vq, combos) {
    if (is.null(vq) || is.null(vq$quantiles)) {
        return(NULL)
    }
    cat(sprintf("  raster: %d x %d = %d pixels\n", vq$nr, vq$nc, vq$nr * vq$nc))
    for (key in intersect(combos, names(vq$quantiles))) {
        medians <- vq$quantiles[[key]][, 3] # column 3 = the 50th percentile (QUANTILE_PROBS)
        stats <- summary_xy(medians)
        cat(sprintf(
            "  %-24s min=%.4f max=%.4f range=%.4f sd=%.4f\n",
            key, stats["min"], stats["max"], stats["range"], stats["sd"]
        ))
    }
}

inspect_height <- function(height, microenv,
                           combos = c("1_day_temp", "1_night_temp", "1_day_relhum",
                                      "1_day_swdown", "annual_both_windspeed")) {
    cat("Inspecting height:", height, "\n")
    h <- load_height(microenv, height)
    if (is.null(h)) {
        cat("  height not found\n")
        return(NULL)
    }
    if (is.null(h$voxel_quantiles)) {
        cat("  no voxel_quantiles -- height file predates the 2026-08-05 quantile-reduction move, regenerate via run_microclimate_site.R\n")
        return(NULL)
    }
    inspect_voxel_quantiles(h$voxel_quantiles, combos)
    invisible(TRUE)
}

main <- function() {
    opts <- parse_args()
    log_msg("Starting test.R")
    log_msg("Loading paths from ", paths_file)
    source_engine_files()
    log_msg("Loaded engine files and helper scripts")

    microenv_path <- opts$microenv
    if (is.null(microenv_path)) microenv_path <- find_microenv()
    if (is.null(microenv_path)) stop("Cannot find La Elenita microenv RDS. Pass --microenv=/path/to/file.rds")
    log_msg("Using microenv: ", microenv_path)
    microenv <- readRDS(microenv_path)

    log_msg("Inspecting raw spatial variance for height 0.5")
    inspect_height(0.5, microenv)
    log_msg("Finished raw variance inspection")

    if (!opts$compare) {
        log_msg("Skipping comparison because --compare was not supplied")
        return(invisible(NULL))
    }

    log_msg("Loading parameters from ", file.path(PARAMS_DIR, "realistic_75founders.rds"))
    params <- readRDS(file.path(PARAMS_DIR, "realistic_75founders.rds"))
    if (!is.null(opts$n_founders)) params$n_founders <- opts$n_founders

    log_msg("Building flat-mean climate cache")
    clim_cache_avg <- build_clim_cache(microenv)
    # The per-voxel stochastic cache is NOT precomputed here (unlike
    # clim_cache_avg above): it needs the landscape's voxel->pixel footprint,
    # which only exists once init_colonization() has computed xDim/yDim/lon_min
    # etc, so runcolonization() below builds it internally per call instead.
    log_msg("Finished building climate cache")

    niches <- load_observations()
    niches <- niches[
        !is.na(niches$lat) &
            !is.na(niches$lon) &
            !is.na(niches$Height_m) &
            !is.na(niches$FinalID),
    ]

    mean_canopy <- mean(
        niches$CanopyHeight_m[niches$Area_or_Site == opts$site],
        na.rm = TRUE
    )
    canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)

    site <- list(Site = opts$site)

    # Both runs use the real engine now -- the old comparison drove a
    # separate, diverging `_spatial` prototype implementation with its own
    # scale-mismatch bug; the engine itself does real per-voxel stochastic
    # climate now (see get_colonization.R's get_clim_voxel()/
    # build_clim_cache_voxel()/.clim_voxel_slice()), so the comparison is just
    # `stochastic=FALSE` (today's deterministic flat-mean-equivalent
    # behavior) vs `stochastic=TRUE` (per-voxel/day-night quantile sampling).
    log_msg("Running deterministic colonization (stochastic=FALSE, timesteps=", opts$timesteps, ", spinup=", opts$spinup, ")")
    res_det <- bench("deterministic", quote({
        log_msg("Launching deterministic colonization")
        runcolonization(
            site, niches, canopy_grid, microenv,
            timesteps = opts$timesteps, spinup = opts$spinup,
            parameters = params, clim_cache = clim_cache_avg,
            stochastic = FALSE, seed = opts$seed
        )
    }))

    log_msg("Running stochastic colonization (stochastic=TRUE, timesteps=", opts$timesteps, ", spinup=", opts$spinup, ")")
    res_stoch <- bench("stochastic", quote({
        log_msg("Launching stochastic colonization")
        runcolonization(
            site, niches, canopy_grid, microenv,
            timesteps = opts$timesteps, spinup = opts$spinup,
            parameters = params, clim_cache = clim_cache_avg,
            stochastic = TRUE, seed = opts$seed
        )
    }))

    if (!dir.exists(opts$output)) dir.create(opts$output, recursive = TRUE)
    saveRDS(list(deterministic = res_det, stochastic = res_stoch), file = file.path(opts$output, "la_elenita_stochastic_vs_deterministic.rds"))
    log_msg("Saved comparison results to ", file.path(opts$output, "la_elenita_stochastic_vs_deterministic.rds"))
}

main()
