# get_colonization.R
# Population dynamics simulation functions for canopy colonization.
# Vital rate functions follow IPM convention (Merow et al. 2013) and are named
# after the quantity they compute. Each cites the source for its form and
# parameter values. Literature parameters are used throughout — individual size
# within each stage class is drawn from a Uniform distribution over the stage
# range (disclosed approximation; we have no individual size census data).
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

# ── Core equation primitives ──────────────────────────────────────────────────
# Each named after the quantity it computes. Used inside the IPM vital rate
# functions below. Extracting them makes the biology explicit and testable.

# Logistic survival probability for one stage given size and monthly climate.
# logit(s) = β₀_stage + β₁·z + climate_penalties
#
# Applied once per month and compounded over 12 months (run_pass3_survive_
# grow()) — survival must succeed every month, so beta0_S/J/A are calibrated
# such that p_month = p_annual^(1/12), i.e. 12 reference-condition months
# compound back to the literature annual target (S: 0.44, J: 0.60, A: 0.85
# at default beta0_A). Passing an annual-target-calibrated intercept here
# directly (logit^-1(beta0) == p_annual) would silently compound it 12x too
# many times and crash the population almost every year — see methods.tex
# for the derivation and the corresponding shifted parameter table.
#
survival_logit <- function(stage, s, temp, relhum, swdown_rel,
                           beta0_S, beta0_J, beta0_A, beta1,
                           beta_light_benefit = 0.30,
                           beta_light_stress_penalty = 0.20) {
  # pmin()/pmax() (not min()/max()) throughout: temp/relhum/swdown_rel are now
  # per-voxel arrays (see .clim_voxel_slice()), and min()/max() on a vector
  # collapse it to a single scalar instead of comparing elementwise -- pmin/
  # pmax give the same result as min/max when the inputs are plain scalars,
  # so this is unchanged behavior for any remaining scalar call site.
  light_effect <- beta_light_benefit * pmin(swdown_rel, 1.0) -
    beta_light_stress_penalty * pmax(0, swdown_rel - 1.5)
  eta <- switch(stage, # eta is the linear predictor
    S = beta0_S + beta1 * s - ((100 - relhum) / 100) + light_effect,
    J = beta0_J + beta1 * s + 0.5 * light_effect,
    A = beta0_A + beta1 * s - pmax(0, (temp - 23) * 0.05),
    stop("survival_logit: stage must be 'S', 'J', or 'A'")
  )
  return(1 / (1 + exp(-eta)))
} # here, eta is calculated depending on the stage

# Logistic stage-transition probability given annual precipitation and RH.
# logit(p) = ψ₀_stage + β_precip·P + β_rh·RH
# Sources: Zotz & Schmidt (2006), Mondragón et al. (2009), Izuddin et al. (2018)
transition_logit <- function(stage, precip_annual_mm, relhum_mean,
                             psi0S = -5.877, # S→J intercept (monthly-compounded)
                             psi0J = -5.319, # J→A intercept (monthly-compounded)
                             beta_precip = 3e-4,
                             beta_rh = 0.010) {
  psi0 <- switch(stage,
    S = psi0S,
    J = psi0J,
    stop("transition_logit: stage must be 'S' or 'J'")
  )
  return(1 / (1 + exp(-(psi0 + beta_precip * precip_annual_mm +
    beta_rh * relhum_mean))))
}


# ── Vital rate functions ──────────────────────────────────────────────────────
# Named following IPM notation: s(z,e), g(s'|s,e), p_r(s), f_s(s), p_est(e)
# s  = pseudobulb length (cm) — state variable
# e  = environmental vector [T, RH, precip, radiation] from microclimf

# ── p_r(s): Flowering probability ─────────────────────────────────────────────
# Logistic in pseudobulb size s (cm).
# Reproductive minimum PsbL ~7 cm → P(flower | z=7) ≈ 0.50.
# "Reproductive effort increased strongly with plant size" — Zotz (1998).
# Form follows Raventós (2015) and Jacquemyn et al. (2010): logistic on size.
# At s=7 cm: P ≈ 0.50 | z=10 cm: P ≈ 0.82 | z=15 cm: P ≈ 0.98
# Sources: Raventós (2015), Zotz & Schmidt (2006), Zotz (1998), Jacquemyn (2010)
flowering_prob <- function(s,
                           alpha0 = -3.5, # intercept: threshold near 7 cm
                           alpha1 = 0.5) { # size slope
  1 / (1 + exp(-(alpha0 + alpha1 * s)))
}

# ── f_s(s): Fruit production (conditional on flowering) ───────────────────────
# Poisson mean as function of pseudobulb size s (cm).
# "Larger individuals produced fruits in larger numbers" — Zotz (1998).
# Form: E[fruits | s] = exp(rho_0 + rho_1·s) — Poisson regression on size (Raventós 2015).
# At s=7 cm: ~0.9 fruits | s=12 cm: ~1.7 fruits | s=20 cm: ~4.1 fruits
# Sources: Raventós (2015) — Poisson regression of fruit number on size alone
#          Zotz (1998) — positive size effect on fruit number and fruit size
fruit_number <- function(s, rho_0 = -1.0, rho_1 = 0.12, stochastic = FALSE) {
  # rho_0 <- starting level of fecundity ; rho_1 <- how quickly fecundity rises with size
  mu <- exp(rho_0 + rho_1 * s) # exp ensures positive results
  if (stochastic) rpois(1, lambda = max(mu, 0)) else mu
}

# ── F(z', z | e): Fecundity kernel ────────────────────────────────────────────
# Full annual seed output: flowering × fruits × pollination × germination × yr-1 survival.
# p_poll: probability a flower is pollinated. Epiphytic orchids are strongly
#   pollinator-limited. "Complete pollination would raise λ to persistence threshold"
#   — Zotz & Schmidt (2006). Default 0.30 (30% of flowers pollinated).
# p_germ: probability a seed germinates. Dust seeds require mycorrhizal partner.
#   "Dispersal not limiting at landscape scale; microsite (mycorrhizal) is"
#   — McCormick & Jacquemyn (2014). Default 0.001 (McCormick & Jacquemyn 2014).
# p_s1: first-year seedling survival. "Fewer than 50% of seedlings survived
#   the first dry season" — Zotz (1998). Default 0.45.
# z is the tracked mean pseudobulb size (cm) for the cell, updated annually
# by the adult size-update step in run_pass3_survive_grow().
# Sources: Zotz (1998), Zotz & Schmidt (2006), Raventós (2015),
#          McCormick & Jacquemyn (2014)
reproduce <- function(N, s,
                      p_poll = 0.30, # pollination probability
                      p_germ = 0.001, # germination probability
                      p_s1 = 0.45, # first-year seelding survival
                      stochastic = FALSE) {
  if (is.na(N) || N == 0) {
    return(0L)
  }
  p_flower <- flowering_prob(s)
  n_fruits <- fruit_number(s, stochastic = stochastic)
  seeds <- N * p_flower * n_fruits * p_poll * p_germ * p_s1
  if (stochastic) {
    rpois(1, lambda = max(seeds, 0))
  } else {
    floor(seeds) + rbinom(1, 1, seeds - floor(seeds))
  }
}

# ── d(x'|x): Dispersal kernel ─────────────────────────────────────────────────
# Wind-mediated exponential dispersal with canopy attenuation.
# Mean dispersal distance follows Murren & Ellison (1998), adapted from the
# ballistic model of Cremer (1977) and the canopy wind-decay term of
# Cionco (1972): Uc = Uf * exp(a*(h - z)/z), where Uf is free-stream wind
# speed above canopy, h is release height, z is canopy height, and a is the
# canopy-openness coefficient (Cionco 1972: a ~ 0.02 open pine forest to
# a ~ 4 dense grass; Murren & Ellison 1998, calibrated for epiphytic orchid
# seed dispersal specifically: a = 2.14 +/- 0.155).
# CORRECTED 2026-07-28: the exponent was previously coded as -(a*height/canopy_z),
# which gives full wind speed at the ground (h=0) and maximal attenuation at
# the canopy top (h=z) -- backwards relative to both the physical process and
# Cionco (1972)/Murren & Ellison (1998) eqn 7. Corrected to a*(h - canopy_z)/canopy_z.
# `a` itself (passed in from init_colonization()) was rescaled from
# a = 23*(1 - difrac) to a = 4*difrac -- bounded to Cionco's actual 0-4
# range instead of an uncited 0-23 range, and the direction was flipped:
# denser canopy (higher diffuse radiation fraction, difrac) -> higher a
# (more wind attenuation), consistent with Cionco's own dense-canopy
# examples, rather than the previous (backwards) open-canopy-at-high-difrac
# assumption. At difrac = 0.5, a = 2, close to Murren & Ellison's a = 2.14.
# Motzer (2005), measuring wind attenuation in a southern-Ecuador montane
# forest similar to our study system, found comparatively efficient
# turbulent mixing and low wind deceleration within the canopy --
# qualitatively consistent with a value on the lower-moderate end of
# Cionco's range rather than the dense-canopy extreme (a -> 4).
# TODO: still not a TMCF-specific calibrated value (Murren & Ellison's
# a=2.14 is from a mangrove system); revisit if a wind-specific source for
# cloud forest canopies turns up. A more mechanistic alternative (LAD-based
# first-order closure, Song et al. 2021) exists but requires vertical leaf
# area density data we don't currently have -- flagged for future work,
# not adopted here.
# Sources: Murren & Ellison (1998), Cionco (1972) [wind attenuation term];
#          Motzer (2005) [qualitative TMCF wind context];
#          Winkler et al. (2009) [dispersal-fecundity tradeoff context]
# Draw all indN seeds at once — vectorized over seeds, no per-seed loop.
.ind_disperse <- function(x, y, z, indN, winddir, meanDisp, Disp, pad,
                          maxDispZ = 5) {
  wind_rad <- (winddir + 180) %% 360 * pi / 180 # wind direction in radians
  dist <- pmin(round(rexp(indN, rate = 1 / max(meanDisp, 0.1))), pad) # distance from exponential distribution
  angle <- wind_rad + runif(indN, -pi / 4, pi / 4)
  tx <- x + round(dist * sin(angle)) + pad
  ty <- y + round(dist * cos(angle)) + pad
  tz <- z + sample(-maxDispZ:maxDispZ, indN, replace = TRUE) + pad
  D <- dim(Disp)
  ok <- tx >= 1L & tx <= D[1] & ty >= 1L & ty <= D[2] & tz >= 1L & tz <= D[3]
  if (any(ok)) {
    idx <- (tx[ok] - 1L) * D[2] * D[3] + (ty[ok] - 1L) * D[3] + tz[ok]
    Disp <- Disp + array(tabulate(idx, nbins = prod(D)), dim = D)
  }
  Disp
}

disperse <- function(x, y, z, seeds, clim, height, canopy_z, a,
                     lambda = 1, Ut = 1, maxDisp = 5, maxDispZ = 5, Disp,
                     stochastic = FALSE, wind_override = NULL) {
  # wind_override: a per-voxel stochastic windspeed draw (see
  # run_pass1_disperse()'s wind_by_height / .clim_voxel_slice()), used in
  # place of clim's flat height-mean when available; winddir has no per-pixel
  # raster structure (it's a single site-wide ERA5 annual mean -- see
  # lookup_climate_by_height()), so it always stays the clim-table mean.
  wind <- if (!is.null(wind_override) && is.finite(wind_override)) wind_override else mean(clim$windspeed, na.rm = TRUE)
  winddir <- mean(clim$winddir, na.rm = TRUE)
  meanDisp <- max(1, min(round((wind * exp(a * (height - canopy_z) / canopy_z)) /
    (lambda * Ut)), maxDisp))
  .ind_disperse(x, y, z,
    indN = seeds, winddir = winddir,
    meanDisp = meanDisp, Disp = Disp, pad = maxDisp,
    maxDispZ = maxDispZ
  )
}
# ── Climate helpers ───────────────────────────────────────────────────────────

# List of usable height tiers for a microenv object, regardless of storage
# format (see load_height() below for the two formats).
microenv_heights <- function(microenv) {
  if (!is.null(microenv$.height_dir)) {
    return(sort(as.numeric(microenv$.heights)))
  }
  h_keys <- names(microenv)[!names(microenv) %in% c(".spatial", ".weather")]
  sort(as.numeric(sub("h", "", h_keys)))
}

# Load one height tier's raw tmax/tmin rasters.
# New-format microenv objects (run_microclimate_site.R) only carry a manifest
# — heights (.heights) and a directory (.height_dir) of per-height RDS files
# on scratch, saved this way because loading every height into memory at once
# needs hundreds of GB. Old-format microenv objects still embed each height's
# data directly as a list element named "h<value>".
load_height <- function(microenv, height) {
  if (!is.null(microenv$.height_dir)) {
    avail <- microenv$.heights
    h_near <- avail[which.min(abs(avail - height))]
    h_path <- file.path(microenv$.height_dir, sprintf("h%.2f.rds", h_near))
    if (!file.exists(h_path)) {
      return(NULL)
    }
    return(readRDS(h_path))
  }
  h_keys <- names(microenv)[!names(microenv) %in% c(".spatial", ".weather")]
  avail <- as.numeric(sub("h", "", h_keys))
  h_key <- h_keys[which.min(abs(avail - height))]
  microenv[[h_key]]
}

# Build the climate data frame for one height tier.
# Spatial mean across the raster at that height → one value per timestep.
#
# Two on-disk height-file formats exist, detected via presence of h$tme:
#   NEW (run_microclimate_site.R, full-year rewrite): flat list with a
#     genuine POSIXct timestamp vector h$tme (one per hourly array slice,
#     ~8760 for a full year) alongside Tz/relhum/windspeed/Rdirdown/Rdifdown.
#     No day-type split -- calendar month comes directly from h$tme.
#   OLD (pre-rewrite microenv objects still on disk / in scratch as of
#     2026-07-31): each height nests two day-type blocks, h$tmax and h$tmin,
#     produced by microclimf::subsetpointmodela(tstep="month") + runmicro()
#     -- 12 monthly representative days x 24 hourly values = 288 steps,
#     concatenated in month order (Jan..Dec), no real timestamps. Kept here
#     because every microenv_*.rds manifest that exists on disk right now
#     still points to old-format scratch files; this branch stays live until
#     they're all regenerated.
# Older-still microenv objects with just one representative day (24 steps,
# no month structure) fall back to reusing that day for every month. Macro
# variables (precip, winddir) come from microenv$.weather (ERA5 hourly),
# replicated across every row (annual mean, not month-resolved).
lookup_climate_by_height <- function(height, microenv) {
  h <- load_height(microenv, height)
  if (is.null(h)) {
    return(NULL)
  }

  .smean <- function(arr) {
    if (length(dim(arr)) == 3) {
      apply(arr, 3, mean, na.rm = TRUE)
    } else {
      rep(mean(arr, na.rm = TRUE), 24)
    }
  }

  if (!is.null(h$tme)) {
    # New format: one flat annual hourly series, real POSIXct timestamps.
    # Month is read directly off each hour's own timestamp -- no positional
    # inference, no day-type/representative-day bookkeeping needed.
    temp <- .smean(h$Tz)
    relhum <- .smean(h$relhum)
    windspeed <- .smean(h$windspeed)
    swdown <- .smean(h$Rdirdown) + .smean(h$Rdifdown)
    difrad <- .smean(h$Rdifdown)
    month <- as.integer(format(h$tme, "%m"))
    df <- data.frame(
      day_type = "annual", month = month, temp = temp, relhum = relhum,
      windspeed = windspeed, swdown = swdown, difrad = difrad
    )
  } else {
    make_df <- function(slot, day_type) {
      temp <- .smean(slot$Tz)
      relhum <- .smean(slot$relhum)
      windspeed <- .smean(slot$windspeed)
      swdown <- .smean(slot$Rdirdown) + .smean(slot$Rdifdown)
      difrad <- .smean(slot$Rdifdown)
      n <- length(temp)
      if (n %% 12 == 0 && n > 24) {
        month <- rep(1:12, each = n / 12)
      } else {
        # single representative day (old format) — reuse it for every month
        month <- rep(1:12, each = n)
        temp <- rep(temp, 12)
        relhum <- rep(relhum, 12)
        windspeed <- rep(windspeed, 12)
        swdown <- rep(swdown, 12)
        difrad <- rep(difrad, 12)
      }
      data.frame(
        day_type = day_type, month = month, temp = temp, relhum = relhum,
        windspeed = windspeed, swdown = swdown, difrad = difrad
      )
    }
    df <- rbind(make_df(h$tmax, "tmax"), make_df(h$tmin, "tmin"))
  }

  # Append macro variables from the ERA5 weather record (site-wide, not
  # height-resolved). precip is binned by calendar month from the weather
  # record's own hourly timestamps and mapped onto each row via its month
  # tag -- like temp/relhum/swdown above -- rather than collapsed to one
  # site-wide mean replicated across every row regardless of season.
  # Downstream, precip is only ever consumed as an annual total
  # (precip_annual in run_pass3_survive_grow(), which averages this column
  # back across the full year and multiplies by 8760h -- see methods.tex,
  # Stage transitions), so this changes HOW that annual figure is built (a
  # genuine seasonal average instead of one blanket mean) without changing
  # what it represents. winddir stays a single site-wide annual mean: it
  # only feeds Pass 1's dispersal step (run_pass1_disperse()), which has no
  # monthly loop, so a month-resolved wind direction would have nothing to
  # attach to.
  w <- microenv$.weather
  if (!is.null(w) && !is.null(w$obs_time)) {
    precip_by_month <- tapply(w$precip, format(w$obs_time, "%m"), mean, na.rm = TRUE)
    df$precip <- as.numeric(precip_by_month[sprintf("%02d", df$month)])
    df$winddir <- mean(w$winddir, na.rm = TRUE)
  } else if (!is.null(w)) {
    # Older microenv objects (or any weather record without a usable
    # timestamp column) -- fall back to the previous single-annual-mean
    # behavior for precip too, rather than erroring.
    df$precip <- mean(w$precip, na.rm = TRUE)
    df$winddir <- mean(w$winddir, na.rm = TRUE)
  } else {
    df$precip <- NA_real_
    df$winddir <- NA_real_
  }
  df
}

# Return the rows of clim (both tmax and tmin blocks) for a given calendar month.
lookup_climate_by_month <- function(clim, month) {
  if (is.null(clim) || nrow(clim) == 0) {
    return(NULL)
  }
  clim[clim$month == month, ]
}

# Build the full per-height (and per-height-per-month) climate lookup table
# for one microenv object. This is the same for every run against a given
# site/microenv, no matter what biological parameters are being tested — so
# callers that run many simulations against the same microenv (e.g.
# run_experiment()'s parameter sweep) should build this once up front and
# pass it into runcolonization()/init_colonization() via `clim_cache` instead
# of letting each run reload every height's climate raster from disk.
build_clim_cache <- function(microenv) {
  heights <- microenv_heights(microenv)
  # lookup_climate_by_height() reads one independent per-height RDS file from
  # scratch each call (see load_height()) -- no shared state, so this
  # parallelizes safely the same way the height-generation loop already does
  # (run_microclimate_site.R). Left sequential (mc.cores=1)
  n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))
  clim_by_height <- parallel::mclapply(heights, function(h) lookup_climate_by_height(h, microenv), mc.cores = n_cores)
  clim_month_by_height <- lapply(clim_by_height, function(clim) {
    lapply(1:12, function(m) lookup_climate_by_month(clim, m))
  })
  list(clim_by_height = clim_by_height, clim_month_by_height = clim_month_by_height)
}

# ── Per-voxel stochastic climate ("voxel" cache) ──────────────────────────────
# Companion to lookup_climate_by_height()/build_clim_cache() above. Instead of collapsing every
# hourly raster slice into ONE spatially-averaged scalar shared by every voxel
# at a height, this keeps per-pixel resolution: for each raster pixel actually
# covered by the landscape footprint, each calendar month, and each daypart
# (day/night, split per pixel per hour via Rdirdown+Rdifdown > 0 -- this
# already reflects topographic/canopy shading computed by microclimf, so it's
# more physically correct than a single global sunrise/sunset cutoff), it
# stores the empirical 5/25/50/75/95th percentiles of that pixel's pooled
# tmax+tmin hourly values. Pass functions draw a fresh stochastic value from
# these quantiles via inverse-CDF each simulated year/month
# (.sample_quantile()), instead of reading one fixed mean shared by every
# voxel at that height.
#
# Degrades to a single pooled (whole-raster) quantile table per (height,
# month, daypart, variable) when microenv$.spatial isn't usable for pixel
# lookup (.spatial_extent_usable()) -- the current state of every existing
# microenv manifest, whose .spatial$ext was saved as a raw terra::ext() S4
# object and does not survive saveRDS()/readRDS() (see run_microclimate_
# site.R's manifest fix). Real per-voxel resolution activates automatically,
# with no further code changes, once microenv is regenerated with that fix.

QUANTILE_PROBS <- c(0.05, 0.25, 0.50, 0.75, 0.95)

# TRUE if microenv$.spatial carries a plain-vector extent + raster dimensions
# (the fixed manifest format) rather than a raw terra::ext() S4 object.
.spatial_extent_usable <- function(microenv) {
  sp <- microenv$.spatial
  if (is.null(sp) || is.null(sp$ext) || is.null(sp$nrow) || is.null(sp$ncol)) {
    return(FALSE)
  }
  is.numeric(sp$ext) && length(sp$ext) == 4 && !is.null(names(sp$ext)) &&
    all(c("xmin", "xmax", "ymin", "ymax") %in% names(sp$ext))
}

# Nearest raster (row, col) for a voxel's lon/lat, using the plain-vector
# extent .spatial_extent_usable() checks for. Row 1 is the top (max lat),
# matching terra's row-major, north-up raster convention.
.lonlat_to_pixel <- function(lon, lat, sp) {
  ext <- sp$ext
  col <- floor((lon - ext[["xmin"]]) / (ext[["xmax"]] - ext[["xmin"]]) * sp$ncol) + 1L
  row <- floor((ext[["ymax"]] - lat) / (ext[["ymax"]] - ext[["ymin"]]) * sp$nrow) + 1L
  list(
    row = as.integer(pmin(pmax(row, 1L), sp$nrow)),
    col = as.integer(pmin(pmax(col, 1L), sp$ncol))
  )
}

# Empirical-quantile inverse-CDF draw -- same interpolation trick
# niche_axis_score() (below) already uses for its own lookup-table
# interpolation. `q` is a length-5 vector aligned to QUANTILE_PROBS; NA in,
# NA out (insufficient data at that voxel/month/daypart -- callers already
# handle a missing/NA climate value the same way elsewhere in this file).
# stochastic=FALSE returns the median (q[3]), i.e. today's deterministic
# behavior, so existing callers are unaffected until they opt in.
.sample_quantile <- function(q, stochastic = FALSE) {
  if (is.null(q) || anyNA(q)) {
    return(NA_real_)
  }
  if (!stochastic) {
    return(q[3])
  }
  approx(QUANTILE_PROBS, q, xout = runif(1), rule = 2)$y
}

# Per-pixel quantiles of one variable, for one (month, daypart) selection, for
# the small set of footprint pixels actually used by the landscape -- NOT the
# whole raster. `var_arr`/`day_arr` are the full (nrow, ncol, ntime) arrays --
# for a new-format height file, straight from its flat top level (h$Tz etc.);
# for an old-format one, from one of its two day-type blocks (h$tmax or
# h$tmin). `month_sel` is a logical
# mask over the time dimension for one calendar month; `want_day` selects the
# daytime (Rdirdown+Rdifdown>0), nighttime, or "both" (no daypart filter --
# used for windspeed, which isn't day/night-specific) half of those hours.
# `pixels` is a data frame of unique (row, col) pairs. Returns a named list
# keyed by "<row>_<col>", each a length-5 quantile vector (or NA x5 if that
# pixel/month/daypart has no finite observations).
.daypart_mask <- function(day_arr, daypart) {
  switch(daypart,
    day = day_arr > 0,
    night = day_arr <= 0,
    both = rep(TRUE, length(day_arr)),
    stop("daypart must be 'day', 'night', or 'both'")
  )
}

.pixel_quantiles <- function(var_arr, day_arr, month_sel, daypart, pixels) {
  out <- vector("list", nrow(pixels))
  names(out) <- paste(pixels$row, pixels$col, sep = "_")
  for (i in seq_len(nrow(pixels))) {
    r <- pixels$row[i]
    cl <- pixels$col[i]
    daypart_sel <- .daypart_mask(day_arr[r, cl, ], daypart)
    sel <- month_sel & daypart_sel
    v <- var_arr[r, cl, sel]
    v <- v[is.finite(v)]
    out[[i]] <- if (length(v) == 0) rep(NA_real_, 5) else stats::quantile(v, QUANTILE_PROBS, na.rm = TRUE, names = FALSE)
  }
  out
}

# Pooled (whole-raster) quantiles of one variable, for one (month, daypart)
# selection -- the fallback used when per-pixel resolution isn't available.
# `var_arr`/`day_arr` may be 3-D (row, col, time) or, for old single-
# representative-day microenv objects, flat (time only).
.pooled_quantiles <- function(var_arr, day_arr, month_sel, daypart) {
  if (length(dim(var_arr)) == 3) {
    time_sel <- which(month_sel)
    if (length(time_sel) == 0) {
      return(rep(NA_real_, 5))
    }
    sub_var <- var_arr[, , time_sel, drop = FALSE]
    sub_day <- day_arr[, , time_sel, drop = FALSE]
    daypart_mask <- .daypart_mask(sub_day, daypart)
    v <- sub_var[daypart_mask]
  } else {
    daypart_mask <- .daypart_mask(day_arr, daypart)
    v <- var_arr[month_sel & daypart_mask]
  }
  v <- v[is.finite(v)]
  if (length(v) == 0) rep(NA_real_, 5) else stats::quantile(v, QUANTILE_PROBS, na.rm = TRUE, names = FALSE)
}

# Per-height voxel climate cache: empirical quantiles per (month, daypart,
# variable), per-pixel when `footprint` (a data frame of unique row/col pixels
# actually covered by the landscape, from init_colonization()) is supplied and
# microenv$.spatial is usable, otherwise one pooled table for the whole
# raster. Mirrors lookup_climate_by_height()'s new-vs-old format detection
# (h$tme present -> flat full-year series with real timestamps; absent ->
# old tmax/tmin day-type blocks) but keeps day/night and (when possible)
# spatial resolution instead of collapsing them.
get_clim_voxel <- function(height, microenv, footprint = NULL) {
  h <- load_height(microenv, height)
  if (is.null(h)) {
    return(NULL)
  }

  if (!is.null(h$tme)) {
    # New format: single flat annual hourly array, real timestamps -- no
    # tmax/tmin split, so no day-type loop and no "_tmax"/"_tmin" key
    # suffix (build_clim_cache_voxel() detects this and skips combining).
    pixel_mode <- .spatial_extent_usable(microenv) && !is.null(footprint) &&
      nrow(footprint) > 0 && length(dim(h$Tz)) == 3
    quantiles <- list()
    day_arr <- h$Rdirdown + h$Rdifdown
    ntime <- if (length(dim(h$Tz)) == 3) dim(h$Tz)[3] else length(h$Tz)
    month_of_t <- as.integer(format(h$tme, "%m"))
    months_present <- c(as.character(1:12), "annual")
    swdown_arr <- h$Rdirdown + h$Rdifdown
    daynight_vars <- list(swdown = swdown_arr, difrad = h$Rdifdown, temp = h$Tz, relhum = h$relhum)
    both_vars <- list(windspeed = h$windspeed)

    for (month in months_present) {
      month_sel <- if (month == "annual") rep(TRUE, ntime) else month_of_t == as.integer(month)
      for (daypart in c("day", "night")) {
        for (vname in names(daynight_vars)) {
          if (pixel_mode) {
            pix_q <- .pixel_quantiles(daynight_vars[[vname]], day_arr, month_sel, daypart, footprint)
            for (pk in names(pix_q)) {
              quantiles[[sprintf("%s_%s_%s_%s", pk, month, daypart, vname)]] <- pix_q[[pk]]
            }
          } else {
            quantiles[[sprintf("pooled_%s_%s_%s", month, daypart, vname)]] <- .pooled_quantiles(
              daynight_vars[[vname]], day_arr, month_sel, daypart
            )
          }
        }
      }
      for (vname in names(both_vars)) {
        if (pixel_mode) {
          pix_q <- .pixel_quantiles(both_vars[[vname]], day_arr, month_sel, "both", footprint)
          for (pk in names(pix_q)) {
            quantiles[[sprintf("%s_%s_%s_%s", pk, month, "both", vname)]] <- pix_q[[pk]]
          }
        } else {
          quantiles[[sprintf("pooled_%s_%s_%s", month, "both", vname)]] <- .pooled_quantiles(
            both_vars[[vname]], day_arr, month_sel, "both"
          )
        }
      }
    }
    return(list(mode = if (pixel_mode) "pixel" else "pooled", quantiles = quantiles))
  }

  # Old format: tmax/tmin day-type blocks, no real timestamps.
  pixel_mode <- .spatial_extent_usable(microenv) && !is.null(footprint) &&
    nrow(footprint) > 0 && length(dim(h$tmax$Tz)) == 3

  quantiles <- list()

  for (day_type in c("tmax", "tmin")) {
    slot <- h[[day_type]]
    day_arr <- slot$Rdirdown + slot$Rdifdown
    n <- length(slot$Tz)
    ntime <- if (length(dim(slot$Tz)) == 3) dim(slot$Tz)[3] else n
    if (ntime %% 12 == 0 && ntime > 24) {
      month_of_t <- rep(1:12, each = ntime / 12)
    } else {
      # single representative day (old format) -- reused for every month
      month_of_t <- rep(1, ntime)
    }
    # "annual" is a 13th pseudo-month pooling every hour across the whole year
    # -- Pass 1 (dispersal) and Pass 2 (establishment) use one shared annual
    # table, matching their existing use of clim_by_height (not clim_month_
    # by_height) for those same variables.
    months_present <- c(as.character(1:12), "annual")

    for (month in months_present) {
      month_sel <- if (month == "annual") {
        rep(TRUE, ntime)
      } else if (ntime %% 12 == 0 && ntime > 24) {
        month_of_t == as.integer(month)
      } else {
        rep(TRUE, ntime) # old single-day format: every month reuses the one day
      }
      swdown_arr <- slot$Rdirdown + slot$Rdifdown
      # temp/relhum/swdown/difrad: day- and night-specific quantiles (the
      # variables that actually feed day/night-sensitive survival/establishment
      # terms). windspeed: pooled day+night only ("both") -- wind-driven seed
      # release (Pass 1 dispersal) isn't day-specific, per the plan.
      daynight_vars <- list(swdown = swdown_arr, difrad = slot$Rdifdown, temp = slot$Tz, relhum = slot$relhum)
      both_vars <- list(windspeed = slot$windspeed)

      for (daypart in c("day", "night")) {
        for (vname in names(daynight_vars)) {
          if (pixel_mode) {
            pix_q <- .pixel_quantiles(daynight_vars[[vname]], day_arr, month_sel, daypart, footprint)
            for (pk in names(pix_q)) {
              quantiles[[sprintf("%s_%s_%s_%s_%s", pk, month, daypart, vname, day_type)]] <- pix_q[[pk]]
            }
          } else {
            quantiles[[sprintf("pooled_%s_%s_%s_%s", month, daypart, vname, day_type)]] <- .pooled_quantiles(
              daynight_vars[[vname]], day_arr, month_sel, daypart
            )
          }
        }
      }
      for (vname in names(both_vars)) {
        if (pixel_mode) {
          pix_q <- .pixel_quantiles(both_vars[[vname]], day_arr, month_sel, "both", footprint)
          for (pk in names(pix_q)) {
            quantiles[[sprintf("%s_%s_%s_%s_%s", pk, month, "both", vname, day_type)]] <- pix_q[[pk]]
          }
        } else {
          quantiles[[sprintf("pooled_%s_%s_%s_%s", month, "both", vname, day_type)]] <- .pooled_quantiles(
            both_vars[[vname]], day_arr, month_sel, "both"
          )
        }
      }
    }
  }

  list(mode = if (pixel_mode) "pixel" else "pooled", quantiles = quantiles)
}

# Combine the tmax/tmin quantile entries written by get_clim_voxel() (one per
# day_type) into a single quantile per (pixel-or-pooled, month, daypart,
# variable), pooling both day-types' raw observations before taking
# quantiles would require keeping every raw value in memory; instead this
# pools the two day-types' quantile *sets* by concatenating the underlying
# probabilities and re-interpolating -- a reasonable approximation given
# tmax/tmin already both feed the same representative-day design, and this
# keeps the cache to 5 numbers/entry instead of retaining raw hourly data.
.combine_daytype_quantiles <- function(q_tmax, q_tmin) {
  if (is.null(q_tmax) || anyNA(q_tmax)) {
    return(q_tmin)
  }
  if (is.null(q_tmin) || anyNA(q_tmin)) {
    return(q_tmax)
  }
  stats::quantile(c(q_tmax, q_tmin), QUANTILE_PROBS, na.rm = TRUE, names = FALSE)
}

# Build the full voxel climate cache across every height tier. `footprint`
# (optional) is a data frame of unique (row, col) raster pixels the landscape
# actually overlaps -- see init_colonization(), which computes it from the
# landscape's lon/lat bounding box padded by maxDisp. Omitting it (the default,
# preserving existing call sites like resolution_diagnostics.R) always falls
# back to pooled (whole-raster) mode.
build_clim_cache_voxel <- function(microenv, footprint = NULL) {
  heights <- microenv_heights(microenv)
  n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))
  raw <- parallel::mclapply(heights, function(h) get_clim_voxel(h, microenv, footprint), mc.cores = n_cores)

  combined_by_height <- lapply(raw, function(entry) {
    if (is.null(entry)) {
      return(NULL)
    }
    keys_tmax <- grep("_tmax$", names(entry$quantiles), value = TRUE)
    if (length(keys_tmax) == 0) {
      # New-format height file (get_clim_voxel()'s h$tme branch): quantile
      # keys carry no tmax/tmin suffix -- there's only one day-type to begin
      # with, so nothing needs combining.
      return(list(mode = entry$mode, quantiles = entry$quantiles))
    }
    out <- list()
    for (k_tmax in keys_tmax) {
      base_key <- sub("_tmax$", "", k_tmax)
      k_tmin <- paste0(base_key, "_tmin")
      out[[base_key]] <- .combine_daytype_quantiles(entry$quantiles[[k_tmax]], entry$quantiles[[k_tmin]])
    }
    list(mode = entry$mode, quantiles = out)
  })
  names(combined_by_height) <- as.character(heights)
  list(heights = heights, clim_voxel_by_height = combined_by_height)
}

# Look up (and stochastically sample) a per-voxel-or-pooled climate value.
# `state` must carry `state$clim_cache_voxel` (see init_colonization()) and,
# for pixel mode, `state$clim_pixel_row`/`state$clim_pixel_col` (xDim x yDim
# matrices mapping each landscape voxel to its raster pixel). `zi` is the
# height-tier index; `month` is 1:12 or "annual" (Pass 1/2 use one shared
# annual table, matching their existing use of clim_by_height rather than
# clim_month_by_height); `daypart` is "day" or "night"; `var` is one of temp/
# relhum/windspeed/swdown/difrad. Returns an (xDim x yDim) matrix (pixel mode)
# or a scalar (pooled mode), broadcasting naturally in the pass functions'
# elementwise arithmetic either way.
.clim_voxel_slice <- function(state, zi, month, daypart, var, stochastic = FALSE) {
  cache <- state$clim_cache_voxel
  if (is.null(cache)) {
    return(NA_real_)
  }
  height_key <- names(cache$clim_voxel_by_height)[zi]
  entry <- cache$clim_voxel_by_height[[height_key]]
  if (is.null(entry)) {
    return(NA_real_)
  }
  month_key <- if (identical(month, "annual")) "annual" else as.character(month)

  if (entry$mode == "pooled" || is.null(state$clim_pixel_row)) {
    q <- entry$quantiles[[sprintf("pooled_%s_%s_%s", month_key, daypart, var)]]
    return(.sample_quantile(q, stochastic))
  }

  xDim <- state$xDim
  yDim <- state$yDim
  out <- matrix(NA_real_, xDim, yDim)
  pixel_keys <- paste(state$clim_pixel_row, state$clim_pixel_col, sep = "_")
  unique_keys <- unique(pixel_keys)
  # One stochastic draw per raster pixel, broadcast to every voxel sharing it
  # -- voxels sharing a pixel share the same real hourly climate record, so
  # each should draw once, not independently (they'd otherwise decorrelate
  # climate between voxels the real raster never distinguished).
  draws <- vapply(unique_keys, function(pk) {
    q <- entry$quantiles[[sprintf("%s_%s_%s_%s", pk, month_key, daypart, var)]]
    .sample_quantile(q, stochastic)
  }, numeric(1))
  out[] <- draws[pixel_keys]
  out
}

# Broadcast a .clim_voxel_slice() result (scalar in pooled mode, an xDim x
# yDim matrix in pixel mode) to a full xDim x yDim matrix, replacing any NA/
# non-finite entries (a pixel/month/daypart with no observations) with
# `fallback` -- keeps downstream arithmetic finite instead of propagating NA
# into rbinom()/survival_logit().
.fill_na <- function(mat_or_scalar, fallback, xDim, yDim) {
  if (length(mat_or_scalar) > 1) {
    out <- mat_or_scalar
    out[!is.finite(out)] <- fallback
    out
  } else {
    matrix(if (is.finite(mat_or_scalar)) mat_or_scalar else fallback, xDim, yDim)
  }
}

# ── Species climate niche ──────────────────────────────────────────────────────
# Realized niche per species, one axis per climate variable (temp, relhum,
# swdown — as many axes as variables, per prof feedback 2026-07-15). Each
# axis gets a continuous 0-100 suitability score: a kernel-density ratio of
# the species' observed ("presence") values for that variable against
# "background" — the full pooled distribution of that variable actually
# available across the landscape (every height tier x month, every site) —
# rescaled so the axis's own peak is 100.
#
# The three axis scores are combined in two steps (prof feedback
# 2026-07-15 + follow-up 2026-07-15):
#   1. Geometric mean of the three axis scores — the scale-consistent way to
#      "multiply, then normalize back to 100" (unlike dividing by a fixed
#      100^2, which just relabels the units without undoing the shrinkage:
#      three axes at 80/100 would divide down to 51.2, not stay near 80).
#      Still enforces "bad on any one axis tanks the whole score" (a 0 on
#      any axis still gives 0 overall), just without extra punishment for
#      being merely imperfect everywhere.
#   2. Rescale by a per-species, per-site "ceiling" — the best geometric-mean
#      score this species actually achieves at its own observed presence
#      heights on THIS site (niche_ceiling(), attached to each niche in
#      init_colonization() below). This exists because the three axes'
#      individual optima essentially never coincide at one real height tier
#      (temp, RH, and light don't peak together in a vertical profile), so
#      even a species' best real height might only geometric-mean to ~50-60
#      without this step — which would cap suitability well below 100
#      everywhere and risk exactly the "I will never get a suitable niche"
#      failure the multiply-then-normalize request was meant to avoid. The
#      realized niche is, by definition, where the species is actually
#      found, so those exact spots are rescaled to read ~100.
#
# Takes the already-computed heights/clim_by_height (from init_colonization,
# via build_clim_cache) rather than re-deriving climate from microenv, since
# re-reading every height tier's raster from disk here would otherwise
# duplicate the same expensive I/O.
NICHE_VARS <- c("temp", "relhum", "swdown")
NICHE_GRID_N <- 512

# Mean of variable `v` in climate table `cl`. swdown drops zero/night-time
# readings first, consistent with every other swdown use in this file
# (run_pass2_establish(), survival_logit(), .niche_matching_canopy()) — a
# mean that included every dark hour would just be diluted by a fixed
# day/night ratio rather than reflecting daytime irradiance. temp/relhum use
# every reading.
.niche_var_mean <- function(cl, v) {
  if (identical(v, "swdown")) {
    return(mean(cl[[v]][cl[[v]] > 0], na.rm = TRUE))
  }
  mean(cl[[v]], na.rm = TRUE)
}

# Kernel density of `vals` on a shared grid spanning [from, to]. density()'s
# own default bandwidth extension already tapers close to zero near both
# ends, so no separate hard threshold is needed on top of it.
.density_grid <- function(vals, from, to, n = NICHE_GRID_N) {
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2 || diff(range(vals)) == 0) {
    # Degenerate (every value identical, or a single point): fall back to a
    # one-cell spike at that value so the ratio below stays well-defined.
    x <- seq(from, to, length.out = n)
    y <- as.numeric(abs(x - mean(vals)) <= (to - from) / n)
    return(list(x = x, y = y))
  }
  d <- density(vals, from = from, to = to, n = n)
  list(x = d$x, y = d$y)
}

# Background: the pooled distribution of each climate variable across every
# voxel-month actually present in the landscape (bg_vals is a data.frame/
# list with one column per variable) — the "available but not necessarily
# occupied" reference every species' presence values are scored against.
build_background_density <- function(bg_vals, vars = NICHE_VARS, n = NICHE_GRID_N) {
  setNames(lapply(vars, function(v) {
    x <- bg_vals[[v]][is.finite(bg_vals[[v]])]
    .density_grid(x, min(x), max(x), n)
  }), vars)
}

# Per-species niche model: for each variable, the presence/background
# density ratio on the background's own grid, normalized so its own peak is
# 100. clim_vals: matrix of observed values, rows = observations, one column
# per variable in `vars`.
niche_density_model <- function(clim_vals, bg_density, vars = NICHE_VARS) {
  axes <- setNames(lapply(vars, function(v) {
    bg <- bg_density[[v]]
    pres <- .density_grid(clim_vals[, v], min(bg$x), max(bg$x), length(bg$x))
    ratio <- pres$y / pmax(bg$y, 1e-8)
    ratio[!is.finite(ratio)] <- 0
    score <- if (max(ratio) > 0) 100 * ratio / max(ratio) else rep(0, length(ratio))
    list(x = bg$x, score = score)
  }), vars)
  list(axes = axes)
}

# Mean climate per height tier (one row per height, one column per variable)
# — the per-height lookup shared by get_niche()'s this-site fallback and by
# niche_ceiling()'s per-site rescale below.
height_clim_scalars <- function(clim_by_height, vars = NICHE_VARS) {
  do.call(rbind, lapply(clim_by_height, function(cl) {
    if (is.null(cl)) {
      return(setNames(rep(NA_real_, length(vars)), vars))
    }
    vapply(vars, function(v) .niche_var_mean(cl, v), numeric(1))
  }))
}

# This-site-only fallback for a species missing from the pooled cross-site
# cache (see init_colonization()) — same density-ratio model, but scored
# against this one site's own background instead of the pooled one.
get_niche <- function(site_obs, heights, clim_by_height, bg_density, vars = NICHE_VARS) {
  clim_scalars <- height_clim_scalars(clim_by_height, vars)
  obs_clim <- function(h) clim_scalars[which.min(abs(heights - h)), , drop = TRUE]

  species_ids <- sort(unique(site_obs$FinalID))
  niches <- lapply(species_ids, function(sp) {
    obs_sp <- site_obs[site_obs$FinalID == sp & !is.na(site_obs$Height_m), ]
    if (nrow(obs_sp) == 0) {
      return(NULL)
    }
    clim_vals <- do.call(rbind, lapply(obs_sp$Height_m, obs_clim))
    niche_density_model(clim_vals, bg_density, vars)
  })
  names(niches) <- species_ids
  niches
}

# 0-100 suitability of a single variable's value on one niche axis (linear
# interpolation over the precomputed grid; constant beyond the grid's edges,
# where the score is already ~0).
niche_axis_score <- function(value, axis) {
  if (is.null(axis) || is.null(value) || is.na(value)) {
    return(NA_real_)
  }
  approx(axis$x, axis$score, xout = value, rule = 2)$y
}

# 0-100 score per axis for a set of climate values (e.g. c(temp=.., relhum=..,
# swdown=..)) under a species' niche. NULL niche (no usable observations)
# returns NULL — caller decides how to treat "no basis to score".
niche_axis_scores <- function(clim_values, niche) {
  if (is.null(niche)) {
    return(NULL)
  }
  vapply(names(niche$axes), function(v) {
    niche_axis_score(clim_values[[v]], niche$axes[[v]])
  }, numeric(1))
}

# Geometric mean of a set of 0-100 scores — the scale-consistent way to
# "multiply, then normalize back to 100" (see header note above). A 0 on any
# axis still gives 0 overall (can't take log(0)); NAs are dropped.
.geomean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(100)
  }
  if (any(x <= 0)) {
    return(0)
  }
  exp(mean(log(x)))
}

# Combined 0-100 suitability BEFORE the per-site ceiling rescale (see
# niche_ceiling()) — the geometric mean of the axis scores. A species with no
# usable observations (niche = NULL) isn't gated at all (returns 100), since
# there's no basis to restrict it.
niche_raw_score <- function(clim_values, niche) {
  scores <- niche_axis_scores(clim_values, niche)
  if (is.null(scores) || all(is.na(scores))) {
    return(100)
  }
  .geomean(scores)
}

# The best niche_raw_score() this species achieves on THIS site — the
# per-site "ceiling" that niche_overall_score() rescales against, so that
# the best-matching height reads ~100 regardless of whether the three axes'
# individual optima ever coincide at one real height tier (they usually
# don't — see header note above).
#
# Two modes, chosen by whether the species has any local observations:
#   - Species observed at this site: ceiling = best score among its OWN
#     observed presence heights (obs_clim_vals) — the Hutchinsonian
#     "realized niche is where it's actually found" framing.
#   - Species with NO observations at this site (e.g. named via
#     params$species_subset to ask "how would this species, characterized
#     elsewhere, do in a landscape it's never been recorded in?"):
#     obs_clim_vals is empty, so this falls back to landscape_clim_vals —
#     every height tier actually present at this site — and the ceiling
#     becomes "the best this landscape's own vertical profile could offer,"
#     since there's no real presence to anchor to instead. This is a
#     genuinely different, weaker claim than the observed-presence case (it
#     says nothing about whether the species could actually establish here,
#     only how the landscape's best height compares to its niche elsewhere)
#     — see init_colonization()'s species_subset log message, which reports
#     which species used which mode.
#
# obs_clim_vals / landscape_clim_vals: matrices of climate rows (rows =
# observations or height tiers, one column per NICHE_VARS), as produced by
# height_clim_scalars() + obs_clim() lookups. Falls back to 100 (no rescale)
# if neither has anything to compute a ceiling from.
niche_ceiling <- function(niche, obs_clim_vals, landscape_clim_vals = NULL) {
  candidates <- if (!is.null(obs_clim_vals) && nrow(obs_clim_vals) > 0) {
    obs_clim_vals
  } else {
    landscape_clim_vals
  }
  if (is.null(niche) || is.null(candidates) || nrow(candidates) == 0) {
    return(100)
  }
  raw <- apply(candidates, 1, function(row) niche_raw_score(as.list(row), niche))
  m <- max(raw, na.rm = TRUE)
  if (!is.finite(m) || m <= 0) 100 else m
}

# Combined 0-100 suitability AFTER the per-site ceiling rescale — the value
# actually used to gate establishment. Clamped at 100 since a voxel other
# than the species' own best-observed height can still exceed that height's
# raw score. niche$ceiling is attached per site in init_colonization() (via
# niche_ceiling() above); niches without one (e.g. evaluated outside a site
# context) fall back to no rescale.
niche_overall_score <- function(clim_values, niche) {
  raw <- niche_raw_score(clim_values, niche)
  ceiling <- if (!is.null(niche) && !is.null(niche$ceiling)) niche$ceiling else 100
  min(100, 100 * raw / ceiling)
}

# Fraction-of-1 form used to gate establishment probability (p_est <-
# p_est * niche_match(...)) — see run_pass2_establish().
niche_match <- function(clim_values, niche) {
  niche_overall_score(clim_values, niche) / 100
}

# ── Array-capable niche scoring (per-voxel) ───────────────────────────────────
# Same math as niche_axis_score()..niche_match() above, but `clim_values`'
# entries may be arrays (one climate value per voxel, e.g. from
# .clim_voxel_slice()) instead of single scalars -- used by run_pass2_
# establish() now that climate is sampled per-voxel rather than once per
# height tier. approx() itself already vectorizes over `xout`; these wrappers
# just add elementwise NA-safety (a bare `if (is.na(value))` errors on a
# vector/array of length > 1) and keep the array's dim attribute through the
# geometric-mean/rescale steps.
niche_axis_score_array <- function(value_arr, axis) {
  if (is.null(axis)) {
    return(array(NA_real_, dim = dim(value_arr)))
  }
  out <- approx(axis$x, axis$score, xout = as.vector(value_arr), rule = 2)$y
  out[!is.finite(as.vector(value_arr))] <- NA_real_
  dim(out) <- dim(value_arr)
  out
}

niche_axis_scores_array <- function(clim_values, niche) {
  if (is.null(niche)) {
    return(NULL)
  }
  lapply(names(niche$axes), function(v) niche_axis_score_array(clim_values[[v]], niche$axes[[v]]))
}

# Elementwise geometric mean across a list of same-shape score arrays --
# array analogue of .geomean() above (NAs dropped per-voxel; a voxel with all
# axes NA gets 100/no-restriction, matching .geomean()'s scalar behavior).
.geomean_array <- function(score_list) {
  d <- dim(score_list[[1]])
  stacked <- simplify2array(lapply(score_list, as.vector)) # (n_voxel, n_axis)
  if (is.null(dim(stacked))) dim(stacked) <- c(length(stacked), 1L)
  out <- apply(stacked, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(100)
    }
    if (any(x <= 0)) {
      return(0)
    }
    exp(mean(log(x)))
  })
  dim(out) <- d
  out
}

niche_raw_score_array <- function(clim_values, niche) {
  scores <- niche_axis_scores_array(clim_values, niche)
  d <- dim(clim_values[[1]])
  if (is.null(scores)) {
    return(array(100, dim = d))
  }
  .geomean_array(scores)
}

niche_overall_score_array <- function(clim_values, niche) {
  raw <- niche_raw_score_array(clim_values, niche)
  ceiling <- if (!is.null(niche) && !is.null(niche$ceiling)) niche$ceiling else 100
  pmin(100, 100 * raw / ceiling)
}

# Fraction-of-1 form, array analogue of niche_match() -- used by
# run_pass2_establish() to gate per-voxel establishment probability.
niche_match_array <- function(clim_values, niche) {
  niche_overall_score_array(clim_values, niche) / 100
}

# ── Forest structure ──────────────────────────────────────────────────────────

# Populate the landscape array with randomly placed trees and derive per-voxel
# carrying capacity from bark surface area (Myster 2017, Johansson 1974).
#
# forestparams must contain:
#   stems_per_ha         tree density for trees ≥10 cm dsh (Myster 2017: ~298/ha)
#   mean_hgt / sd_hgt    tree height distribution (m)
#   mean_crown_r / sd_crown_r  crown radius distribution (m)
#   trunk_r              mean trunk radius in m (Myster 2017: mean dsh 22.7 cm → 0.114 m)
#   branch_density       m² branch surface per m² projected crown area (literature: 2–5)
#   epiphyte_footprint_m2  bark area per individual Maxillariinae (~0.02 m²)
#
# Crown shape: bell curve peaking at 75% of tree height (widest) tapering to
# point at top — approximates tropical montane cloud forest crown architecture.
# Johansson zones 1–2 = trunk only; zones 3–5 = expanding horizontal crown.
build_forest <- function(landscape, heights, forestparams, site_obs, resolution) {
  dims <- dim(landscape)
  xDim <- dims[1]
  yDim <- dims[2]
  zDim <- dims[3]

  lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
  lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 *
    cos(mean(site_obs$lat) * pi / 180)
  area_ha <- (lat_range_m * lon_range_m) / 10000
  nTree <- max(1L, round(area_ha * forestparams$stems_per_ha))

  trees <- data.frame(
    x       = sample(1:xDim, nTree, replace = TRUE),
    y       = sample(1:yDim, nTree, replace = TRUE),
    height  = pmax(1.0, rnorm(nTree, forestparams$mean_hgt, forestparams$sd_hgt)),
    crown_r = pmax(0.5, rnorm(nTree, forestparams$mean_crown_r, forestparams$sd_crown_r))
  )
  trees$crown_r_cells <- trees$crown_r / resolution

  # voxel height thickness (m) — used for bark surface area calculation
  vox_heights <- diff(c(0, heights))

  landscape[] <- FALSE
  zone <- array(0L, dim = c(xDim, yDim, zDim))
  carCap_voxel <- array(1L, dim = c(xDim, yDim, zDim))

  for (ti in seq_len(nTree)) {
    tx <- trees$x[ti]
    ty <- trees$y[ti]
    th <- trees$height[ti]
    cr <- trees$crown_r_cells[ti]

    x_range <- max(1, tx - ceiling(cr)):min(xDim, tx + ceiling(cr))
    y_range <- max(1, ty - ceiling(cr)):min(yDim, ty + ceiling(cr))

    for (x in x_range) {
      for (y in y_range) {
        horiz_dist <- sqrt((x - tx)^2 + (y - ty)^2)

        for (z in seq_len(zDim)) {
          h <- heights[z]
          if (h > th || h < 0.5) next

          rel_h <- h / th

          jzone <- if (rel_h < 0.10) {
            1L
          } else if (rel_h < 0.30) {
            2L
          } else if (rel_h < 0.50) {
            3L
          } else if (rel_h < 0.80) {
            4L
          } else {
            5L
          }

          # Trunk zones: only the single column cell; crown zones (3-5, starting
          # at rel_h=0.30): bell-shaped radius. Domain must start at zone 3's
          # actual lower boundary (0.30), not 0.5 — using 0.5 here previously
          # made crown_fraction negative (and effective_r therefore negative,
          # i.e. never satisfied by any horiz_dist >= 0) for the entire
          # 0.30-0.50 sub-range, silently excluding that whole band from valid
          # canopy habitat on every tree.
          in_tree <- if (jzone <= 2) {
            horiz_dist == 0
          } else {
            crown_fraction <- (rel_h - 0.30) / 0.70 # 0 at zone 3 base, 1 at top
            effective_r <- cr * sin(crown_fraction * pi) # peaks at 65% height
            horiz_dist <= effective_r
          }
          if (!in_tree) next

          landscape[x, y, z] <- TRUE
          if (jzone > zone[x, y, z]) zone[x, y, z] <- jzone

          # Bark surface area (m²) per voxel → carrying capacity
          vox_h_m <- vox_heights[z]
          bark_area <- if (jzone <= 2) {
            2 * pi * forestparams$trunk_r * vox_h_m
          } else {
            crown_fraction <- (rel_h - 0.30) / 0.70
            eff_r_m <- cr * resolution * sin(crown_fraction * pi) # grid cells → m
            pi * eff_r_m^2 * forestparams$branch_density * vox_h_m
          }
          cap <- max(1L, as.integer(floor(bark_area / forestparams$epiphyte_footprint_m2)))
          if (cap > carCap_voxel[x, y, z]) carCap_voxel[x, y, z] <- cap
        }
      }
    }
  }

  log_msg(sprintf(
    "build_forest: %d trees | %.1f ha | valid voxels: %d | mean carCap: %.1f",
    nTree, area_ha, sum(landscape), mean(carCap_voxel[landscape])
  ))
  list(
    landscape = landscape, zone = zone, carCap_voxel = carCap_voxel,
    trees = trees, n_trees = nTree, area_ha = area_ha
  )
}

# ── Simulation setup ──────────────────────────────────────────────────────────

init_colonization <- function(site, niches, canopy_grid, microenv,
                              resolution = 10, carCap = 1, maxDisp = 5, params,
                              forestparams = NULL, allsites = FALSE,
                              clim_cache = NULL, clim_cache_voxel = NULL) {
  site_name <- site$Site
  heights <- microenv_heights(microenv)
  site_obs <- if (allsites) niches else niches[niches$Area_or_Site == site_name, ]
  zDim <- length(heights)
  lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
  lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 *
    cos(mean(site_obs$lat) * pi / 180)
  xDim <- max(round(lon_range_m / resolution), 10) + 4
  yDim <- max(round(lat_range_m / resolution), 10) + 4

  # Species actually observed at this site. params$species_subset (optional
  # -- see run_colonization.R's species_file arg) REPLACES this list
  # rather than narrowing it, so it can also name a species never observed
  # at this site at all -- e.g. "how would species X, characterized from
  # other sites, do in a landscape it's never been recorded in?" (a species
  # in the subset still needs *some* niche to score against -- either from
  # species_niches.rds, characterize_niches.R's cross-site cache, or from
  # this site's own get_niche() fallback below -- a species with neither is
  # simply unrestricted, same as any species with no usable observations).
  # site_obs itself is NOT filtered by the subset -- it still drives the
  # landscape's physical spatial extent (lat/lon range above) and stays
  # available for other species' per-observation lookups, so a species list
  # that adds or removes species never changes the modeled landscape itself.
  # See niche_ceiling() below for how a species with no local observations
  # at this site gets its per-site ceiling instead.
  species_ids <- sort(unique(site_obs$FinalID))
  if (!is.null(params$species_subset)) {
    requested <- sort(unique(params$species_subset))
    n_foreign <- sum(!requested %in% species_ids)
    species_ids <- requested
    log_msg(sprintf(
      "params$species_subset: modeling %d species (%d observed at this site, %d not -- scored against this landscape's best-available height instead of a local ceiling)",
      length(species_ids), length(species_ids) - n_foreign, n_foreign
    ))
  }
  n_species <- length(species_ids)
  sp_index <- setNames(seq_along(species_ids), species_ids)

  log_msg(sprintf(
    "Landscape: %d x %d x %d cells (%.0f x %.0f m, %.1f-%.1fm) | %d species",
    xDim, yDim, zDim,
    xDim * resolution, yDim * resolution,
    min(heights), max(heights), n_species
  ))

  resolution_deg <- resolution / 111000
  lon_min <- min(site_obs$lon) - resolution_deg
  lat_min <- min(site_obs$lat) - resolution_deg

  coord_to_idx <- function(lon, lat, h) {
    c(
      x = min(xDim, max(1, round((lon - lon_min) / resolution_deg) + 1)),
      y = min(yDim, max(1, round((lat - lat_min) / resolution_deg) + 1)),
      z = which.min(abs(heights - h))
    )
  }

  landscape <- array(FALSE, dim = c(xDim, yDim, zDim))

  if (!is.null(forestparams)) {
    # Tree-based landscape: place stems randomly, derive zone and carCap from geometry
    forest <- build_forest(landscape, heights, forestparams, site_obs, resolution)
    landscape <- forest$landscape
    zone <- forest$zone
    carCap_voxel <- forest$carCap_voxel
  } else {
    # Fallback: flat canopy grid ceiling, uniform carrying capacity
    cg_xDim <- nrow(canopy_grid)
    cg_yDim <- ncol(canopy_grid)
    for (x in 1:xDim) {
      for (y in 1:yDim) {
        for (z in 1:zDim) {
          cx <- min(x, cg_xDim)
          cy <- min(y, cg_yDim)
          landscape[x, y, z] <- heights[z] <= canopy_grid[cx, cy] && heights[z] >= 0.5
        }
      }
    }
    zone <- array(0L, dim = c(xDim, yDim, zDim))
    carCap_voxel <- array(as.integer(carCap), dim = c(xDim, yDim, zDim))
  }

  if (is.null(clim_cache)) {
    log_msg("Pre-computing climate lookup table from microenv...")
    clim_cache <- build_clim_cache(microenv)
    log_msg("Climate lookup ready.")
  } else {
    log_msg("Using precomputed climate lookup table.")
  }
  clim_by_height <- clim_cache$clim_by_height
  clim_month_by_height <- clim_cache$clim_month_by_height

  # Voxel -> raster-pixel mapping, for the per-voxel stochastic climate cache
  # (see .clim_voxel_slice()). Only meaningful once microenv$.spatial is a
  # usable, plain-vector extent (run_microclimate_site.R's manifest fix) --
  # on today's existing (extent-broken) manifests this stays NULL and every
  # pass function transparently falls back to pooled (whole-raster) sampling
  # instead of per-voxel (see get_clim_voxel()/build_clim_cache_voxel()).
  clim_pixel_row <- NULL
  clim_pixel_col <- NULL
  footprint <- NULL
  if (.spatial_extent_usable(microenv)) {
    sp <- microenv$.spatial
    lon_grid <- lon_min + (seq_len(xDim) - 1) * resolution_deg
    lat_grid <- lat_min + (seq_len(yDim) - 1) * resolution_deg
    lon_mat <- matrix(lon_grid, nrow = xDim, ncol = yDim)
    lat_mat <- matrix(lat_grid, nrow = xDim, ncol = yDim, byrow = TRUE)
    px <- .lonlat_to_pixel(as.vector(lon_mat), as.vector(lat_mat), sp)
    clim_pixel_row <- matrix(px$row, xDim, yDim)
    clim_pixel_col <- matrix(px$col, xDim, yDim)
    footprint <- unique(data.frame(row = as.vector(clim_pixel_row), col = as.vector(clim_pixel_col)))
    log_msg(sprintf("Per-voxel climate footprint: %d landscape voxels map to %d raster pixels", xDim * yDim, nrow(footprint)))
  }

  if (is.null(clim_cache_voxel)) {
    log_msg("Pre-computing per-voxel stochastic climate cache...")
    clim_cache_voxel <- build_clim_cache_voxel(microenv, footprint = footprint)
    log_msg(sprintf("Per-voxel climate cache ready (mode: %s).",
      if (is.null(footprint)) "pooled" else "pixel"))
  }

  valid_clim <- which(!sapply(clim_by_height, is.null))
  mid_zi <- valid_clim[which.min(abs(heights[valid_clim] -
    max(min(heights), mean(canopy_grid, na.rm = TRUE) / 2)))]
  if (length(mid_zi) == 0) mid_zi <- valid_clim[1]
  clim_mid <- clim_by_height[[mid_zi]]
  difrac <- mean(clim_mid$difrad / (clim_mid$swdown + 0.001), na.rm = TRUE)
  a <- 4 * difrac

  med_zi <- valid_clim[which.min(abs(heights[valid_clim] -
    median(heights[valid_clim])))]
  mean_swdown_site <- mean(clim_by_height[[med_zi]]$swdown, na.rm = TRUE)

  log_msg(sprintf(
    "Valid climate heights: %d/%d | mean_swdown: %.1f",
    length(valid_clim), zDim, mean_swdown_site
  ))

  # Per-species realized climate niche (temp/relhum/swdown density-ratio
  # model — see niche_density_model()/niche_match() above), gating
  # establishment in run_pass2_establish().
  #
  # Prefer the pooled cross-site cache from characterize_niches.R (every
  # observation of a species across all 5 sites, not just this one) if it
  # exists and covers every species observed at this site — a species with
  # only 2-3 records at one site may have several more elsewhere. Falls back
  # to this-site-only get_niche() for any species the cache doesn't cover
  # (e.g. added to the field data since the cache was last regenerated), or
  # entirely if the cache doesn't exist at all — scored against this site's
  # own background distribution in that case.
  # species_ids already computed above (and filtered by params$species_subset
  # if set) -- reused here rather than re-deriving from site_obs, so a
  # restricted subset also skips the niche-scoring work below for every
  # species that isn't in it.
  niche_cache_path <- file.path(PROCESSED_DIR, "species_niches.rds")
  niches_by_species <- setNames(vector("list", length(species_ids)), species_ids)

  cached_geom <- if (file.exists(niche_cache_path)) readRDS(niche_cache_path) else NULL
  missing_from_cache <- character(0)
  for (sp in species_ids) {
    if (!is.null(cached_geom) && !is.null(cached_geom[[sp]])) {
      niches_by_species[[sp]] <- cached_geom[[sp]]
    } else {
      missing_from_cache <- c(missing_from_cache, sp)
    }
  }
  if (length(missing_from_cache) > 0) {
    # Background: every height tier's monthly climate at this site (day-
    # filtered swdown, via .niche_var_mean) — same granularity as the
    # cross-site background characterize_niches.R builds, just for one site.
    site_bg_rows <- do.call(rbind, lapply(clim_month_by_height, function(by_month) {
      do.call(rbind, lapply(by_month, function(cm) {
        if (is.null(cm) || nrow(cm) == 0) {
          return(NULL)
        }
        as.data.frame(as.list(setNames(
          vapply(NICHE_VARS, function(v) .niche_var_mean(cm, v), numeric(1)), NICHE_VARS
        )))
      }))
    }))
    site_bg <- build_background_density(site_bg_rows)
    fallback <- get_niche(
      site_obs[site_obs$FinalID %in% missing_from_cache, ],
      heights, clim_by_height, site_bg
    )
    niches_by_species[missing_from_cache] <- fallback[missing_from_cache]
  }

  n_cached <- length(species_ids) - length(missing_from_cache)
  log_msg(sprintf(
    "%d/%d species have a usable observed niche (%d from cross-site cache, %d this-site-only)",
    sum(!sapply(niches_by_species, is.null)), n_species,
    n_cached, length(missing_from_cache)
  ))

  # Attach each niche's per-site ceiling (see niche_ceiling()/
  # niche_overall_score() above): the best score this species actually
  # achieves at its own observed presence heights on THIS site, so those
  # exact spots rescale to ~100 regardless of which cache the niche came
  # from. A species with no local observations here (e.g. a foreign species
  # named via params$species_subset -- see above) falls back to the best
  # score among every height tier actually present in this landscape
  # instead, so it still gets a meaningful 0-100 answer to "how well does
  # this landscape's best available height match my niche?" rather than no
  # rescale at all.
  height_scalars <- height_clim_scalars(clim_by_height)
  landscape_clim_vals <- height_scalars[stats::complete.cases(height_scalars), , drop = FALSE]
  obs_clim_at <- function(h) height_scalars[which.min(abs(heights - h)), , drop = TRUE]
  n_landscape_ceiling <- 0L
  for (sp in species_ids) {
    if (is.null(niches_by_species[[sp]])) next
    obs_sp <- site_obs[site_obs$FinalID == sp & !is.na(site_obs$Height_m), ]
    obs_clim_vals <- if (nrow(obs_sp) > 0) do.call(rbind, lapply(obs_sp$Height_m, obs_clim_at)) else NULL
    if (is.null(obs_clim_vals)) n_landscape_ceiling <- n_landscape_ceiling + 1L
    niches_by_species[[sp]]$ceiling <- niche_ceiling(niches_by_species[[sp]], obs_clim_vals, landscape_clim_vals)
  }
  if (n_landscape_ceiling > 0) {
    log_msg(sprintf(
      "%d species have no local observations at this site -- ceiling based on this landscape's best-available height instead",
      n_landscape_ceiling
    ))
  }

  params$a <- a
  params$mean_swdown_site <- mean_swdown_site
  # ensure size-tracking defaults exist if caller omitted them
  if (is.null(params$s_A_min)) params$s_A_min <- 7.0
  if (is.null(params$s_A_max)) params$s_A_max <- 20.0
  if (is.null(params$delta_s_base)) params$delta_s_base <- 0.80
  if (is.null(params$cost_repro)) params$cost_repro <- 0.50

  list(
    site_name = site_name,
    site_obs = site_obs,
    heights = heights,
    xDim = xDim, yDim = yDim, zDim = zDim,
    n_species = n_species,
    species_ids = species_ids,
    niches_by_species = niches_by_species,
    sp_index = sp_index,
    landscape = landscape,
    zone = zone,
    carCap_voxel = carCap_voxel,
    coord_to_idx = coord_to_idx,
    clim_by_height = clim_by_height,
    clim_month_by_height = clim_month_by_height,
    clim_cache_voxel = clim_cache_voxel,
    clim_pixel_row = clim_pixel_row,
    clim_pixel_col = clim_pixel_col,
    params = params,
    carCap = carCap,
    maxDisp = maxDisp,
    dispersalmatrix = array(0L, dim = c(
      xDim + 2 * maxDisp,
      yDim + 2 * maxDisp,
      zDim + 2 * maxDisp,
      n_species
    ))
  )
}

# ── Simulation passes ─────────────────────────────────────────────────────────

# Pass 1: Adults reproduce and seeds are dispersed through the 3D landscape.
# Uses tracked pseudobulb size (size_A[x,y,z,t,sp]) for the F kernel.
# Returns both the dispersal matrix and a fruited boolean array so pass 3
# can apply the cost-of-reproduction penalty to the right cells.
run_pass1_disperse <- function(state, abundanceA, size_A, t, stochastic = FALSE) {
  p <- state$params
  Disp <- array(0L, dim = c(
    state$xDim + 2 * state$maxDisp,
    state$yDim + 2 * state$maxDisp,
    state$zDim + 2 * state$maxDisp,
    state$n_species
  ))
  fruited <- array(FALSE, dim = c(state$xDim, state$yDim, state$zDim, state$n_species))
  total_seeds <- 0L

  # One windspeed draw per raster pixel per height (pooled day+night, annual
  # table -- wind-driven seed release isn't day-specific), reused by every
  # adult voxel dispersing from that height this timestep: voxels sharing a
  # pixel share the same real hourly wind record, so they draw once, not
  # independently (see .clim_voxel_slice()). Falls back to NULL (disperse()'s
  # existing clim-table mean) wherever the voxel cache has no data.
  wind_by_height <- lapply(seq_len(state$zDim), function(zi) {
    .clim_voxel_slice(state, zi, "annual", "both", "windspeed", stochastic)
  })

  for (sp in 1:state$n_species) {
    nonzero <- which(abundanceA[, , , t, sp] > 0, arr.ind = TRUE)
    for (i in seq_len(nrow(nonzero))) {
      x <- nonzero[i, 1]
      y <- nonzero[i, 2]
      z <- nonzero[i, 3]
      if (!state$landscape[x, y, z]) next
      N <- abundanceA[x, y, z, t, sp]
      if (is.na(N) || N == 0) next
      psb_s <- size_A[x, y, z, t, sp] # tracked pseudobulb size (cm)
      seeds <- reproduce(N,
        s = psb_s,
        p_poll = p$p_poll, p_germ = p$p_germ, p_s1 = p$p_s1
      )
      total_seeds <- total_seeds + seeds
      if (seeds == 0) next
      fruited[x, y, z, sp] <- TRUE
      wind_zi <- wind_by_height[[z]]
      wind_here <- if (length(wind_zi) > 1) wind_zi[x, y] else wind_zi
      Disp[, , , sp] <- disperse(
        x = x, y = y, z = z, seeds = seeds,
        clim = state$clim_by_height[[z]],
        height = state$heights[z], canopy_z = p$canopy_z,
        a = p$a, lambda = p$lambda, Ut = p$Ut,
        maxDisp = state$maxDisp, Disp = Disp[, , , sp],
        wind_override = wind_here
      )
    }
  }
  message("Pass1: seeds produced=", total_seeds, " dispersed=", sum(Disp))
  list(Disp = Disp, fruited = fruited)
}

# Pass 2: Dispersed seeds try to establish in new cells.
# Vectorized per height: p_establish is a scalar per height tier (climate is
# spatially uniform at each z), so one rbinom() call handles the full [xDim,yDim]
# slice instead of looping voxel-by-voxel.
run_pass2_establish <- function(state, abundanceS, abundanceJ, abundanceA,
                                size_S, dispersalmatrix, t, stochastic = FALSE) {
  p <- state$params
  pad <- state$maxDisp
  tnext <- t + 1L
  if (tnext > dim(abundanceS)[4]) {
    return(list(S = abundanceS, size_S = size_S))
  }

  # Total occupancy across all species — needed for carCap check
  total_occ <- apply(abundanceS[, , , t, , drop = FALSE], 1:3, sum) +
    apply(abundanceJ[, , , t, , drop = FALSE], 1:3, sum) +
    apply(abundanceA[, , , t, , drop = FALSE], 1:3, sum)
  total_seeds_seen <- 0L
  total_established <- 0L

  for (sp in seq_len(state$n_species)) {
    for (zi in seq_len(state$zDim)) {
      clim <- state$clim_by_height[[zi]]
      if (is.null(clim)) next

      # Extract the [xDim, yDim] seed-rain slice (trim padding)
      seeds_slice <- dispersalmatrix[
        (pad + 1L):(pad + state$xDim),
        (pad + 1L):(pad + state$yDim),
        zi + pad, sp
      ]
      if (sum(seeds_slice) == 0L) next

      # p_establish: per-voxel now (was a single scalar for the whole height
      # tier). Daytime quantile draw, annual table (matches Pass 1/2's
      # existing use of clim_by_height rather than the monthly table) --
      # falls back to the old flat clim-table mean wherever the voxel cache
      # has no data for a pixel/month/daypart (keeps p_est finite instead of
      # propagating NA into rbinom()). No p_germ term here — mycorrhizal
      # germination potential is already accounted for once, in Pass 1's
      # fecundity kernel (reproduce()), so every dispersed seed has already
      # "passed" that gate. This gate is purely site suitability (humidity/
      # light) for a seed that already has germination potential.
      xDim <- state$xDim
      yDim <- state$yDim
      relhum_mat <- .fill_na(
        .clim_voxel_slice(state, zi, "annual", "day", "relhum", stochastic),
        mean(clim$relhum, na.rm = TRUE), xDim, yDim
      )
      swdown_mat <- .fill_na(
        .clim_voxel_slice(state, zi, "annual", "day", "swdown", stochastic),
        mean(clim$swdown[clim$swdown > 0], na.rm = TRUE), xDim, yDim
      )
      temp_mat <- .fill_na(
        .clim_voxel_slice(state, zi, "annual", "day", "temp", stochastic),
        mean(clim$temp, na.rm = TRUE), xDim, yDim
      )
      p_est <- (relhum_mat / 100) * pmin(swdown_mat / p$mean_swdown_site, 1)

      # Gate by this species' realized climate niche (temp/relhum/swdown
      # density-ratio model — see niche_density_model()/niche_match_array()
      # in init_colonization()) — a height whose climate falls outside where
      # the species was actually observed is less likely to be colonized.
      niche_sp <- state$niches_by_species[[state$species_ids[sp]]]
      p_est <- p_est * niche_match_array(
        list(temp = temp_mat, relhum = relhum_mat, swdown = swdown_mat),
        niche_sp
      )

      # Mask: valid landscape, not at capacity, has seeds
      can_est <- state$landscape[, , zi] &
        total_occ[, , zi] < state$carCap_voxel[, , zi] &
        seeds_slice > 0L

      if (!any(can_est)) next

      total_seeds_seen <- total_seeds_seen + sum(seeds_slice[can_est])
      n <- sum(can_est)
      established <- rbinom(n, as.integer(seeds_slice[can_est]), p_est[can_est])
      # Clamp to remaining capacity
      space <- pmax(0L, state$carCap_voxel[, , zi][can_est] - total_occ[, , zi][can_est])
      established <- pmin(as.integer(established), space)
      abundanceS[, , zi, tnext, sp][can_est] <-
        abundanceS[, , zi, tnext, sp][can_est] + established
      # Newly established seedlings start at s_S_min (a freshly germinated
      # dust seed). Slot tnext is guaranteed untouched before this call --
      # nothing else writes to next year's slot earlier in the timestep --
      # so this can SET rather than blend; the accumulation step in
      # runcolonization()/run_spinup() blends this year's surviving
      # seedlings into slot tnext afterward (see methods.tex, Population
      # state and stage structure).
      newly_est <- can_est
      newly_est[can_est] <- established > 0L
      size_S[, , zi, tnext, sp][newly_est] <- p$s_S_min
      total_established <- total_established + sum(established)
    }
  }
  message("Pass2: seeds seen=", total_seeds_seen, " established=", total_established)
  list(S = abundanceS, size_S = size_S)
}

# Pass 3: Survival and stage transitions, vectorized per height tier.
#
# Survival is genuinely per-voxel for all three stages: size_S/size_J/size_A
# each track a persistent, per-voxel mean pseudobulb size (not just size_A
# as before), so survival_logit() is evaluated once per voxel using that
# voxel's own tracked size. S/J previously used a single scalar drawn fresh
# from the stage's size range every month and thrown away immediately
# (never remembered, never incremented); A previously used the SPATIAL MEAN
# size across the whole height tier rather than each voxel's own value.
# Both were replaced so that (a) a voxel's seedlings/juveniles have a real,
# growing size rather than a random stand-in, and (b) an individual
# promoted S->J or J->A carries its actual accumulated size forward instead
# of resetting to the new stage's floor (see methods.tex, Population state
# and stage structure, and Stage transitions). rbinom() still handles the
# whole [xDim x yDim] slice in one call per stage — passing a same-shape
# array of probabilities instead of a scalar is natively vectorized, so
# this doesn't reintroduce a per-voxel loop.
run_pass3_survive_grow <- function(state, abundanceS, abundanceJ, abundanceA,
                                   size_S, size_J, size_A, fruited, t, stochastic = FALSE) {
  p <- state$params
  xDim <- state$xDim
  yDim <- state$yDim

  # Helper: apply rbinom to a 2D slice with a scalar probability. Still used
  # for stage transitions, which remain climate-driven rather than
  # size-driven (see methods.tex, Stage transitions) — unchanged from before.
  .surv_slice <- function(sl, prob) {
    # prob can come out NaN (transition_logit() on a
    # precip_annual that's itself NaN, when a height tier's clim_year$precip
    # has no valid values) -- "prob <= 0" on NaN is NA, not FALSE, which
    # crashed if() with "missing value where TRUE/FALSE needed"
    # (MindoMirador/Yanayacu, 2026-07-27 and reportedly since 2026-07-16).
    # Treat an undefined probability as no transition this month, same as
    # prob <= 0, rather than propagating the NaN into rbinom() or crashing.
    if (is.na(prob) || prob <= 0 || sum(sl) == 0L) {
      return(sl * 0L)
    }
    if (prob >= 1) {
      return(sl)
    }
    array(rbinom(length(sl), as.integer(sl), prob), dim = dim(sl))
  }
  # Same idea, but for a per-voxel probability array rather than one shared
  # scalar — used for survival now that it depends on each voxel's own
  # tracked size. Cells outside the landscape are 0 already — no masking
  # needed (rbinom with size=0 is always 0 regardless of prob).
  .surv_slice_vec <- function(sl, prob_arr) {
    if (sum(sl) == 0L) {
      return(sl * 0L)
    }
    prob_arr <- pmin(pmax(prob_arr, 0), 1)
    # An undefined per-voxel probability (e.g. a pixel/month/daypart with no
    # climate observations) is treated as no transition/no death this month,
    # same convention .surv_slice() already uses for a scalar NaN prob.
    prob_arr[!is.finite(prob_arr)] <- 0
    array(rbinom(length(sl), as.integer(sl), prob_arr), dim = dim(sl))
  }

  for (sp in seq_len(state$n_species)) {
    # ── Monthly loop: survival, transitions, and growth are all evaluated per
    # month and accumulated across the year, interleaved (survive -> maybe
    # transition -> accrue growth, each month in sequence) rather than
    # survival looping monthly while transitions/growth used one annual
    # value — an individual must survive a given month to be eligible for
    # anything else that month.
    #
    # Growth's monthly increment uses the ANNUAL precip ratio (P_ann/P_ref)
    # by design — the literature-calibrated growth response is to annual
    # precipitation, not a given month's (see methods.tex, Adult size
    # dynamics) — even though precip is itself now binned by month upstream
    # (lookup_climate_by_height()) before being re-aggregated back into that annual figure.
    #
    # All three vital rates are calibrated (Table params_stable/varied) for
    # what a SINGLE evaluation should reproduce ANNUALLY, so each is
    # converted to its monthly-equivalent rate before being applied 12x:
    #   survival: p_month = p_annual^(1/12) — must succeed every month
    #     (multiplicative), so a straight 12th root recovers the annual rate.
    #   transitions: q_month = 1-(1-p_annual)^(1/12) — "at least once this
    #     year" event, so the complement (not-yet-transitioned probability)
    #     is what compounds multiplicatively across months.
    #   growth: delta_s_base/12 per month (additive, not compounded) — see
    #     size update below.
    # The corresponding intercepts (beta0_S/J/A, psi0S/J) already encode this
    # conversion; see methods.tex for the derivation.
    for (zi in seq_len(state$zDim)) {
      if (sum(abundanceS[, , zi, t, sp]) + sum(abundanceJ[, , zi, t, sp]) +
        sum(abundanceA[, , zi, t, sp]) == 0L) {
        next
      }

      clim_year <- state$clim_by_height[[zi]]
      if (is.null(clim_year)) next
      precip_annual <- mean(clim_year$precip, na.rm = TRUE) * 8760

      slS <- abundanceS[, , zi, t, sp]
      slJ <- abundanceJ[, , zi, t, sp]
      slA <- abundanceA[, , zi, t, sp]
      sS_slice <- size_S[, , zi, t, sp]
      sJ_slice <- size_J[, , zi, t, sp]
      sa_slice <- size_A[, , zi, t, sp]
      delta_s_total_S <- array(0, dim = c(xDim, yDim))
      delta_s_total_J <- array(0, dim = c(xDim, yDim))
      delta_s_total_A <- array(0, dim = c(xDim, yDim))

      for (month in 1:12) {
        cm <- state$clim_month_by_height[[zi]][[month]]
        if (is.null(cm) || nrow(cm) == 0) next

        # Per-voxel daytime quantile draws (see .clim_voxel_slice()), falling
        # back to this month's flat clim-table mean wherever the voxel cache
        # has no data for a pixel/daypart -- keeps every downstream term
        # finite instead of propagating NA into rbinom().
        temp_mat <- .fill_na(
          .clim_voxel_slice(state, zi, month, "day", "temp", stochastic),
          mean(cm$temp, na.rm = TRUE), xDim, yDim
        )
        relhum_mat <- .fill_na(
          .clim_voxel_slice(state, zi, month, "day", "relhum", stochastic),
          mean(cm$relhum, na.rm = TRUE), xDim, yDim
        )
        swdown_mat <- .fill_na(
          .clim_voxel_slice(state, zi, month, "day", "swdown", stochastic),
          mean(cm$swdown[cm$swdown > 0], na.rm = TRUE), xDim, yDim
        )
        swdown_rel <- if (p$mean_swdown_site > 0) swdown_mat / p$mean_swdown_site else matrix(1.0, xDim, yDim)
        swdown_rel[!is.finite(swdown_rel)] <- 1.0

        # ── s(z,e): monthly survival, per voxel using each voxel's own
        # tracked size AND now its own tracked climate draw (see function
        # header note above) ─────────────────────────────────────────────
        s_S_arr <- survival_logit("S", sS_slice, temp_mat, relhum_mat, swdown_rel, p$beta0_S, p$beta0_J, p$beta0_A, p$beta1)
        s_J_arr <- survival_logit("J", sJ_slice, temp_mat, relhum_mat, swdown_rel, p$beta0_S, p$beta0_J, p$beta0_A, p$beta1)
        s_A_arr <- survival_logit("A", sa_slice, temp_mat, relhum_mat, swdown_rel, p$beta0_S, p$beta0_J, p$beta0_A, p$beta1)

        slS <- .surv_slice_vec(slS, s_S_arr)
        slJ <- .surv_slice_vec(slJ, s_J_arr)
        slA <- .surv_slice_vec(slA, s_A_arr)

        # ── g(s'|s,e): monthly stage transitions (still climate-driven, not
        # size-driven — unchanged from before, now per-voxel via relhum_mat) ─
        p_StoJ <- transition_logit("S", precip_annual, relhum_mat, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
        p_JtoA <- transition_logit("J", precip_annual, relhum_mat, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
        epsilon <- rnorm(1, 0, p$sigma)
        p_StoJ <- pmin(1, pmax(0, p_StoJ + epsilon))
        p_JtoA <- pmin(1, pmax(0, p_JtoA + epsilon))

        n_StoJ <- .surv_slice_vec(slS, p_StoJ)
        n_JtoA <- .surv_slice_vec(slJ, p_JtoA)

        # Size carry-over: promoted individuals bring their CURRENT tracked
        # size with them instead of resetting to the destination stage's
        # floor, blended (count-weighted) with whatever's already in that
        # stage this voxel. J->A promotions are floored at s_A_min (the
        # adult stage's defined lower bound) in case a transition fires
        # while the juvenile's tracked size is still below it. Order
        # matters: A's blend must use sJ_slice's value BEFORE J's blend
        # below overwrites it. Explicit mask + index (rather than ifelse())
        # to match this file's existing style and avoid any doubt about
        # dim-attribute preservation on a matrix.
        sa_denom <- slA + n_JtoA
        sa_mask <- sa_denom > 0L
        sa_blend <- (slA * sa_slice + n_JtoA * pmax(p$s_A_min, sJ_slice)) / pmax(sa_denom, 1)
        sa_slice[sa_mask] <- sa_blend[sa_mask]

        J_stayers <- slJ - n_JtoA
        sJ_denom <- J_stayers + n_StoJ
        sJ_mask <- sJ_denom > 0L
        sJ_blend <- (J_stayers * sJ_slice + n_StoJ * sS_slice) / pmax(sJ_denom, 1)
        sJ_slice[sJ_mask] <- sJ_blend[sJ_mask]

        slS <- slS - n_StoJ
        slJ <- slJ + n_StoJ - n_JtoA
        slA <- slA + n_JtoA

        # ── Growth increment: this month's share, own noise draw per stage.
        # Deferred and applied once at year end (below), exactly like
        # adults already worked — only the transition-time carry-over above
        # needs to happen progressively within the year, since a newly-
        # promoted individual needs *some* valid tracked size immediately.
        # sigma/sqrt(12) keeps total annual variance matched to the original
        # (pre-monthly) calibration, since variances of independent draws
        # sum. All three stages currently share delta_s_base (literature-
        # calibrated for ADULTS specifically, Zotz 1998) — a disclosed
        # simplification in the absence of stage-specific growth-rate
        # literature values for seedlings/juveniles; see methods.tex.
        ds_base <- (p$delta_s_base / 12) * (precip_annual / 2500) * (relhum_mat / 85)
        noise_sd <- p$sigma * 0.5 / sqrt(12)
        ds_S <- ds_base + rnorm(xDim * yDim, 0, noise_sd)
        dim(ds_S) <- c(xDim, yDim)
        ds_J <- ds_base + rnorm(xDim * yDim, 0, noise_sd)
        dim(ds_J) <- c(xDim, yDim)
        ds_A <- ds_base + rnorm(xDim * yDim, 0, noise_sd)
        dim(ds_A) <- c(xDim, yDim)
        delta_s_total_S <- delta_s_total_S + ds_S
        delta_s_total_J <- delta_s_total_J + ds_J
        delta_s_total_A <- delta_s_total_A + ds_A
      }

      abundanceS[, , zi, t, sp] <- slS
      abundanceJ[, , zi, t, sp] <- slJ
      abundanceA[, , zi, t, sp] <- slA

      # Cost of reproduction: reduce the year's total growth where the cell
      # fruited (adults only — S/J don't reproduce, so this doesn't apply
      # to their growth totals).
      fr_slice <- fruited[, , zi, sp]
      delta_s_total_A[fr_slice] <- delta_s_total_A[fr_slice] * p$cost_repro

      sS_new <- pmin(p$s_S_max, pmax(p$s_S_min, sS_slice + delta_s_total_S))
      sJ_new <- pmin(p$s_J_max, pmax(p$s_J_min, sJ_slice + delta_s_total_J))
      sA_new <- pmin(p$s_A_max, pmax(p$s_A_min, sa_slice + delta_s_total_A))
      # NOTE: an individual promoted mid-year still accrues this whole
      # year's growth total on top of its carried-over starting size,
      # rather than only its fraction of the year since promotion — a
      # small, disclosed over-count traded off here against the complexity
      # of tracking which month each voxel's promotion happened in.

      # Only write back where that stage still has individuals this year,
      # so a tracked size from a previous year isn't lost while a voxel is
      # temporarily unoccupied at that stage.
      has_S <- abundanceS[, , zi, t, sp] > 0L
      has_J <- abundanceJ[, , zi, t, sp] > 0L
      has_A <- abundanceA[, , zi, t, sp] > 0L
      sS_slice[has_S] <- sS_new[has_S]
      sJ_slice[has_J] <- sJ_new[has_J]
      sa_slice[has_A] <- sA_new[has_A]
      size_S[, , zi, t, sp] <- sS_slice
      size_J[, , zi, t, sp] <- sJ_slice
      size_A[, , zi, t, sp] <- sa_slice
    }
  }
  list(
    S = abundanceS, J = abundanceJ, A = abundanceA,
    size_S = size_S, size_J = size_J, size_A = size_A
  )
}

# ── Internal live visualisation ───────────────────────────────────────────────

.plot_live <- function(state, abundanceS, abundanceJ, abundanceA,
                       totalS, totalJ, totalA, t, carCap, sleeptime = 0.2) {
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols <- if (state$n_species == 1) {
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  } else {
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)
  }
  total <- totalS + totalJ + totalA
  n_sp <- state$n_species
  par(mfrow = c(1, n_sp + 1), mar = c(4, 4, 3, 2))
  for (sp in 1:n_sp) {
    ts_S <- sapply(1:t, function(i) sum(abundanceS[, , , i, sp]))
    ts_J <- sapply(1:t, function(i) sum(abundanceJ[, , , i, sp]))
    ts_A <- sapply(1:t, function(i) sum(abundanceA[, , , i, sp]))
    ts_total <- ts_S + ts_J + ts_A
    plot(ts_total,
      type = "b", col = sp_cols[sp], lwd = 2,
      ylim = c(0, max(ts_total, 1)), xlab = "Year", ylab = "Abundance",
      main = paste0(state$species_ids[sp], " (t=", t, ")"), las = 1
    )
    lines(ts_S, type = "b", col = stage_cols[1], pch = 16, lty = 2)
    lines(ts_J, type = "b", col = stage_cols[2], pch = 17, lty = 2)
    lines(ts_A, type = "b", col = stage_cols[3], pch = 15, lty = 2)
    legend("topleft",
      legend = c("Total", "S", "J", "A"),
      col = c(sp_cols[sp], stage_cols), lty = c(1, 2, 2, 2),
      pch = c(NA, 16, 17, 15), cex = 0.6
    )
  }
  plot(total[1:t],
    type = "b", col = "black", lwd = 2,
    ylim = c(0, max(total, 1)), xlab = "Year", ylab = "Abundance",
    main = paste0("All species (t=", t, ")"), las = 1
  )
  lines(totalS[1:t], type = "b", col = stage_cols[1], pch = 16)
  lines(totalJ[1:t], type = "b", col = stage_cols[2], pch = 17)
  lines(totalA[1:t], type = "b", col = stage_cols[3], pch = 15)
  abline(h = carCap * state$xDim * state$yDim * state$zDim, col = "red", lty = 2)
  legend("topleft",
    legend = c("Total", "S", "J", "A"),
    col = c("black", stage_cols), lty = 1, pch = c(NA, 16, 17, 15), cex = 0.6
  )
  dev.flush()
  Sys.sleep(sleeptime)
}

# ── Post-hoc abundance plot ───────────────────────────────────────────────────

plot_abundance <- function(result, t = NULL, species_specific = TRUE) {
  state <- result$state
  t_max <- if (is.null(t)) length(result$totalabundanceA) else t
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols <- if (state$n_species == 1) {
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  } else {
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)
  }
  totalS <- result$totalabundanceS
  totalJ <- result$totalabundanceJ
  totalA <- result$totalabundanceA
  total <- totalS + totalJ + totalA
  n_panels <- if (species_specific) state$n_species + 1L else 1L
  par(mfrow = c(1, n_panels), mar = c(4, 4, 3, 2))
  if (species_specific) {
    for (sp in 1:state$n_species) {
      ts_S <- sapply(1:t_max, function(i) sum(result$abundanceS[, , , i, sp]))
      ts_J <- sapply(1:t_max, function(i) sum(result$abundanceJ[, , , i, sp]))
      ts_A <- sapply(1:t_max, function(i) sum(result$abundanceA[, , , i, sp]))
      ts_total <- ts_S + ts_J + ts_A
      plot(ts_total,
        type = "b", col = sp_cols[sp], lwd = 2,
        ylim = c(0, max(ts_total, 1)), xlab = "Year", ylab = "Abundance",
        main = paste0(state$species_ids[sp], " (t=", t_max, ")"), las = 1
      )
      lines(ts_S, type = "b", col = stage_cols[1], pch = 16, lty = 2)
      lines(ts_J, type = "b", col = stage_cols[2], pch = 17, lty = 2)
      lines(ts_A, type = "b", col = stage_cols[3], pch = 15, lty = 2)
      legend("topleft",
        legend = c("Total", "S", "J", "A"),
        col = c(sp_cols[sp], stage_cols), lty = c(1, 2, 2, 2),
        pch = c(NA, 16, 17, 15), cex = 0.6
      )
    }
  }
  plot(total[1:t_max],
    type = "b", col = "black", lwd = 2,
    ylim = c(0, max(total[1:t_max], 1)), xlab = "Year", ylab = "Abundance",
    main = paste0(state$site_name, " — All species (t=", t_max, ")"), las = 1
  )
  lines(totalS[1:t_max], type = "b", col = stage_cols[1], pch = 16)
  lines(totalJ[1:t_max], type = "b", col = stage_cols[2], pch = 17)
  lines(totalA[1:t_max], type = "b", col = stage_cols[3], pch = 15)
  abline(h = state$carCap * state$xDim * state$yDim * state$zDim, col = "red", lty = 2)
  legend("topleft",
    legend = c("Total", "S", "J", "A"),
    col = c("black", stage_cols), lty = 1, pch = c(NA, 16, 17, 15), cex = 0.6
  )
}

# ── 3D post-hoc visualisation ─────────────────────────────────────────────────

# Translucent point-cloud underlay shared by plot_3d_abundance() and
# plot_3d_abundance_animated(): state$landscape is the full boolean canopy-
# occupancy array already computed for the run, reused here purely as a
# visual backdrop so abundance markers read as embedded in the canopy rather
# than floating in empty space. Subsampled (default cap 20,000 voxels) --
# this is a shape cue, not a faithful full-resolution render, and plotly gets
# slow/heavy (especially the animated HTML export) well before every valid
# voxel is actually needed to convey "there is canopy here."
.canopy_context_trace <- function(state, max_points = 20000, seed = 1) {
  idx <- which(state$landscape, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(NULL)
  }
  if (nrow(idx) > max_points) {
    set.seed(seed)
    idx <- idx[sample.int(nrow(idx), max_points), , drop = FALSE]
  }
  data.frame(x = idx[, 1], y = idx[, 2], z = idx[, 3])
}

plot_3d_abundance <- function(result, t = NULL, show_canopy = TRUE,
                              canopy_opacity = 0.05, canopy_max_points = 20000) {
  state <- result$state
  if (is.null(t)) t <- dim(result$abundanceA)[4]
  sp_cols <- if (state$n_species == 1) {
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  } else {
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)
  }
  rows <- list()
  for (sp in 1:state$n_species) {
    sp_name <- state$species_ids[sp]
    col <- sp_cols[sp]
    for (stg in list(
      list(arr = result$abundanceS, nm = "S", sym = "circle"),
      list(arr = result$abundanceJ, nm = "J", sym = "diamond"),
      list(arr = result$abundanceA, nm = "A", sym = "square")
    )) {
      idx <- which(stg$arr[, , , t, sp] > 0, arr.ind = TRUE)
      if (nrow(idx) > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          x = idx[, 1], y = idx[, 2], z = idx[, 3],
          species = sp_name, stage = stg$nm, color = col, symbol = stg$sym,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) {
    message("No individuals to plot at t=", t)
    return(invisible(NULL))
  }
  df <- do.call(rbind, rows)
  traces <- split(df, paste0(df$species, "_", df$stage))
  fig <- plotly::plot_ly()
  if (show_canopy) {
    canopy_df <- .canopy_context_trace(state, max_points = canopy_max_points)
    if (!is.null(canopy_df)) {
      fig <- plotly::add_trace(fig,
        data = canopy_df, x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        name = "Canopy", showlegend = TRUE,
        marker = list(
          color = "#6b4423", size = 2,
          opacity = canopy_opacity
        )
      )
    }
  }
  for (tr in traces) {
    fig <- plotly::add_trace(fig,
      data = tr, x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      name = paste0(tr$species[1], " ", tr$stage[1]),
      marker = list(
        symbol = tr$symbol[1], color = tr$color[1],
        size = 6, opacity = 0.85
      )
    )
  }
  fig <- plotly::layout(fig,
    title = paste0("Abundance (t=", t, ")"),
    scene = list(
      xaxis = list(title = "x"), yaxis = list(title = "y"),
      zaxis = list(title = "height tier")
    )
  )
  print(fig)
  invisible(fig)
}

# Same idea as plot_3d_abundance() but across every timestep, using plotly's
# built-in frame/animation support (play button + slider) instead of a
# single static scatter. Saved as a self-contained HTML if out_path is
# given. Not yet run against real output — verify once you have a result
# worth animating (e.g. from a best_case replicate that actually persists).
# The canopy-context trace (show_canopy=TRUE default, added 2026-07-24) in
# particular needs a visual check: mixing an unframed static trace with a
# framed animated one in the same figure is standard plotly behavior, but
# wasn't exercised against a real result in this session -- open the saved
# HTML once and confirm the canopy points stay put while the slider moves.
plot_3d_abundance_animated <- function(result, out_path = NULL, show_canopy = TRUE,
                                       canopy_opacity = 0.05, canopy_max_points = 20000) {
  state <- result$state
  n_t <- dim(result$abundanceA)[4]
  sp_cols <- if (state$n_species == 1) {
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  } else {
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)
  }

  rows <- list()
  for (t in seq_len(n_t)) {
    for (sp in seq_len(state$n_species)) {
      sp_name <- state$species_ids[sp]
      col <- sp_cols[sp]
      for (stg in list(
        list(arr = result$abundanceS, nm = "S", sym = "circle"),
        list(arr = result$abundanceJ, nm = "J", sym = "diamond"),
        list(arr = result$abundanceA, nm = "A", sym = "square")
      )) {
        idx <- which(stg$arr[, , , t, sp] > 0, arr.ind = TRUE)
        if (nrow(idx) > 0) {
          rows[[length(rows) + 1]] <- data.frame(
            x = idx[, 1], y = idx[, 2], z = idx[, 3], t = t,
            species = sp_name, stage = stg$nm, color = col,
            symbol = stg$sym, trace = paste0(sp_name, " ", stg$nm),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(rows) == 0) {
    message("No individuals to plot across any timestep")
    return(invisible(NULL))
  }
  df <- do.call(rbind, rows)

  # trace -> color lookup (one row per unique trace, in matching order —
  # safer than pairing two independently-deduplicated vectors)
  trace_lu <- df[!duplicated(df$trace), c("trace", "color")]
  colors_named <- setNames(trace_lu$color, trace_lu$trace)

  # Canopy trace is added first, with no `frame` mapping, so it renders as a
  # static backdrop that persists unchanged across every animation frame
  # (plotly supports mixing framed and unframed traces in one figure) --
  # only the abundance trace below actually animates by year.
  fig <- plotly::plot_ly()
  if (show_canopy) {
    canopy_df <- .canopy_context_trace(state, max_points = canopy_max_points)
    if (!is.null(canopy_df)) {
      fig <- plotly::add_trace(fig,
        data = canopy_df, x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        name = "Canopy", showlegend = TRUE,
        marker = list(
          color = "#6b4423", size = 2,
          opacity = canopy_opacity
        )
      )
    }
  }
  fig <- plotly::add_trace(
    fig,
    data = df, x = ~x, y = ~y, z = ~z, frame = ~t, color = ~trace,
    colors = colors_named,
    symbol = ~symbol, symbols = c(circle = "circle", diamond = "diamond", square = "square"),
    type = "scatter3d", mode = "markers",
    marker = list(size = 6, opacity = 0.85)
  )
  fig <- fig |>
    plotly::layout(
      title = "Abundance over time",
      scene = list(
        xaxis = list(title = "x"), yaxis = list(title = "y"),
        zaxis = list(title = "height tier")
      )
    ) |>
    plotly::animation_opts(frame = 400, transition = 200, redraw = TRUE) |>
    plotly::animation_slider(currentvalue = list(prefix = "Year: "))

  if (!is.null(out_path)) {
    # selfcontained=TRUE needs pandoc (not installed on the cluster); FALSE
    # writes a small "<name>_files/" dependency folder alongside the HTML
    # instead -- keep the two together when copying/viewing elsewhere.
    htmlwidgets::saveWidget(fig, out_path, selfcontained = FALSE)
    message("Saved: ", out_path)
  }
  invisible(fig)
}

# ── Spin-up ───────────────────────────────────────────────────────────────────

# Founders: n_founders individuals per species (a fixed, decoupled count —
# not tied to however many field observations that species happens to have),
# placed at random canopy voxels whose local climate matches that species'
# realized niche (see niche_density_model()/niche_match()). Previously
# founders were placed only at the exact voxel of an observed individual,
# which failed whenever the stochastic forest didn't happen to mark that
# exact voxel as canopy — with a handful of observations per site, that
# could (and did) place zero founders across every replicate, guaranteeing
# extinction before the simulation even started, independent of any
# vital-rate parameter.
#
# "Matches" here means overall score > 0, not the historical exact-1.0
# match — with a continuous density-ratio score across three axes, a voxel
# essentially never hits the joint maximum on all three at once (that would
# require every axis to peak simultaneously at that height), so requiring
# score >= 1 would systematically exclude almost every voxel and defeat the
# whole point of moving off the old hard-threshold box. Any voxel with a
# nonzero score is inside the plausible climate envelope on every axis.
# Falls back to unrestricted canopy placement if no voxel matches the niche
# (species with too few observations to build one, or a niche too narrow for
# this stochastic forest) so founder placement never silently fails.
.niche_matching_canopy <- function(state, niche_sp) {
  zDim <- state$zDim
  height_ok <- vapply(seq_len(zDim), function(zi) {
    clim <- state$clim_by_height[[zi]]
    if (is.null(clim)) {
      return(FALSE)
    }
    niche_match(
      list(
        temp = mean(clim$temp, na.rm = TRUE),
        relhum = mean(clim$relhum, na.rm = TRUE),
        swdown = mean(clim$swdown[clim$swdown > 0], na.rm = TRUE)
      ),
      niche_sp
    ) > 0
  }, logical(1))

  candidate_list <- lapply(which(height_ok), function(zi) {
    xy <- which(state$landscape[, , zi], arr.ind = TRUE)
    if (nrow(xy) == 0) {
      return(NULL)
    }
    cbind(xy, z = zi)
  })
  do.call(rbind, Filter(Negate(is.null), candidate_list))
}

run_spinup <- function(state, n_gens = 5, Visualize = TRUE,
                       carCap = 1, sleeptime = 0.2, visualize_dispersion = FALSE,
                       stochastic = FALSE) {
  p <- state$params
  xDim <- state$xDim
  yDim <- state$yDim
  zDim <- state$zDim
  n_species <- state$n_species
  spinupS <- array(0L, dim = c(xDim, yDim, zDim, n_gens, n_species))
  spinupJ <- array(0L, dim = c(xDim, yDim, zDim, n_gens, n_species))
  spinupA <- array(0L, dim = c(xDim, yDim, zDim, n_gens, n_species))
  # size_S/size_J/size_A: mean pseudobulb length (cm) per seedling/juvenile/
  # adult cell; initialised at each stage's own floor.
  size_S <- array(p$s_S_min, dim = c(xDim, yDim, zDim, n_gens, n_species))
  size_J <- array(p$s_J_min, dim = c(xDim, yDim, zDim, n_gens, n_species))
  size_A <- array(p$s_A_min, dim = c(xDim, yDim, zDim, n_gens, n_species))
  totalS <- numeric(n_gens)
  totalJ <- numeric(n_gens)
  totalA <- numeric(n_gens)

  n_founders <- if (!is.null(p$n_founders)) p$n_founders else 30

  for (sp in seq_len(n_species)) {
    sp_name <- state$species_ids[sp]
    niche_sp <- state$niches_by_species[[sp_name]]
    candidates <- .niche_matching_canopy(state, niche_sp)
    if (is.null(candidates) || nrow(candidates) == 0) {
      log_msg(sprintf(
        "Spin-up: no niche-matching canopy for %s, falling back to unrestricted placement", sp_name
      ))
      candidates <- which(state$landscape, arr.ind = TRUE)
    }
    if (nrow(candidates) == 0) {
      log_msg(sprintf("Spin-up: no canopy at all for %s -- 0 founders placed", sp_name))
      next
    }

    chosen <- candidates[sample.int(nrow(candidates), n_founders,
      replace = nrow(candidates) < n_founders
    ), , drop = FALSE]
    for (k in seq_len(nrow(chosen))) {
      x <- chosen[k, 1]
      y <- chosen[k, 2]
      z <- chosen[k, 3]
      spinupA[x, y, z, 1, sp] <- spinupA[x, y, z, 1, sp] + 1L
      size_A[x, y, z, 1, sp] <- p$s_A_min
    }
  }
  log_msg(sprintf(
    "Spin-up: placed %d founders/species x %d species = %d total",
    n_founders, n_species, n_founders * n_species
  ))

  fruited <- array(FALSE, dim = c(xDim, yDim, zDim, n_species))
  Disp <- NULL

  for (gen in 1:n_gens) {
    log_msg(sprintf("Spin-up generation %d/%d", gen, n_gens))
    pass1 <- run_pass1_disperse(state, spinupA, size_A, gen, stochastic = stochastic)
    Disp <- pass1$Disp
    fruited <- pass1$fruited
    if (visualize_dispersion) {
      fields::image.plot(apply(Disp, c(1, 2), sum),
        col = scico::scico(25, palette = "lajolla"),
        main = paste0("Dispersal (gen ", gen, ") | seeds: ", sum(Disp)),
        xlab = "x", ylab = "y"
      )
      dev.flush()
      Sys.sleep(sleeptime)
    }
    pass2 <- run_pass2_establish(state, spinupS, spinupJ, spinupA, size_S, Disp, gen, stochastic = stochastic)
    spinupS <- pass2$S
    size_S <- pass2$size_S
    result <- run_pass3_survive_grow(
      state, spinupS, spinupJ, spinupA,
      size_S, size_J, size_A, fruited, gen,
      stochastic = stochastic
    )
    spinupS <- result$S
    spinupJ <- result$J
    spinupA <- result$A
    size_S <- result$size_S
    size_J <- result$size_J
    size_A <- result$size_A
    if (gen < n_gens) {
      # size_S needs a count-weighted blend, not a plain carry-forward: slot
      # gen+1 already holds this generation's newly-established seedlings
      # (from run_pass2_establish() above, size s_S_min) BEFORE the
      # surviving/grown population from slot gen is added on top -- see
      # methods.tex, Population state and stage structure. size_J/size_A
      # have no such second contributor here (nothing establishes directly
      # into J or A), so they still just carry forward.
      n_new <- spinupS[, , , gen + 1, ]
      n_carry <- spinupS[, , , gen, ]
      tot_n <- n_new + n_carry
      blended <- (n_new * size_S[, , , gen + 1, ] + n_carry * size_S[, , , gen, ]) / pmax(tot_n, 1)
      blended[tot_n == 0] <- p$s_S_min
      size_S[, , , gen + 1, ] <- blended

      spinupS[, , , gen + 1, ] <- spinupS[, , , gen + 1, ] + spinupS[, , , gen, ]
      spinupJ[, , , gen + 1, ] <- spinupJ[, , , gen + 1, ] + spinupJ[, , , gen, ]
      spinupA[, , , gen + 1, ] <- spinupA[, , , gen + 1, ] + spinupA[, , , gen, ]
      size_J[, , , gen + 1, ] <- size_J[, , , gen, ] # carry size forward
      size_A[, , , gen + 1, ] <- size_A[, , , gen, ] # carry size forward
    }
    totalS[gen] <- sum(spinupS[, , , gen, ])
    totalJ[gen] <- sum(spinupJ[, , , gen, ])
    totalA[gen] <- sum(spinupA[, , , gen, ])
    if (Visualize) {
      .plot_live(
        state, spinupS, spinupJ, spinupA,
        totalS, totalJ, totalA, gen, carCap, sleeptime
      )
    }
    log_msg(sprintf(
      "Gen %d: S=%d J=%d A=%d total=%d", gen,
      totalS[gen], totalJ[gen], totalA[gen],
      totalS[gen] + totalJ[gen] + totalA[gen]
    ))
  }
  log_msg(sprintf(
    "Spin-up complete: S=%d J=%d A=%d",
    sum(spinupS[, , , n_gens, ]), sum(spinupJ[, , , n_gens, ]),
    sum(spinupA[, , , n_gens, ])
  ))
  list(
    S = spinupS[, , , n_gens, ], J = spinupJ[, , , n_gens, ], A = spinupA[, , , n_gens, ],
    size_S = size_S[, , , n_gens, ], size_J = size_J[, , , n_gens, ],
    size_A = size_A[, , , n_gens, ], last_disp = Disp
  )
}

# ── Main wrapper ──────────────────────────────────────────────────────────────

runcolonization <- function(site, niches, canopy_grid, microenv,
                            timesteps = 50, resolution = 10, carCap = 1,
                            maxDisp = 5, stochastic = FALSE,
                            Visualize = TRUE, sleeptime = 0.2,
                            visualize_dispersion = FALSE, spinup = 5,
                            parameters, forestparams = NULL, allsites = FALSE,
                            train_frac = 0.70, seed = 42, clim_cache = NULL,
                            clim_cache_voxel = NULL) {
  set.seed(seed)

  # Per-species train/validation split — train_frac of observations per species
  # are used for spin-up; the remainder are held out for validation.
  site_obs_all <- if (allsites) niches else niches[niches$Area_or_Site == site, ]
  train_idx <- unlist(lapply(
    split(
      seq_len(nrow(site_obs_all)),
      site_obs_all$FinalID
    ),
    function(idx) sample(idx, max(1L, round(length(idx) * train_frac)))
  ))
  niches_train <- site_obs_all[train_idx, ]
  niches_val <- site_obs_all[-train_idx, ]
  log_msg(sprintf(
    "Train/val split: %d train | %d val observations (train_frac=%.2f)",
    nrow(niches_train), nrow(niches_val), train_frac
  ))

  state <- init_colonization(site, niches_train, canopy_grid, microenv,
    resolution, carCap, maxDisp,
    params = parameters,
    forestparams = forestparams, allsites = allsites,
    clim_cache = clim_cache, clim_cache_voxel = clim_cache_voxel
  )
  xDim <- state$xDim
  yDim <- state$yDim
  zDim <- state$zDim
  n_species <- state$n_species
  abundanceS <- array(0L, dim = c(xDim, yDim, zDim, timesteps, n_species))
  abundanceJ <- array(0L, dim = c(xDim, yDim, zDim, timesteps, n_species))
  abundanceA <- array(0L, dim = c(xDim, yDim, zDim, timesteps, n_species))
  size_S <- array(parameters$s_S_min, dim = c(xDim, yDim, zDim, timesteps, n_species))
  size_J <- array(parameters$s_J_min, dim = c(xDim, yDim, zDim, timesteps, n_species))
  size_A <- array(parameters$s_A_min, dim = c(xDim, yDim, zDim, timesteps, n_species))
  totalabundanceS <- numeric(timesteps)
  totalabundanceJ <- numeric(timesteps)
  totalabundanceA <- numeric(timesteps)

  sp <- run_spinup(state,
    Visualize = Visualize, n_gens = spinup,
    carCap = carCap, sleeptime = sleeptime,
    visualize_dispersion = visualize_dispersion,
    stochastic = stochastic
  )
  abundanceS[, , , 1, ] <- sp$S
  abundanceJ[, , , 1, ] <- sp$J
  abundanceA[, , , 1, ] <- sp$A
  size_S[, , , 1, ] <- sp$size_S
  size_J[, , , 1, ] <- sp$size_J
  size_A[, , , 1, ] <- sp$size_A
  totalabundanceS[1] <- sum(abundanceS[, , , 1, ])
  totalabundanceJ[1] <- sum(abundanceJ[, , , 1, ])
  totalabundanceA[1] <- sum(abundanceA[, , , 1, ])
  log_msg(sprintf(
    "Starting population: %d S, %d J, %d A",
    totalabundanceS[1], totalabundanceJ[1], totalabundanceA[1]
  ))

  fruited <- array(FALSE, dim = c(xDim, yDim, zDim, n_species))
  Disp <- sp$last_disp

  for (t in 1:(timesteps - 1)) {
    pass1 <- run_pass1_disperse(state, abundanceA, size_A, t, stochastic = stochastic)
    Disp <- pass1$Disp
    fruited <- pass1$fruited
    pass2 <- run_pass2_establish(state, abundanceS, abundanceJ, abundanceA, size_S, Disp, t, stochastic = stochastic)
    abundanceS <- pass2$S
    size_S <- pass2$size_S
    result <- run_pass3_survive_grow(
      state, abundanceS, abundanceJ, abundanceA,
      size_S, size_J, size_A, fruited, t,
      stochastic = stochastic
    )
    abundanceS <- result$S
    abundanceJ <- result$J
    abundanceA <- result$A
    size_S <- result$size_S
    size_J <- result$size_J
    size_A <- result$size_A

    # size_S needs a count-weighted blend, not a plain carry-forward: slot
    # t+1 already holds this year's newly-established seedlings (from
    # run_pass2_establish() above, size s_S_min) BEFORE the surviving/grown
    # population from slot t is added on top — see methods.tex, Population
    # state and stage structure. size_J/size_A have no such second
    # contributor here (nothing establishes directly into J or A), so they
    # still just carry forward.
    n_new <- abundanceS[, , , t + 1, ]
    n_carry <- abundanceS[, , , t, ]
    tot_n <- n_new + n_carry
    blended <- (n_new * size_S[, , , t + 1, ] + n_carry * size_S[, , , t, ]) / pmax(tot_n, 1)
    blended[tot_n == 0] <- parameters$s_S_min
    size_S[, , , t + 1, ] <- blended

    abundanceS[, , , t + 1, ] <- abundanceS[, , , t + 1, ] + abundanceS[, , , t, ]
    abundanceJ[, , , t + 1, ] <- abundanceJ[, , , t + 1, ] + abundanceJ[, , , t, ]
    abundanceA[, , , t + 1, ] <- abundanceA[, , , t + 1, ] + abundanceA[, , , t, ]
    size_J[, , , t + 1, ] <- size_J[, , , t, ] # carry size to next timestep
    size_A[, , , t + 1, ] <- size_A[, , , t, ] # carry size to next timestep
    totalabundanceS[t + 1] <- sum(abundanceS[, , , t + 1, ])
    totalabundanceJ[t + 1] <- sum(abundanceJ[, , , t + 1, ])
    totalabundanceA[t + 1] <- sum(abundanceA[, , , t + 1, ])
    if (Visualize) {
      .plot_live(
        state, abundanceS, abundanceJ, abundanceA,
        totalabundanceS, totalabundanceJ, totalabundanceA,
        t + 1, carCap, sleeptime
      )
    }
    log_msg(sprintf(
      "t=%d | S=%d J=%d A=%d total=%d", t + 1,
      totalabundanceS[t + 1], totalabundanceJ[t + 1], totalabundanceA[t + 1],
      totalabundanceS[t + 1] + totalabundanceJ[t + 1] + totalabundanceA[t + 1]
    ))
  }
  list(
    landscape = state$landscape,
    abundanceS = abundanceS, abundanceJ = abundanceJ, abundanceA = abundanceA,
    size_S = size_S, size_J = size_J, size_A = size_A,
    totalabundanceS = totalabundanceS, totalabundanceJ = totalabundanceJ,
    totalabundanceA = totalabundanceA,
    heights = state$heights, species_ids = state$species_ids,
    xDim = xDim, yDim = yDim, zDim = zDim, n_species = n_species,
    state = state, last_disp = Disp,
    obs_train = niches_train, obs_val = niches_val
  )
}

run_one <- function(params, tag = "run", timesteps = 20, spinup = 3, clim_cache = NULL, seed = 42) {
  tryCatch(
    runcolonization(
      site         = site,
      niches       = niches,
      canopy_grid  = canopy_grid,
      microenv     = microenv,
      timesteps    = timesteps,
      resolution   = 10,
      carCap       = 5,
      maxDisp      = 10,
      spinup       = spinup,
      Visualize    = FALSE,
      parameters   = params,
      forestparams = forestparams,
      clim_cache   = clim_cache,
      seed         = seed
    ),
    error = function(e) {
      # run_experiment() wraps this call in suppressMessages(), so a plain
      # message() here would vanish with no trace anywhere. Write directly
      # to the shared log file (defined by the caller, e.g.
      # run_colonization.R) so failed workers are still visible.
      err_line <- sprintf("ERROR [%s]: %s", tag, e$message)
      message(err_line)
      if (exists("log_file", inherits = TRUE)) {
        cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", err_line, "\n"),
          file = log_file, append = TRUE
        )
      }
      NULL
    }
  )
}

# ── Experiment runner ─────────────────────────────────────────────────────────
# Varies param_name across values; n_reps replicates per value.
# Returns tidy data frame of S, J, A totals over time.

run_experiment <- function(param_name, values, base = base_params,
                           n_reps = 2, timesteps = 20, spinup = 3) {
  jobs <- expand.grid(
    val = values, rep = seq_len(n_reps),
    stringsAsFactors = FALSE
  )
  cat(sprintf(
    "\n====== %s (%d jobs on %d cores) ======\n",
    param_name, nrow(jobs), N_CORES
  ))

  # Climate is identical across every value/rep in this sweep (only the
  # biological parameter differs) — build the per-height lookup table once
  # here, before mclapply forks, so all workers inherit it via copy-on-write
  # instead of every one of them re-reading every height's raster from disk.
  cat("Pre-computing shared climate lookup table for the sweep...\n")
  clim_cache <- build_clim_cache(microenv)

  rows <- mclapply(seq_len(nrow(jobs)), function(i) {
    val <- jobs$val[i]
    rep <- jobs$rep[i]
    p <- base
    p[[param_name]] <- val
    # Suppress log_msg inside workers
    suppressMessages(
      r <- run_one(p,
        tag = sprintf("%s=%s rep%d", param_name, val, rep),
        timesteps = timesteps, spinup = spinup, clim_cache = clim_cache,
        seed = i
      )
    )
    if (is.null(r)) {
      return(NULL)
    }
    T <- timesteps
    data.frame(
      param_value = as.character(val), rep = rep, t = 1:T,
      totalS = r$totalabundanceS, totalJ = r$totalabundanceJ,
      totalA = r$totalabundanceA,
      total = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
      extinct = all(r$totalabundanceA[(T %/% 2):T] == 0)
    )
  }, mc.cores = N_CORES)

  df <- do.call(rbind, Filter(Negate(is.null), rows))
  df$param_value <- factor(df$param_value, levels = as.character(values))
  cat(sprintf(
    "  Done: %d/%d runs succeeded\n",
    length(Filter(Negate(is.null), rows)), nrow(jobs)
  ))
  df
}

# ── Factorial experiment runner ────────────────────────────────────────────────
# Crosses several parameters at once (unlike run_experiment(), which varies
# only one). `param_values` is a named list, e.g.
#   list(p_poll = c(...), p_germ = c(...), p_s1 = c(...))
# giving length(p_poll) x length(p_germ) x length(p_s1) x n_reps jobs.
# Returns a tidy data frame of S/J/A totals over time, with one column per
# swept parameter recording the value used in that run.
#
# If checkpoint_var/checkpoint_path are given, jobs are processed in blocks
# — one block per unique value of checkpoint_var — and the accumulated
# results so far are saveRDS()'d to checkpoint_path after every block. This
# bounds how much work is lost if the process is killed partway through a
# long factorial (e.g. a SLURM walltime limit): at most one block's worth,
# instead of the entire sweep, since without this the only saveRDS() call
# happens after ALL jobs finish. mclapply's built-in mc.cores parallelism is
# unaffected within each block; blocks just run sequentially relative to
# each other, so total wall-clock time is essentially unchanged.
run_factorial_experiment <- function(param_values, base = base_params,
                                     n_reps = 1, timesteps = 20, spinup = 3,
                                     checkpoint_var = NULL, checkpoint_path = NULL) {
  param_names <- names(param_values)
  jobs <- do.call(
    expand.grid,
    c(param_values, list(rep = seq_len(n_reps)), stringsAsFactors = FALSE)
  )
  cat(sprintf(
    "\n====== factorial %s (%d combos x %d reps = %d jobs on %d cores) ======\n",
    paste(param_names, collapse = " x "),
    nrow(jobs) / n_reps, n_reps, nrow(jobs), N_CORES
  ))

  # Climate doesn't depend on any of the swept parameters — build it once,
  # before mclapply forks, same reasoning as run_experiment().
  cat("Pre-computing shared climate lookup table for the factorial...\n")
  clim_cache <- build_clim_cache(microenv)

  run_job <- function(i) {
    p <- base
    for (nm in param_names) p[[nm]] <- jobs[[nm]][i]
    tag <- paste(
      sprintf(
        "%s=%s", param_names,
        sapply(param_names, function(nm) jobs[[nm]][i])
      ),
      collapse = " "
    )
    tag <- paste0(tag, sprintf(" rep%d", jobs$rep[i]))
    suppressMessages(
      r <- run_one(p,
        tag = tag, timesteps = timesteps, spinup = spinup,
        clim_cache = clim_cache, seed = i
      )
    )
    if (is.null(r)) {
      return(NULL)
    }
    T <- timesteps
    df <- data.frame(
      rep = jobs$rep[i], t = 1:T,
      totalS = r$totalabundanceS, totalJ = r$totalabundanceJ,
      totalA = r$totalabundanceA,
      total = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
      extinct = all(r$totalabundanceA[(T %/% 2):T] == 0)
    )
    for (nm in param_names) df[[nm]] <- jobs[[nm]][i]
    df
  }

  use_checkpoints <- !is.null(checkpoint_var) && !is.null(checkpoint_path) &&
    checkpoint_var %in% param_names

  if (!use_checkpoints) {
    rows <- mclapply(seq_len(nrow(jobs)), run_job, mc.cores = N_CORES)
    df <- do.call(rbind, Filter(Negate(is.null), rows))
    cat(sprintf(
      "  Done: %d/%d runs succeeded\n",
      length(Filter(Negate(is.null), rows)), nrow(jobs)
    ))
    return(df)
  }

  # Resume support: if checkpoint_path already holds results from a previous
  # (e.g. walltime-killed) attempt, skip any checkpoint_var block it already
  # covers instead of redoing it. A 625-combo factorial can take longer than
  # a single SLURM walltime limit (MindoTarabita's reproduction_factorial_v3
  # run, 2026-07-16: timed out at 8h with 500/625 jobs already checkpointed,
  # 4 of 5 n_founders blocks complete) -- without this, simply resubmitting
  # would throw away every already-completed block just to redo the one that
  # didn't finish in time.
  prior_df <- NULL
  done_blocks <- character(0)
  if (file.exists(checkpoint_path)) {
    prior_df <- readRDS(checkpoint_path)
    if (!is.null(prior_df) && nrow(prior_df) > 0 && checkpoint_var %in% names(prior_df)) {
      done_blocks <- as.character(unique(prior_df[[checkpoint_var]]))
      cat(sprintf(
        "  Resuming from existing checkpoint: %d rows already done, %d/%d %s block(s) complete\n",
        nrow(prior_df), length(done_blocks),
        length(unique(jobs[[checkpoint_var]])), checkpoint_var
      ))
    }
  }

  blocks <- split(seq_len(nrow(jobs)), jobs[[checkpoint_var]])
  all_rows <- if (!is.null(prior_df)) list(prior_df) else list()
  n_done <- if (!is.null(prior_df)) nrow(prior_df) else 0L
  for (block_val in names(blocks)) {
    if (block_val %in% done_blocks) {
      cat(sprintf("  -- block %s=%s: already in checkpoint, skipping --\n", checkpoint_var, block_val))
      next
    }
    idx <- blocks[[block_val]]
    cat(sprintf("  -- block %s=%s: %d jobs --\n", checkpoint_var, block_val, length(idx)))
    rows <- Filter(Negate(is.null), mclapply(idx, run_job, mc.cores = N_CORES))
    n_done <- n_done + length(rows)
    all_rows <- c(all_rows, rows)
    df_so_far <- do.call(rbind, all_rows)
    saveRDS(df_so_far, checkpoint_path)
    cat(sprintf(
      "  Checkpoint saved (%d/%d jobs done so far): %s\n",
      n_done, nrow(jobs), checkpoint_path
    ))
  }
  cat(sprintf("  Done: %d/%d runs succeeded\n", n_done, nrow(jobs)))
  do.call(rbind, all_rows)
}

# ── Replicated runner for a single (non-swept) parameter set ───────────────────
# Runs n_reps independent replicates of one fixed params set in parallel —
# e.g. to check whether an outcome (like establishment never succeeding) is
# genuinely blocked or just one unlucky stochastic draw. Unlike
# run_experiment()/run_factorial_experiment() (which discard each run's full
# spatial arrays down to S/J/A totals, since they cover hundreds of combos),
# this keeps every replicate's complete runcolonization() output — reasonable
# since n_reps is typically a handful, not hundreds — so results stay usable
# for spatial/animation plotting (see plot_3d_abundance_animated()), not just
# aggregate totals. Returns list(runs = <one runcolonization() output per
# replicate>, summary = <tidy S/J/A-over-time data frame, like
# run_experiment()'s output but without a param_value column>).
run_replicated <- function(params, n_reps = 1, timesteps = 20, spinup = 3,
                           clim_cache = NULL) {
  cat(sprintf("\n====== %d replicate(s) on %d cores ======\n", n_reps, N_CORES))
  if (is.null(clim_cache)) clim_cache <- build_clim_cache(microenv)

  runs <- mclapply(seq_len(n_reps), function(rep) {
    suppressMessages(
      r <- run_one(params,
        tag = sprintf("rep%d", rep), timesteps = timesteps,
        spinup = spinup, clim_cache = clim_cache, seed = rep
      )
    )
    r
  }, mc.cores = min(N_CORES, n_reps))

  ok <- !vapply(runs, is.null, logical(1))
  cat(sprintf("  Done: %d/%d replicates succeeded\n", sum(ok), n_reps))

  # A run where every replicate failed (OOM-killed workers, a climate-data
  # bug, etc.) must not silently return an empty-but-valid-looking result --
  # 2026-07-27, Saloya's realistic_273founders run hit exactly this: all 5
  # replicates were OOM-killed (mclapply's own "did not deliver results"
  # warning), yet the caller still saved a "Done" .rds with a NULL summary.
  # Same class of bug as the microenv height-loop fix earlier today --
  # failing loudly here lets SLURM/the pipeline's dependency chain react
  # instead of downstream code silently consuming garbage.
  if (sum(ok) == 0) {
    stop(sprintf("run_replicated(): 0/%d replicates succeeded -- aborting without returning a result.", n_reps))
  }

  summary_df <- do.call(rbind, lapply(seq_along(runs), function(i) {
    r <- runs[[i]]
    if (is.null(r)) {
      return(NULL)
    }
    T <- timesteps
    data.frame(
      rep = i, t = 1:T,
      totalS = r$totalabundanceS, totalJ = r$totalabundanceJ,
      totalA = r$totalabundanceA,
      total = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
      extinct = all(r$totalabundanceA[(T %/% 2):T] == 0)
    )
  }))

  list(runs = runs, summary = summary_df)
}

# ── Plot: S, J, A panels ──────────────────────────────────────────────────────

plot_experiment <- function(df, param_name, title = NULL) {
  title <- title %||% sprintf("Effect of %s — Maquipucuna", param_name)

  make_panel <- function(y_var, y_lab, col) {
    mean_df <- aggregate(as.formula(paste(y_var, "~ param_value + t")),
      data = df, FUN = mean
    )
    ggplot(df, aes(
      x = t, y = .data[[y_var]], colour = param_value,
      group = interaction(param_value, rep)
    )) +
      geom_line(alpha = 0.20, linewidth = 0.4) +
      geom_line(
        data = mean_df,
        aes(
          x = t, y = .data[[y_var]], colour = param_value,
          group = param_value
        ),
        linewidth = 1.2, inherit.aes = FALSE
      ) +
      scale_colour_brewer(palette = "RdYlBu", direction = -1, name = param_name) +
      labs(x = "Year", y = y_lab) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "right")
  }

  (make_panel("totalS", "Seedlings", "#4dac26") +
    make_panel("totalJ", "Juveniles", "#f1a340") +
    make_panel("totalA", "Adults", "#08519c")) +
    plot_annotation(
      title   = title,
      caption = "Thick = mean of replicates; thin = individual runs"
    )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
