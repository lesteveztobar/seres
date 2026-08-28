# synthetic_e2e_test.R
#
# Purely synthetic, fully in-memory end-to-end smoke test of the niche-
# scoring + establishment pipeline, written to confirm the whole chain
# (microenv-shaped object -> per-voxel quantile cache -> species niches ->
# establishment probability) produces sane, finite numbers now that
# .voxel_point_quantile()'s key-format bug is fixed (2026-08-0X).
#
# NOT a unit test harness and not meant to be run all at once with
# Rscript. It's a series of numbered, independent-ish sections separated
# by "# -- Step N: ... --" banners, each ending in a cat()/print() of its
# key result, meant to be run ONE STEP AT A TIME in an interactive R
# session (source() up to a point, or copy-paste/Ctrl+Enter section by
# section in RStudio/VS Code). Look at each step's printed output before
# moving to the next -- later steps depend on objects earlier steps leave
# in your environment (microenv, niches, state, ...).
#
# No SLURM, no real cluster data, no disk writes to Lustre scratch --
# everything lives in R objects in memory for the whole script. (The
# microenv object deliberately has NO .height_dir field -- see Step 1's
# comment -- which routes load_height()/microenv_heights() through their
# existing in-memory fallback branch instead of reading per-height .rds
# files from scratch, exactly the branch a real run never takes.)
#
# Usage:
#   source("scripts/02_model/engine/get_colonization.R")   # do this first, once
#   source("scripts/tests/synthetic_e2e_test.R")            # then either:
#     - source the whole file at once, or
#     - open it and run section by section (Ctrl+Enter / Cmd+Enter per line,
#       or select a "# -- Step N --" block and run just that selection)
#
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

if (!exists("get_niche_voxel")) {
  stop("get_colonization.R doesn't look sourced yet -- run:\n",
       '  source("scripts/02_model/engine/get_colonization.R")\n',
       "first, then source/step through this file.")
}

set.seed(20260806)

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Build synthetic microenv inputs
# ─────────────────────────────────────────────────────────────────────────────
# 3 height tiers (understory / mid-canopy / upper-canopy), a small 4x4
# raster, and a FULL SYNTHETIC YEAR of hourly data (not just a few days).
# A full year matters here specifically because get_niche_voxel() and
# voxel_background_table() are called downstream (Step 4, via
# init_colonization()) with their default months = 1:12 -- with only a
# few days of single-month data, 11 of 12 monthly lookups per observation
# would come back NA (not wrong, just a less clean demonstration). A full
# year keeps every month's quantile matrix genuinely populated, so later
# steps show real numbers, not partially-NA ones.
#
# Array size stays tiny (4*4*8760 * 5 vars * 8 bytes =~ 5.6 MB) so this is
# still fast interactively (a few seconds for Step 3's quantile reduction).

nr <- 4L
nc <- 4L
npix <- nr * nc
n_hours <- 24L * 365L

tme <- seq(as.POSIXct("2024-01-01 00:00", tz = "UTC"), by = "hour", length.out = n_hours)
hour_of_day <- as.integer(format(tme, "%H"))
month_of_t  <- as.integer(format(tme, "%m"))
is_day <- hour_of_day >= 6 & hour_of_day <= 18

heights_vec <- c(1, 8, 18)  # understory, mid-canopy, upper-canopy

# Per-pixel "shade" gradient, 0 (sunny/warm corner, row=1/col=1) to 1
# (shaded/cool corner, row=nr/col=nc) -- deliberately spatially structured,
# not just noise, so later steps can show real per-PIXEL differentiation
# (the whole point of the per-voxel refactor), not just per-height.
pix_shade <- outer(1:nr, 1:nc, function(r, c) ((r - 1) + (c - 1)) / ((nr - 1) + (nc - 1)))

# Builds one height tier's full raw (nr, nc, n_hours) arrays -- the same
# shape microclimf::runmicro()'s real output has -- plus the spatially-
# flattened *_mean fields the OLD (still-active) lookup_climate_by_height()
# pathway needs, plus the per-pixel quantile reduction the NEW pathway
# needs. Vectorized (no explicit per-hour R loop) purely for speed; the
# actual values/ranges are what matter for this test, not how they're
# generated.
build_height_entry <- function(height_m) {
  seasonal  <- 3 * sin(2 * pi * (month_of_t - 7) / 12)     # warmer mid-year
  diurnal   <- 4 * sin(2 * pi * (hour_of_day - 14) / 24)   # warm mid-afternoon
  base_temp <- 22 - 0.15 * height_m + seasonal + diurnal    # length n_hours, cooler higher up

  shade_bcast <- array(rep(as.vector(pix_shade), times = n_hours), dim = c(nr, nc, n_hours))
  temp_bcast  <- array(rep(base_temp, each = npix), dim = c(nr, nc, n_hours))
  day_bcast   <- array(rep(as.numeric(is_day), each = npix), dim = c(nr, nc, n_hours))
  sun_bcast   <- 1 - shade_bcast

  Tz <- temp_bcast - 2 * shade_bcast + array(rnorm(npix * n_hours, 0, 0.4), dim = c(nr, nc, n_hours))
  relhum_raw <- 85 - 1.5 * (Tz - (22 - 0.15 * height_m)) + array(rnorm(npix * n_hours, 0, 3), dim = c(nr, nc, n_hours))

  relhum <- pmin(pmax(relhum_raw, 30), 100)
  windspeed <- pmax(array(rnorm(npix * n_hours, 1.5, 0.6), dim = c(nr, nc, n_hours)), 0)

  Rdirdown <- pmax(day_bcast * (300 * sun_bcast +
  array(rnorm(npix * n_hours, 0, 30), dim = c(nr, nc, n_hours))), 0)

  Rdifdown <- pmax(day_bcast * (60 * sun_bcast +
  array(rnorm(npix * n_hours, 0, 10), dim = c(nr, nc, n_hours))), 0)

  h_raw <- list(Tz = Tz, relhum = relhum, windspeed = windspeed,
                Rdirdown = Rdirdown, Rdifdown = Rdifdown, tme = tme)

  list(
    # raw per-pixel arrays (kept for transparency/debugging, not required
    # downstream once voxel_quantiles below is computed)
    Tz = Tz, relhum = relhum, windspeed = windspeed,
    Rdirdown = Rdirdown, Rdifdown = Rdifdown, tme = tme,
    # spatially-flattened hourly means -- what the OLD, still-active
    # lookup_climate_by_height()/build_clim_cache() pathway reads (wind
    # attenuation coefficient `a`, mean_swdown_site, Pass 1 dispersal, the
    # annual-precipitation term in transition_logit(), etc.)
    temp_mean = apply(Tz, 3, mean), relhum_mean = apply(relhum, 3, mean),
    windspeed_mean = apply(windspeed, 3, mean),
    Rdirdown_mean = apply(Rdirdown, 3, mean), Rdifdown_mean = apply(Rdifdown, 3, mean),
    # per-pixel/month/daypart quantile reduction -- what the NEW per-voxel
    # pathway (get_clim_voxel(), and via it the niche system + Pass 2/3)
    # reads. This is the exact function run_microclimate_site.R calls at
    # ERA5-downscaling write time in the real pipeline.
    voxel_quantiles = .compute_voxel_quantiles(h_raw, nr, nc)
  )
}

microenv <- list()
for (h in heights_vec) {
  microenv[[sprintf("h%.2f", h)]] <- build_height_entry(h)
}

# Plain-vector spatial extent -- deliberately NOT an S4 terra::ext() object,
# exactly the format .spatial_extent_usable() checks for. Small ~55m x 55m
# tile split into the 4x4 raster above.
microenv$.spatial <- list(
  ext  = c(xmin = -78.7000, xmax = -78.6995, ymin = 0.0000, ymax = 0.0005),
  nrow = nr, ncol = nc,
  crs  = "synthetic, unused by any function this test exercises"
)

# Site-wide weather record (ERA5-equivalent) -- feeds precip/winddir via
# lookup_climate_by_height(); not used by the per-voxel pathway at all.
microenv$.weather <- list(
  obs_time = tme,
  precip   = pmax(0, 150 + 100 * sin(2 * pi * (month_of_t - 7) / 12) / 12 + rnorm(n_hours, 0, 20)),
  winddir  = (180 + rnorm(n_hours, 0, 30)) %% 360
)

# NOTE: microenv$.height_dir is deliberately never set. Both
# microenv_heights() and load_height() branch on
# is.null(microenv$.height_dir) -- leaving it unset routes them through
# their in-memory fallback (read microenv[["h<height>"]] directly), which
# is exactly why this whole test needs zero disk I/O.

cat("=== Step 1 result ===\n")
cat("Heights embedded:", paste(microenv_heights(microenv), collapse = ", "), "\n")
cat(".spatial_extent_usable(microenv):", .spatial_extent_usable(microenv), "\n")
cat("load_height() works with no .height_dir?",
    !is.null(load_height(microenv, heights_vec[1])), "\n")
str(microenv$.spatial)

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Build synthetic observations
# ─────────────────────────────────────────────────────────────────────────────
# 4 fake species. Sp_A and Sp_B are deliberately placed at opposite
# corners/heights (Sp_A: low + sunny/warm corner; Sp_B: high + shaded/cool
# corner) so Step 6's "near vs. far" niche-score check has an unambiguous
# expected answer. Sp_C/Sp_D sit in between, for variety.
#
# Columns match exactly what init_colonization()/get_niche_voxel() read
# (confirmed by grepping every site_obs$/obs_sp$/niches$ reference in
# get_colonization.R): Area_or_Site, FinalID, Height_m, lat, lon.

site_name <- "SyntheticSite"

niches <- data.frame(
  Area_or_Site = site_name,
  FinalID = c(rep("Sp_A", 3), rep("Sp_B", 3), rep("Sp_C", 2), rep("Sp_D", 2)),
  Height_m = c(0.8, 1.0, 1.2,      # Sp_A: understory
               17.0, 18.0, 19.0,   # Sp_B: upper canopy
               8.0, 8.5,           # Sp_C: mid-canopy
               1.0, 8.0),          # Sp_D: broad (understory + mid)
  lon = c(-78.69990, -78.69985, -78.69988,   # Sp_A: near xmin corner (sunny)
          -78.69955, -78.69952, -78.69958,   # Sp_B: near xmax corner (shaded)
          -78.69970, -78.69972,               # Sp_C: middle
          -78.69965, -78.69975),              # Sp_D: middle-ish
  lat = c(0.00046, 0.00048, 0.00044,          # Sp_A: near ymax (sunny corner, row~1)
          0.00006, 0.00004, 0.00008,          # Sp_B: near ymin (shaded corner, row~nrow)
          0.00025, 0.00027,
          0.00030, 0.00020),
  stringsAsFactors = FALSE
)

cat("=== Step 2 result ===\n")
print(niches)
cat("\nn species:", length(unique(niches$FinalID)), "| n observations:", nrow(niches), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Run the quantile/cache layer
# ─────────────────────────────────────────────────────────────────────────────
# .compute_voxel_quantiles() already ran once per height inside Step 1
# (build_height_entry()); this step re-demonstrates it standalone on one
# height's raw arrays for inspection, then runs the caching layer
# (build_clim_cache_voxel()) across the whole synthetic microenv the way
# init_colonization() will in Step 4 (footprint = NULL here -> "pooled"/
# unsubset mode, i.e. every pixel kept, nothing landscape-restricted yet).

h1_raw <- list(
  Tz = microenv[["h1.00"]]$Tz, relhum = microenv[["h1.00"]]$relhum,
  windspeed = microenv[["h1.00"]]$windspeed, Rdirdown = microenv[["h1.00"]]$Rdirdown,
  Rdifdown = microenv[["h1.00"]]$Rdifdown, tme = microenv[["h1.00"]]$tme
)
vq_demo <- .compute_voxel_quantiles(h1_raw, nr, nc)

cat("=== Step 3a: .compute_voxel_quantiles() on height=1m, standalone ===\n")
cat("mode:", vq_demo$mode, "| nr x nc:", vq_demo$nr, "x", vq_demo$nc, "\n")
cat("n quantile keys:", length(vq_demo$quantiles), "\n")
cat("sample keys:\n")
print(head(names(vq_demo$quantiles), 8))
cat("any 'pooled_' key written? (should be FALSE -- confirms the fix's diagnosis)\n")
print(any(grepl("^pooled_", names(vq_demo$quantiles))))
cat("one real quantile row, height=1m, month=6, day, temp, pixel (1,1):\n")
print(vq_demo$quantiles[["6_day_temp"]][1, ])  # 5 numbers: the QUANTILE_PROBS percentiles

clim_cache_voxel <- build_clim_cache_voxel(microenv, footprint = NULL)

cat("\n=== Step 3b: build_clim_cache_voxel() across the whole synthetic microenv ===\n")
cat("heights in cache:", paste(clim_cache_voxel$heights, collapse = ", "), "\n")
str(clim_cache_voxel$clim_voxel_by_height[["1"]], max.level = 1)

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Run init_colonization()
# ─────────────────────────────────────────────────────────────────────────────
# Minimal forestparams (literature-scale values from methods.tex's own
# Table~\ref{tab:params_forest}, mean_hgt/sd_hgt bumped up so trees can
# plausibly reach the 18m upper-canopy tier) and a near-empty params list
# (init_colonization() fills in its own defaults for anything else it
# needs -- s_A_min/max, delta_s_base, cost_repro, a, mean_swdown_site).
# canopy_grid is a required argument but unused whenever forestparams is
# supplied (see build_forest() branch in init_colonization()) -- passed
# as a small dummy matrix.

site <- list(Site = site_name)
canopy_grid_dummy <- matrix(20, nrow = 5, ncol = 5)  # unused; forestparams takes priority

forestparams <- list(
  stems_per_ha = 298,
  mean_hgt = 15, sd_hgt = 6,          # taller/more spread than production (8.4m) so
  mean_crown_r = 2.0, sd_crown_r = 0.5, # some stochastic trees plausibly reach 18m
  trunk_r = 0.114,
  branch_density = 3.0,
  epiphyte_footprint_m2 = 0.02
)

params <- list(
  n_founders = 10,
  s_S_min = 0.0  # NOT auto-defaulted by init_colonization() (only s_A_min/max are) --
                 # run_pass2_establish() writes newly-established seedlings to this size
)

state <- init_colonization(
  site = site, niches = niches, canopy_grid = canopy_grid_dummy, microenv = microenv,
  resolution = 10, carCap = 1, maxDisp = 5, params = params,
  forestparams = forestparams, allsites = FALSE
)

cat("=== Step 4 result ===\n")
cat("init_colonization() completed:", !is.null(state), "\n")
cat("Landscape dims (x,y,z):", state$xDim, state$yDim, state$zDim, "\n")
cat("n species:", state$n_species, "| species_ids:", paste(state$species_ids, collapse = ", "), "\n")
cat("Per-voxel spatial resolution active (clim_pixel_row non-NULL)?",
    !is.null(state$clim_pixel_row), "\n")
cat("Valid canopy voxels per height tier:\n")
for (zi in seq_len(state$zDim)) {
  cat(sprintf("  height %.1fm (zi=%d): %d valid voxels\n",
              state$heights[zi], zi, sum(state$landscape[, , zi])))
}

sp1 <- state$species_ids[1]
cat("\n--- niches_by_species[[\"", sp1, "\"]] ---\n", sep = "")
str(state$niches_by_species[[sp1]], max.level = 2)
cat("\nceiling for", sp1, ":", state$niches_by_species[[sp1]]$ceiling, "\n")
cat("Is it a real finite number (not NA/NULL)?",
    is.finite(state$niches_by_species[[sp1]]$ceiling), "\n")

cat("\nceilings for all species (should all be finite numbers, most near 100):\n")
for (sp in state$species_ids) {
  ceil <- state$niches_by_species[[sp]]$ceiling
  cat(sprintf("  %-6s ceiling = %s\n", sp, if (is.null(ceil)) "NULL" else round(ceil, 2)))
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: Run run_pass2_establish()
# ─────────────────────────────────────────────────────────────────────────────
# One timestep. Abundance/size arrays sized [xDim, yDim, zDim, 2, n_species]
# (t=1 -> tnext=2), all starting empty. Seeds are injected directly into
# state$dispersalmatrix (normally Pass 1's output) at a handful of REAL
# valid-canopy voxels, found by inspecting state$landscape directly rather
# than guessing coordinates blind.

n_timesteps <- 2L
dims_abund <- c(state$xDim, state$yDim, state$zDim, n_timesteps, state$n_species)
abundanceS <- array(0L, dim = dims_abund)
abundanceJ <- array(0L, dim = dims_abund)
abundanceA <- array(0L, dim = dims_abund)
size_S     <- array(NA_real_, dim = dims_abund)

pad <- state$maxDisp
dispersalmatrix <- state$dispersalmatrix  # all zeros, shape (xDim+2pad, yDim+2pad, zDim+2pad, n_species)

# Seed a handful of real valid-canopy voxels per height tier, for the first
# two species, with a modest seed count each.
set.seed(1)
n_seeded <- 0L
for (zi in seq_len(state$zDim)) {
  valid_xy <- which(state$landscape[, , zi], arr.ind = TRUE)
  if (nrow(valid_xy) == 0) {
    cat(sprintf("height tier %d (%.1fm): no valid canopy voxels in this stochastic forest -- skipping seeding here.\n",
                zi, state$heights[zi]))
    next
  }
  pick <- valid_xy[sample.int(nrow(valid_xy), min(3, nrow(valid_xy))), , drop = FALSE]
  for (sp in seq_len(min(2, state$n_species))) {
    for (row in seq_len(nrow(pick))) {
      x <- pick[row, 1]; y <- pick[row, 2]
      dispersalmatrix[x + pad, y + pad, zi + pad, sp] <- sample(5:20, 1)
      n_seeded <- n_seeded + 1L
    }
  }
}
cat("=== Step 5a: seeding ===\n")
cat("Seeded", n_seeded, "(voxel, species) seed-count entries across", state$zDim, "height tiers.\n")
cat("Total seeds injected:", sum(dispersalmatrix), "\n")

pass2_result <- run_pass2_establish(
  state = state, abundanceS = abundanceS, abundanceJ = abundanceJ, abundanceA = abundanceA,
  size_S = size_S, dispersalmatrix = dispersalmatrix, t = 1L, stochastic = FALSE
)

cat("\n=== Step 5b: run_pass2_establish() result ===\n")
newly_established <- pass2_result$S[, , , 2, , drop = FALSE]
cat("Total newly-established seedlings (tnext slot, summed):", sum(newly_established), "\n")
cat("Any NA in the result?", anyNA(pass2_result$S), "\n")
cat("Range of established counts (non-zero voxels only):\n")
nz <- newly_established[newly_established > 0]
if (length(nz) > 0) print(summary(nz)) else cat("  (all zero -- see the seeding/forest-sparsity note below)\n")

cat("\nPer-height-tier established seedling totals:\n")
for (zi in seq_len(state$zDim)) {
  cat(sprintf("  height %.1fm (zi=%d): %d established\n",
              state$heights[zi], zi, sum(pass2_result$S[, , zi, 2, ])))
}
cat("\n(If everything above is 0: this small a stochastic forest can genuinely\n",
    " place zero canopy at a given height tier by chance -- check Step 4's\n",
    " 'Valid canopy voxels per height tier' printout first. Re-running Step 1\n",
    " onward with a different set.seed(), or raising forestparams$stems_per_ha,\n",
    " should fix it if that's the cause.)\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# Step 6: Spot-check niche_match_array()
# ─────────────────────────────────────────────────────────────────────────────
# Uses Sp_A's and Sp_B's niches directly (from Step 4). Sp_A was observed
# warm/sunny/low (Step 2); Sp_B was observed cool/shaded/high. A climate
# point resembling Sp_A's own observations should score much higher under
# Sp_A's niche than a point resembling Sp_B's -- and vice versa.

niche_A <- state$niches_by_species[["Sp_A"]]
niche_B <- state$niches_by_species[["Sp_B"]]

# A handful of explicit (temp, relhum, swdown) points to score, using the
# scalar path (niche_match()) for easy one-at-a-time reading:
points <- list(
  warm_dry_sunny = c(temp = 24, relhum = 65, swdown = 260),   # resembles Sp_A's conditions
  cool_humid_shaded = c(temp = 17, relhum = 92, swdown = 40),  # resembles Sp_B's conditions
  middling = c(temp = 20, relhum = 80, swdown = 140)
)

cat("=== Step 6a: scalar niche_match() spot-checks ===\n")
for (pname in names(points)) {
  pt <- as.list(points[[pname]])
  score_A <- niche_match(pt, niche_A)
  score_B <- niche_match(pt, niche_B)
  cat(sprintf("  point '%s' (T=%.0f RH=%.0f swdown=%.0f): Sp_A match=%.3f | Sp_B match=%.3f\n",
              pname, pt$temp, pt$relhum, pt$swdown, score_A, score_B))
}
cat("\nExpectation: 'warm_dry_sunny' should score higher for Sp_A than Sp_B;\n",
    "'cool_humid_shaded' should score higher for Sp_B than Sp_A.\n", sep = "")

# Array form (niche_match_array()) on the same three points at once, to
# confirm the vectorized path used by run_pass2_establish() agrees with
# the scalar spot-checks above.
temp_arr <- array(sapply(points, function(p) p["temp"]), dim = c(1, length(points)))
relhum_arr <- array(sapply(points, function(p) p["relhum"]), dim = c(1, length(points)))
swdown_arr <- array(sapply(points, function(p) p["swdown"]), dim = c(1, length(points)))

match_A_arr <- niche_match_array(list(temp = temp_arr, relhum = relhum_arr, swdown = swdown_arr), niche_A)
match_B_arr <- niche_match_array(list(temp = temp_arr, relhum = relhum_arr, swdown = swdown_arr), niche_B)

cat("\n=== Step 6b: niche_match_array(), same 3 points, vectorized ===\n")
cat("Sp_A matches:", paste(round(as.vector(match_A_arr), 3), collapse = ", "),
    " (order:", paste(names(points), collapse = ", "), ")\n")
cat("Sp_B matches:", paste(round(as.vector(match_B_arr), 3), collapse = ", "),
    " (order:", paste(names(points), collapse = ", "), ")\n")
cat("Scalar and array paths agree?",
    isTRUE(all.equal(as.vector(match_A_arr), sapply(points, function(p) niche_match(as.list(p), niche_A)))) &&
    isTRUE(all.equal(as.vector(match_B_arr), sapply(points, function(p) niche_match(as.list(p), niche_B)))),
    "\n")

# The agreement check above compares as.vector(match_A_arr) (plain, unnamed
# -- as.vector() strips all attributes off an array) against
# sapply(points, ...) (named by points' own names) -- if all.equal() reports
# unequal here despite the printed values matching, that's almost certainly*
# an attribute (names) mismatch, not a value mismatch. This block confirms
# which one it actually is before trusting/fixing the check above.
scalar_A <- sapply(points, function(p) niche_match(as.list(p), niche_A))
cat("=== Structural diagnostic ===\n")
cat("str(match_A_arr):\n"); str(match_A_arr)
cat("\nstr(scalar_A):\n"); str(scalar_A)
cat("\ndim(match_A_arr):", dim(match_A_arr), "\n")
cat("dim(scalar_A):", if (is.null(dim(scalar_A))) "NULL (plain vector)" else dim(scalar_A), "\n")
cat("class(match_A_arr):", class(match_A_arr), "\n")
cat("class(scalar_A):", class(scalar_A), "\n")
cat("names(match_A_arr):", names(match_A_arr), "\n")
cat("names(scalar_A):", names(scalar_A), "\n")
cat("as.vector(match_A_arr) == unname(scalar_A):", as.vector(match_A_arr) == unname(scalar_A), "\n")

cat("All scores in [0,1] (niche_match) / would be [0,100] before the /100 in niche_overall_score()?",
    all(match_A_arr >= 0 & match_A_arr <= 1, na.rm = TRUE) &&
    all(match_B_arr >= 0 & match_B_arr <= 1, na.rm = TRUE), "\n")

cat("\n=== End of synthetic_e2e_test.R ===\n")
