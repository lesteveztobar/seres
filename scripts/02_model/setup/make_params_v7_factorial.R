# make_params_v7_factorial.R -- v7 reproduction-survival-seed factorial
# design (2026-09-08). Builds on reproduction_factorial_v4.rds's 4-parameter,
# 5-level grid, adding S (seeds/capsule, the v7 fecundity fix) as a 5th
# swept parameter at 3 levels (Phase 4.2's own levels: low/default/high).
# HELD -- not submitted by this script; produces the params files only, so
# the combination/run counts can be reported and approved before anything
# runs.
source("scripts/02_model/config/paths.R")

params_reprofactorial_v7 <- list(
  beta0S = -0.24 + 2.889, beta0J = 0.41 + 2.729, beta0A = 1.73 + 2.563,
  beta1 = 0.10,
  s_S_min = 0.0, s_S_max = 1.0,
  s_J_min = 1.0, s_J_max = 7.0,
  s_A_min = 7.0, s_A_max = 20.0,
  psi0S = -3.30 - 2.577, psi0J = -2.70 - 2.619,
  beta_precip = 3e-4, beta_rh = 0.010,
  sigma = 0.10, delta_s_base = 0.80,
  cost_repro = 0.50,
  p_poll = c(0.375, 0.525, 0.675, 0.825, 1.000),
  p_germ = c(0.002125, 0.004375, 0.006625, 0.008875, 0.012250),
  p_s1 = c(0.506, 0.619, 0.731, 0.844, 1.000),
  S = c(5.0e5, 1.76e6, 3.9e6),   # low / default / high (Phase 4.2 levels)
  # canopy_z omitted -- site-specific, overwritten by run_colonization.R
  lambda = 3.23, Ut = 0.23,
  n_founders = c(53, 113, 225, 400, 500)  # 2026-09-29: rescaled under the 500 cap (was 151..950)
)
saveRDS(params_reprofactorial_v7, file.path(PARAMS_DIR, "reproduction_survival_seed_v7.rds"))

n_combos <- prod(vapply(params_reprofactorial_v7, length, integer(1))[
  vapply(params_reprofactorial_v7, length, integer(1)) > 1])
cat(sprintf("Full factorial: %d combinations (p_poll x p_germ x p_s1 x S x n_founders = 5x5x5x3x5)\n", n_combos))

# ── Noise-floor add-on: 20 combinations spanning the parameter space, 5
# replicates each (not 1) -- separates a parameter effect from stochastic
# noise, which the 1-rep main grid alone cannot do (persistence outcomes
# are close to bimodal -- 0/3 replicates recruiting was already observed
# in earlier runs). Systematic span, not random: cycles each parameter
# through its own levels out of phase with the others (a resolution-V-ish
# design) so all 5 levels of every parameter appear roughly equally often
# across the 20 combos, rather than only sampling one corner of the space.
lv <- function(x, i) x[((i - 1) %% length(x)) + 1]
noise_floor_combos <- do.call(rbind, lapply(0:19, function(i) data.frame(
  p_poll      = lv(params_reprofactorial_v7$p_poll, i + 1),
  p_germ      = lv(params_reprofactorial_v7$p_germ, i + 2),
  p_s1        = lv(params_reprofactorial_v7$p_s1, i + 3),
  S           = lv(params_reprofactorial_v7$S, i + 1),
  n_founders  = lv(params_reprofactorial_v7$n_founders, i + 4)
)))
saveRDS(noise_floor_combos, file.path(PARAMS_DIR, "reproduction_survival_seed_v7_noisefloor_combos.rds"))
cat(sprintf("Noise-floor design: %d combinations x 5 replicates = %d runs (1 site)\n",
            nrow(noise_floor_combos), nrow(noise_floor_combos) * 5))
