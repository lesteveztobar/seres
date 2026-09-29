# make_b4_params.R -- one-off: write the literature-default params (as
# run_colonization.R itself would construct them) plus n_reps=3, for the
# Phase B4 equivalence test. canopy_z/n_founders are overwritten by
# run_colonization.R per-site regardless, so their placeholder values here
# don't matter.
params <- list(
  beta0S = -0.24 + 2.889, beta0J = 0.41 + 2.729, beta0A = 1.73 + 2.563,
  beta1 = 0.10, s_S_min = 0.0, s_S_max = 1.0, s_J_min = 1.0, s_J_max = 7.0,
  s_A_min = 7.0, s_A_max = 20.0,
  psi0S = -3.30 - 2.577, psi0J = -2.70 - 2.619,
  beta_precip = 3e-4, beta_rh = 0.010, sigma = 0.10, delta_s_base = 0.80,
  cost_repro = 0.50,
  S = 1.76e6, p_poll = 0.30, p_germ = 0.001, p_s1 = 0.45,
  canopy_z = 15, lambda = 1, Ut = 1,
  n_founders = 30, n_reps = 3
)
saveRDS(params, "data/params/b4_equivalence_params.rds")
cat("Saved data/params/b4_equivalence_params.rds\n")
