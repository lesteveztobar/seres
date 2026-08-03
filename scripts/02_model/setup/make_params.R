source("scripts/02_model/config/paths.R")

# beta0S/beta0J/beta0A and psi0S/psi0J below are all written as the original
# literature-calibrated value plus an explicit offset (e.g. `-0.24 + 2.889`),
# rather than a single pre-computed number, so the monthly-compounding
# recalibration stays visible and auditable here. survival_logit() is
# applied once per month and compounded over 12 months
# (run_pass3_survive_grow(), get_colonization.R) — the offsets convert each
# intercept from "single evaluation = annual target" (the literature
# calibration) to "single evaluation = monthly-equivalent rate", so 12
# compounded months reproduce the annual target on average instead of
# crashing it by several orders of magnitude. Survival offsets solve
# p_month = p_annual^(1/12) (survival must succeed every month, so it
# compounds multiplicatively); transition offsets solve
# q_month = 1-(1-p_annual)^(1/12) (a transition is an "at least once this
# year" event, so the complement compounds). See methods.tex for the
# full derivation.

# Founders per species (see get_colonization.R::run_spinup()), used as a
# fixed baseline across every experiment below except the n_founders sweep
# itself. Update this once that sweep identifies a value that actually
# persists, then regenerate all the other params files so they pick it up.
N_FOUNDERS_DEFAULT <- 30

params_ppoll <- list(
    # ── s(z, e): survival ────────────────────────────────────────────────────
    beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
    beta1   =  0.10,
    s_S_min =  0.0,  s_S_max =  1.0,
    s_J_min =  1.0,  s_J_max =  7.0,
    s_A_min =  7.0,  s_A_max = 20.0,
    # ── g(z'|z, e): growth / stage transitions ────────────────────────────
    psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_s_base =  0.80,
    cost_repro   =  0.50,
    # ── p_r(s) × f_s(s): fecundity ──────────────────────────────────────
    p_poll  = c(0.05, 0.15, 0.30, 0.50, 0.70),  p_germ  = 0.001,  p_s1 = 0.45,
    # ── d(x'|x): dispersal ────────────────────────────────────────────────
    # canopy_z omitted — it's site-specific and run_colonization.R
    # always overwrites it with the per-site mean canopy height.
    lambda = 1,  Ut = 1,
    # ── spin-up ────────────────────────────────────────────────────────────
    n_founders = N_FOUNDERS_DEFAULT
)

out_path <- file.path(PARAMS_DIR, "p_poll.rds")
saveRDS(params_ppoll, out_path)


params_beta0A <- list(
    # ── s(z, e): survival ────────────────────────────────────────────────────
    beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,
    beta0A  =  c(0.50, 1.00, 1.73, 2.50, 3.50) + 2.563,  # shifted range preserves relative meaning of each level
    beta1   =  0.10,
    s_S_min =  0.0,  s_S_max =  1.0,
    s_J_min =  1.0,  s_J_max =  7.0,
    s_A_min =  7.0,  s_A_max = 20.0,
    # ── g(z'|z, e): growth / stage transitions ────────────────────────────
    psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_s_base =  0.80,
    cost_repro   =  0.50,
    # ── p_r(s) × f_s(s): fecundity ──────────────────────────────────────
    p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
    # ── d(x'|x): dispersal ──────────────────────────────────────────────
    # canopy_z omitted — site-specific, overwritten by run_colonization.R
    lambda = 1,  Ut = 1,
    n_founders = N_FOUNDERS_DEFAULT
)


out_path <- file.path(PARAMS_DIR, "beta0A.rds")
saveRDS(params_beta0A, out_path)

params_pgerm <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
# ── p_r(s) × f_s(s): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = c(0.0001, 0.0005, 0.001, 0.003, 0.005),  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
)

out_path <- file.path(PARAMS_DIR, "p_germ.rds")
saveRDS(params_pgerm, out_path)

params_costrepro <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  c(0.20, 0.40, 0.60, 0.80, 1.00),
# ── p_r(s) × f_s(s): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
)

out_path <- file.path(PARAMS_DIR, "cost_repro.rds")
saveRDS(params_costrepro, out_path)

params_betarh <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  c(0.002, 0.005, 0.010, 0.020, 0.040),
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
# ── p_r(s) × f_s(s): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
)

out_path <- file.path(PARAMS_DIR, "beta_rh.rds")
saveRDS(params_betarh, out_path)

params_betaprecip <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  c(1e-4, 2e-4, 3e-4, 6e-4, 1e-3),  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
# ── p_r(s) × f_s(s): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
)

out_path <- file.path(PARAMS_DIR, "beta_precip.rds")
saveRDS(params_betaprecip, out_path)

# ── Spin-up: founder count ──────────────────────────────────────────────────────
# Sweeps n_founders per species (see get_colonization.R::run_spinup()) —
# decoupled from field observation count. Run this experiment alone first:
# the fecundity math (reproduce(), get_colonization.R:226-234) means expected
# seed output per adult per year is tiny even at the most generous swept
# p_poll/p_germ/p_s1 combination, so persistence may hinge on founder count
# more than on any vital rate. Whatever value here actually persists should
# become the new N_FOUNDERS_DEFAULT above, then regenerate the other params
# files so every other experiment (including the reproduction factorial)
# uses a founder count that's actually viable.
params_nfounders <- list(
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = seq(5, 100, by = 5)
)

out_path <- file.path(PARAMS_DIR, "n_founders.rds")
saveRDS(params_nfounders, out_path)

# ── Reproduction factorial: n_founders x p_poll x p_germ x p_s1 ────────────────
# Full factorial across founder count and the three components of fecundity
# (pollination success, mycorrhizal germination, first-year seed-to-seedling
# survival) — 5^4 = 625 combinations, run via run_factorial_experiment().
# n_founders alone (5-100, at literature-default reproduction) came back
# 100% extinct at every level tested — the fecundity formula makes expected
# seed output per adult per year ~0.00006 at baseline, so no realistic
# founder count rescues it there. Crossing n_founders into the factorial
# instead of testing it in isolation finds the actual joint threshold
# directly, mirroring how the simple model's Fig. 5 (report.pdf sec. 3.3)
# found the p_est x lambda threshold in one factorial rather than
# sequentially. n_founders levels span from what we already know fails
# (50, 100) up to what the math says should be comfortably enough at the
# most generous reproduction corner tested (400, 800): at
# p_poll=0.7/p_germ=0.005/p_s1=0.75, expected seeds/adult/yr ~0.0011, so
# n_founders=200 gives ~5 expected seeds over 30 years, n_founders=800 ~30.
#
# SUPERSEDED by params_reprofactorial_v2 below, after the 2026-07-11
# run_pass2_establish() fix (removed a double-application of p_germ that
# made establishment ~1000x more restrictive than intended -- see
# report/methods.tex, Pass 2: Establishment). The result this design
# produced (archived to data/processed/pre_pgerm_fix/) is not representative
# of current model behavior. Left here for provenance/reproducibility only
# -- don't regenerate reproduction_factorial.rds from this block expecting
# it to mean the same thing anymore.
params_reprofactorial <- list(
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
p_poll  = c(0.05, 0.15, 0.30, 0.50, 0.70),
p_germ  = c(0.0001, 0.0005, 0.001, 0.003, 0.005),
p_s1    = c(0.15, 0.30, 0.45, 0.60, 0.75),
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = c(50, 100, 200, 400, 800)
)

out_path <- file.path(PARAMS_DIR, "reproduction_factorial.rds")
saveRDS(params_reprofactorial, out_path)

# ── Reproduction factorial v2: re-bracketed post-p_germ-fix ────────────────────
# The v1 ranges above were calibrated to rescue persistence under the old,
# ~1000x more restrictive establishment gate. Post-fix, best_case.rds (every
# lever pushed to its most generous tested value, including p_poll=0.90,
# p_germ=0.01, p_s1=0.90, n_founders=1000) persists robustly (5/5 replicates,
# ~2200 adults by year 30 at Maquipucuna h0.25 -- see
# data/processed/colonization_Maquipucuna_best_case_h0.25.rds). So the open
# question isn't "can the model persist at all" anymore, it's "how far below
# that extreme corner does it still persist" -- v1's ranges don't actually
# bracket that: e.g. p_poll topped out at 0.70 and p_s1 at 0.75, both below
# what best_case needed, so a v1-shaped sweep could miss the transition
# entirely. Every level below re-brackets from literature-default up through
# best_case's own value for that parameter, holding everything else (this
# is still ONLY a reproduction/founder-count factorial -- survival, growth,
# and dispersal stay at literature defaults, same as v1) unchanged so the
# result isolates the recruitment threshold specifically:
#   n_founders: 10, 30 (default), 100, 300, 1000 (best_case)
#   p_poll:     0.15, 0.30 (default), 0.50, 0.70, 0.90 (best_case)
#   p_germ:     0.0001, 0.0005, 0.001 (default), 0.005, 0.01 (best_case)
#   p_s1:       0.15, 0.30, 0.45 (default), 0.60, 0.90 (best_case)
# Still 5^4 = 625 combinations, same cost as v1.
#
# SUPERSEDED by params_reprofactorial_v3 below, after realistic.rds and
# best_case.rds were both actually run to completion (2026-07-14):
# realistic.rds (literature defaults throughout) collapses toward extinction
# (total abundance 59.8 -> 17.8 -> 5.6 over 30 years, all-adult, no
# recruitment behind it -- see data/processed/colonization_Maquipucuna_
# realistic_h0.25.rds), while best_case.rds sustains ~2200 adults through
# year 30. v2's levels still included sub-realistic values (e.g. p_poll=0.15
# is below realistic's 0.30) that are now known to be uninformative -- if
# realistic itself dies, anything below it dies too. Left here for
# provenance only.
params_reprofactorial_v2 <- list(
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
p_poll  = c(0.15, 0.30, 0.50, 0.70, 0.90),
p_germ  = c(0.0001, 0.0005, 0.001, 0.005, 0.01),
p_s1    = c(0.15, 0.30, 0.45, 0.60, 0.90),
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = c(10, 30, 100, 300, 1000)
)

out_path <- file.path(PARAMS_DIR, "reproduction_factorial_v2.rds")
saveRDS(params_reprofactorial_v2, out_path)

# ── Reproduction factorial v3: bracketed realistic -> best_case ────────────────
# realistic.rds and best_case.rds are now the two confirmed endpoints of the
# question this factorial is answering: literature-default reproduction
# collapses (realistic), maximally-favourable reproduction sustains ~2200
# adults (best_case). v3 drops every sub-realistic level from v2 and instead
# takes 5 evenly-spaced (linear) levels from each parameter's realistic value
# up to its best_case value, inclusive, for all four factorial parameters --
# still fully crossed, still 5^4 = 625 combinations, same cost as v1/v2, just
# spending every combination inside the range that's actually informative.
# Everything else (survival/growth/dispersal) stays at literature defaults,
# same as v1/v2, so the result still isolates the recruitment threshold:
#   n_founders: 30, 273, 515, 758, 1000
#   p_poll:     0.30, 0.45, 0.60, 0.75, 0.90
#   p_germ:     0.00100, 0.00325, 0.00550, 0.00775, 0.01000
#   p_s1:       0.450, 0.563, 0.675, 0.788, 0.900
# Levels are linearly spaced (not log-spaced) for simplicity, per the
# starting-coarse approach here -- if the persistence transition turns out to
# sit close to the realistic end of a given parameter's range, a follow-up
# sweep can re-space that parameter logarithmically or narrow the range
# further around wherever the transition actually falls.
params_reprofactorial_v3 <- list(
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
p_poll  = c(0.30, 0.45, 0.60, 0.75, 0.90),
p_germ  = c(0.00100, 0.00325, 0.00550, 0.00775, 0.01000),
p_s1    = c(0.450, 0.563, 0.675, 0.788, 0.900),
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = c(30, 273, 515, 758, 1000)
)

out_path <- file.path(PARAMS_DIR, "reproduction_factorial_v3.rds")
saveRDS(params_reprofactorial_v3, out_path)

# ── Best case: everything pushed as favourable as possible ─────────────────────
# Every lever more generous than anything tested so far in the OAT sweeps or
# the factorial — not meant to be realistic, just to answer "can the model
# persist at all?" n_reps replicates (see run_replicated() in
# get_colonization.R) rather than a single run, since a single zero-
# establishment outcome can't distinguish "genuinely blocked" from "just an
# unlucky stochastic draw" — see project memory on the niche/founder/
# capacity-saturation investigation for why this matters here specifically.
params_bestcase <- list(
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  3.50 + 2.563,   # max tested (best survival)
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  1e-3,  beta_rh      =  0.040,          # max tested (best growth)
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.20,                                   # min tested (least reproduction cost)
p_poll  = 0.90,  p_germ  = 0.01,  p_s1 = 0.90,           # beyond anything tested so far
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = 1000,                                       # beyond the 800 max tested
n_reps     = 5
)

out_path <- file.path(PARAMS_DIR, "best_case.rds")
saveRDS(params_bestcase, out_path)

# ── Realistic: literature-default values, nothing pushed to an extreme ─────────
# The complement to best_case.rds: same structure (n_reps replicates via
# run_replicated(), for the same stochastic-bad-luck reasons noted above),
# but every value held at its literature-calibrated default -- i.e. the exact
# fallback params run_colonization.R uses when no params_file is
# given. best_case answers "can the model persist under the most generous
# corner of the tested space?"; realistic answers "does it persist under the
# values the model is actually calibrated to?" Comparing the two tells us
# whether the factorial sweep (params_reprofactorial_v2 above) needs its
# lower bound raised off the literature default, or whether the default
# already sits inside the persisting region.
params_realistic <- list(
beta0S  = -0.24 + 2.889,  beta0J  =  0.41 + 2.729,  beta0A  =  1.73 + 2.563,
beta1   =  0.10,
s_S_min =  0.0,  s_S_max =  1.0,
s_J_min =  1.0,  s_J_max =  7.0,
s_A_min =  7.0,  s_A_max = 20.0,
psi0S        = -3.30 - 2.577,  psi0J        = -2.70 - 2.619,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_s_base =  0.80,
cost_repro   =  0.50,
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# canopy_z omitted — site-specific, overwritten by run_colonization.R
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT,
n_reps     = 5
)

out_path <- file.path(PARAMS_DIR, "realistic.rds")
saveRDS(params_realistic, out_path)

# ── Realistic, more founders: same as realistic.rds but n_founders bumped ──────
# realistic.rds's final populations are tiny (2026-07-17 height_resolution_
# experiment.R run: final total abundance means of 4-7 individuals at t=30,
# sd nearly half the mean) -- too small/noisy to tell a genuine
# height-resolution effect apart from demographic-stochasticity noise (that
# run's Kruskal-Wallis was not significant, p=0.14, but with populations this
# small that's as likely to reflect low power as a true null result). Bumps
# n_founders from realistic's literature-default 30 up to 273 -- the second
# level from the reproduction_factorial_v3 sweep (30->273->515->758->1000),
# reusing an already-tested value rather than an arbitrary new one -- to get
# a bigger, less noise-dominated population while keeping every other
# parameter at its literature-calibrated (realistic) value.
params_realistic_273founders <- modifyList(params_realistic, list(n_founders = 273))

out_path <- file.path(PARAMS_DIR, "realistic_273founders.rds")
saveRDS(params_realistic_273founders, out_path)
