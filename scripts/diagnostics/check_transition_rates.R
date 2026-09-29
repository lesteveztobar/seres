# check_transition_rates.R
# Diagnostic: compute the actual monthly and annual-equivalent S->J and
# J->A transition probabilities under a given params list, to check whether
# a params combination causes individuals to race through the J stage
# faster than the annual census can catch them (see growth_prob() /
# run_pass3_survive_grow() in get_colonization.R -- psi0S/psi0J are
# calibrated as "monthly-equivalent" assuming DEFAULT beta_precip/beta_rh;
# pushing those slopes to extreme values, as best_case.rds does, breaks
# that calibration).
#
# Usage: Rscript scripts/diagnostics/check_transition_rates.R [params_rds]
#   Defaults to data/params/best_case.rds if no path given.
# ─────────────────────────────────────────────────────────────────────────────
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

args        <- commandArgs(trailingOnly = TRUE)
params_path <- if (length(args) >= 1 && nzchar(args[1])) args[1] else
  file.path(PARAMS_DIR, "best_case.rds")
p <- readRDS(params_path)

cat(sprintf("\n== %s ==\n", params_path))
cat(sprintf("psi0S=%.3f  psi0J=%.3f  beta_precip=%.2e  beta_rh=%.4f\n",
            p$psi0S, p$psi0J, p$beta_precip, p$beta_rh))

# Representative precip/RH range to test across -- edit if you want to check
# a specific site's actual values instead.
precip_vals <- c(1500, 2500, 3500)  # mm/yr
rh_vals     <- c(75, 85, 95)        # %

for (precip in precip_vals) {
  for (rh in rh_vals) {
    p_StoJ_month <- growth_prob("S", precip, rh, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
    p_JtoA_month <- growth_prob("J", precip, rh, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)

    # Probability of transitioning at least once across 12 months
    p_StoJ_year <- 1 - (1 - p_StoJ_month)^12
    p_JtoA_year <- 1 - (1 - p_JtoA_month)^12

    # Expected number of months spent in J before promotion to A, given a
    # constant per-month promotion probability (geometric distribution mean)
    expected_months_in_J <- if (p_JtoA_month > 0) 1 / p_JtoA_month else Inf

    cat(sprintf(
      "\nprecip=%dmm/yr, RH=%d%%:\n  p(S->J)/month=%.3f  p(J->A)/month=%.3f\n  p(S->J)/year>=1x=%.3f  p(J->A)/year>=1x=%.3f\n  expected months spent in J before promotion: %.1f\n",
      precip, rh, p_StoJ_month, p_JtoA_month, p_StoJ_year, p_JtoA_year, expected_months_in_J))
  }
}

cat("\nIf 'expected months spent in J' is well under 12, individuals promote\n",
    "out of J faster than the annual census can catch them -- J showing up\n",
    "near-zero in yearly summaries is expected under these params, not a bug.\n")
