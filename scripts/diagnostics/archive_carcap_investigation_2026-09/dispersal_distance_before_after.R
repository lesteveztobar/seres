# dispersal_distance_before_after.R
# Phase 2.2 sanity-gate addition: mean dispersal distance per site, BEFORE
# (canopy_z = MEAN measured CanopyHeight_m, the pre-v7 definition) vs AFTER
# (canopy_z = MAX, Phase 1.6) the ceiling-reconciliation fix. Calls the real
# init_colonization()/disperse() machinery (not a hand-derived formula) so
# the actual `a` (canopy-openness coefficient, `4 * difrac` -- itself
# indirectly sensitive to the ceiling choice via `mid_zi`, see below) is
# used, not assumed constant across the two scenarios.
#
# CORRECTION to Phase 0's canopy_grid claim: init_colonization() (L1411,
# get_colonization.R) sets `mid_zi` (the height at which the diffuse-
# radiation fraction `difrac`/`a` is evaluated) from
# `mean(canopy_grid, na.rm=TRUE) / 2` UNCONDITIONALLY -- not gated behind
# `forestparams=NULL` the way the flat-canopy-grid landscape fallback
# (L1353-1360) is. Phase 0 called canopy_grid "dead in every current
# production run" -- true for the flat-landscape fallback specifically, but
# NOT true for this one line: canopy_grid's value (now MAX-based,
# `canopy_ceiling`) still shapes where `a` gets evaluated in every run.
# Corrected here rather than left standing.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
source("scripts/02_model/lib_logging.R")
log_msg <- make_log_msg()

sites <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")
niches_all <- load_observations()
niches_all <- niches_all[!is.na(niches_all$lat) & !is.na(niches_all$lon) &
                          !is.na(niches_all$Height_m) & !is.na(niches_all$FinalID), ]

params_base <- list(
  beta0S = -0.24 + 2.889, beta0J = 0.41 + 2.729, beta0A = 1.73 + 2.563, beta1 = 0.10,
  s_S_min = 0.0, s_S_max = 1.0, s_J_min = 1.0, s_J_max = 7.0, s_A_min = 7.0, s_A_max = 20.0,
  psi0S = -3.30 - 2.577, psi0J = -2.70 - 2.619, beta_precip = 3e-4, beta_rh = 0.010,
  sigma = 0.10, delta_s_base = 0.80, cost_repro = 0.50,
  S = 1.76e6, p_poll = 0.30, p_germ = 0.001, p_s1 = 0.45,
  lambda = 1, Ut = 1, n_founders = 30
)

for (s in sites) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", s))
  if (!file.exists(microenv_path)) { cat(s, ": no microenv\n"); next }
  microenv <- readRDS(microenv_path)
  heights <- microenv_heights(microenv)
  niches <- niches_all

  ceiling_old <- suppressWarnings(mean(niches$CanopyHeight_m[niches$Area_or_Site == s], na.rm = TRUE))
  if (!is.finite(ceiling_old)) ceiling_old <- max(heights)
  ceiling_new <- site_canopy_ceiling(s, niches, max(heights))

  rep_height <- mean(niches$Height_m[niches$Area_or_Site == s], na.rm = TRUE)
  fp <- site_forestparams(s, ceiling_new)

  run_one <- function(ceiling) {
    canopy_grid <- matrix(ceiling, nrow = 50, ncol = 50)
    params <- params_base
    params$canopy_z <- ceiling
    state <- tryCatch(
      init_colonization(list(Site = s), niches, canopy_grid, microenv,
                        resolution = 10, carCap = 1, maxDisp = 5, params = params,
                        forestparams = fp),
      error = function(e) { message(s, ": init_colonization failed -- ", conditionMessage(e)); NULL }
    )
    if (is.null(state)) return(NULL)
    a <- state$params$a
    zi <- which.min(abs(state$heights - rep_height))
    wind <- mean(state$clim_by_height[[zi]]$windspeed, na.rm = TRUE)
    meanDisp <- max(1, min(round((wind * exp(a * (rep_height - ceiling) / ceiling)) / (params$lambda * params$Ut)), 5))
    list(a = a, meanDisp = meanDisp)
  }

  old <- run_one(ceiling_old)
  new <- run_one(ceiling_new)
  if (is.null(old) || is.null(new)) next

  cat(sprintf(
    "%-15s rep_height=%.2fm | OLD ceiling(mean)=%.1fm a=%.3f meanDisp=%dm | NEW ceiling(max)=%.1fm a=%.3f meanDisp=%dm | ratio=%.2f\n",
    s, rep_height, ceiling_old, old$a, old$meanDisp, ceiling_new, new$a, new$meanDisp,
    new$meanDisp / old$meanDisp))
}
