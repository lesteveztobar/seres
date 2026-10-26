# build_v10_design.R -- 2026-10-02 sensitivity campaign v10 (two calibration arms).
#
# Writes two designs to data/params/:
#   sensitivity_design_v10_main.rds    2 arms x 6 sites x 500 LHS points x 3 seeds = 18,000 runs
#   sensitivity_design_v10_descent.rds 2 arms x 6 sites x 9 p_germ levels x 5 seeds =    540 runs
# plus sensitivity_design_v10_provenance.rds (seed, ranges, solved intercepts).
#
# Every row carries `task` (array index) and the task order is INTERLEAVED by
# arm x site: round r holds task r of every one of the 12 cells, in a cell
# order shuffled per round under DESIGN_SEED. Task r of every cell holds the
# SAME LHS points, so any prefix of completed tasks is a balanced subsample
# on an identical point set across sites and arms.
#
# Arms: `current` = intercepts from realistic.rds; `corrected` = intercepts
# re-solved so that realised annual survival at reference conditions (RH 85%,
# site-mean shortwave, T < 23 C) equals 0.44/0.60/0.85 at the MEAN REALISED
# SIZE per stage (abundance-weighted, all timesteps, the six completed
# realistic runs -- logs/phase1_20261002/p1.rds). delta_z is NOT changed.
suppressPackageStartupMessages({
  .libPaths(c("~/R/marvin_libs", .libPaths()))
  library(lhs)
})
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

DESIGN_SEED <- 20261002L
SITES <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")
ARMS <- c("current", "corrected")
N_POINTS <- 500L; N_SEEDS_MAIN <- 3L; POINTS_PER_TASK <- 10L
PGERM_LEVELS <- c(1e-2, 3.2e-3, 1e-3, 3.2e-4, 1e-4, 3.2e-5, 1e-5, 3.2e-6, 1e-6)
N_SEEDS_DESC <- 5L; LEVELS_PER_TASK <- 3L

base <- readRDS(file.path(PARAMS_DIR, "realistic.rds"))

# ---- corrected intercepts -----------------------------------------------------
diag <- readRDS("logs/phase1_20261002/p1.rds")
sz <- diag$size[diag$size$window == "all_t", ]
zbar <- sapply(c(S = "S", J = "J", A = "A"), function(st) {
  x <- sz[sz$stage == st, ]; sum(x$n_ind * x$mean_z) / sum(x$n_ind)
})
target <- c(S = 0.44, J = 0.60, A = 0.85)
# climate offset e at reference, read off the engine's own survival_logit()
b0_cur <- c(S = base$beta0S, J = base$beta0J, A = base$beta0A)
e_ref <- sapply(names(b0_cur), function(st)
  qlogis(survival_logit(st, 0, 20, 85, 1, base$beta0S, base$beta0J, base$beta0A, 0)) - b0_cur[[st]])
b0_cor <- qlogis(target^(1 / 12)) - base$beta1 * zbar - e_ref
check <- sapply(names(b0_cor), function(st)
  survival_logit(st, zbar[[st]], 20, 85, 1, b0_cor[["S"]], b0_cor[["J"]], b0_cor[["A"]], base$beta1)^12)
stopifnot(max(abs(check - target)) < 1e-9)
realised_cur_at_zbar <- sapply(names(b0_cur), function(st)
  survival_logit(st, zbar[[st]], 20, 85, 1, b0_cur[["S"]], b0_cur[["J"]], b0_cur[["A"]], base$beta1)^12)
arm_b0 <- rbind(current = b0_cur, corrected = b0_cor)

# ---- interleaved task order ---------------------------------------------------
set.seed(DESIGN_SEED)
cells <- expand.grid(arm = ARMS, site = SITES, stringsAsFactors = FALSE)
interleave <- function(n_rounds) {
  do.call(rbind, lapply(seq_len(n_rounds), function(r) {
    o <- cells[sample.int(nrow(cells)), ]
    data.frame(round = r, arm = o$arm, site = o$site, stringsAsFactors = FALSE)
  }))
}
add_arm <- function(d) {
  d$beta0S <- arm_b0[d$arm, "S"]; d$beta0J <- arm_b0[d$arm, "J"]; d$beta0A <- arm_b0[d$arm, "A"]
  d
}
site_abbr <- setNames(SITES, SITES)

# ---- main block: one LHS shared by every site and arm ------------------------
u <- randomLHS(N_POINTS, 5)
pts <- data.frame(
  p_poll = 0.02 + u[, 1] * (0.545 - 0.02),
  p_germ = 10^(log10(1e-4) + u[, 2] * (log10(1.2e-2) - log10(1e-4))),   # log-uniform
  p_s1 = 0.05 + u[, 3] * (0.50 - 0.05),
  S = 5.0e5 + u[, 4] * (3.9e6 - 5.0e5),
  n_founders = as.integer(pmin(499, 31 + floor(u[, 5] * (499 - 31 + 1))))  # uniform integer 31..499
)
pts$point <- seq_len(N_POINTS)
blocks <- split(sample.int(N_POINTS), rep(seq_len(N_POINTS / POINTS_PER_TASK), each = POINTS_PER_TASK))
order_main <- interleave(length(blocks))
order_main$task <- seq_len(nrow(order_main))
main <- do.call(rbind, lapply(seq_len(nrow(order_main)), function(i) {
  o <- order_main[i, ]
  pp <- pts[blocks[[o$round]], ]
  d <- do.call(rbind, lapply(seq_len(N_SEEDS_MAIN), function(s) cbind(pp, rep = s, seed = s)))
  d$task <- o$task; d$arm <- o$arm; d$site <- o$site; d$block <- "main"
  d$point_id <- sprintf("%s_%04d", site_abbr[[o$site]], d$point)
  d$run_id <- sprintf("v10m_%s_%s_%04d_r%d", o$arm, site_abbr[[o$site]], d$point, d$rep)
  d
}))
main <- add_arm(main)

# ---- descent block --------------------------------------------------------------
lvl_groups <- split(seq_along(PGERM_LEVELS), rep(seq_len(LEVELS_PER_TASK), length.out = length(PGERM_LEVELS)))
order_desc <- interleave(length(lvl_groups))
order_desc$task <- seq_len(nrow(order_desc))
desc <- do.call(rbind, lapply(seq_len(nrow(order_desc)), function(i) {
  o <- order_desc[i, ]
  g <- expand.grid(level = lvl_groups[[o$round]], rep = seq_len(N_SEEDS_DESC))
  data.frame(point = g$level, p_poll = base$p_poll, p_germ = PGERM_LEVELS[g$level], p_s1 = base$p_s1,
             S = 1.76e6,   # engine default (realistic.rds carries no S)
             n_founders = as.integer(base$n_founders), rep = g$rep, seed = g$rep,
             task = o$task, arm = o$arm, site = o$site, block = "descent",
             point_id = sprintf("%s_g%d", site_abbr[[o$site]], g$level),
             run_id = sprintf("v10d_%s_%s_g%d_r%d", o$arm, site_abbr[[o$site]], g$level, g$rep),
             stringsAsFactors = FALSE)
}))
desc <- add_arm(desc)

stopifnot(nrow(main) == 18000, nrow(desc) == 540, !anyDuplicated(main$run_id), !anyDuplicated(desc$run_id))
saveRDS(main, file.path(PARAMS_DIR, "sensitivity_design_v10_main.rds"))
saveRDS(desc, file.path(PARAMS_DIR, "sensitivity_design_v10_descent.rds"))
prov <- list(design_seed = DESIGN_SEED, zbar = zbar, e_ref = e_ref, target = target,
             beta0 = arm_b0, realised_current_at_zbar = realised_cur_at_zbar,
             realised_current_at_mid = diag$surv_mid, check_corrected = check,
             order_main = order_main, order_desc = order_desc, lhs_points = pts)
saveRDS(prov, file.path(PARAMS_DIR, "sensitivity_design_v10_provenance.rds"))

cat("zbar (pooled, abundance-weighted):\n"); print(round(zbar, 4))
cat("e_ref:\n"); print(round(e_ref, 4))
cat("intercepts:\n"); print(round(arm_b0, 4))
cat("realised current at zbar:\n"); print(round(realised_cur_at_zbar, 4))
cat("corrected check:\n"); print(check)
cat(sprintf("main: %d runs, %d tasks; descent: %d runs, %d tasks\n",
            nrow(main), max(main$task), nrow(desc), max(desc$task)))
print(sapply(pts[, 1:5], range))
cat(sprintf("p_germ median %.3g; frac <= 1e-3: %.3f\n", median(pts$p_germ), mean(pts$p_germ <= 1e-3)))
print(head(order_main, 13))
