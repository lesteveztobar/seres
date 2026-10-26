# phase1_survival_size_diag.R -- 2026-10-02 pre-campaign diagnostics.
# (a) realised annual survival / transition at reference conditions, computed
#     through the engine's own survival_logit()/transition_logit();
# (b) per-stage realised size z in the completed realistic runs;
# (c) run_seconds by regime across the v9 sensitivity parts.
# Read-only: writes one RDS of results to the path given as argument 1.
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")
out_rds <- commandArgs(trailingOnly = TRUE)[1]
res <- list()

p <- readRDS(file.path(PARAMS_DIR, "realistic.rds"))
res$params <- p[c("beta0S", "beta0J", "beta0A", "beta1", "psi0S", "psi0J",
                  "beta_precip", "beta_rh", "sigma", "delta_s_base",
                  "s_S_min", "s_S_max", "s_J_min", "s_J_max", "s_A_min", "s_A_max")]

# ---- (a) survival at reference: RH 85, swdown_rel 1, temp 20 (< 23) --------
annual_surv <- function(stage, z, b0S = p$beta0S, b0J = p$beta0J, b0A = p$beta0A) {
  survival_logit(stage, z, temp = 20, relhum = 85, swdown_rel = 1,
                 beta0_S = b0S, beta0_J = b0J, beta0_A = b0A, beta1 = p$beta1)^12
}
mid <- c(S = 0.5, J = 4, A = 13.5)
res$surv_mid <- sapply(names(mid), function(s) annual_surv(s, mid[[s]]))
res$surv_z0_e0 <- c(S = plogis(p$beta0S)^12, J = plogis(p$beta0J)^12, A = plogis(p$beta0A)^12)
# linear-predictor climate offset e at reference, per stage
res$e_ref <- sapply(names(mid), function(s) qlogis(survival_logit(s, 0, 20, 85, 1, p$beta0S, p$beta0J, p$beta0A, 0)) -
                      c(S = p$beta0S, J = p$beta0J, A = p$beta0A)[[s]])
# seedling survival at the drift sizes quoted in the brief (5.26, 1.10) and at the class ceiling
res$surv_S_at <- sapply(c(0.5, 1.0, 1.10, 5.26), function(z) annual_surv("S", z))

# ---- transitions at reference RH 85 and P_ref = 2500 mm (engine's growth reference)
tr_annual <- function(stage, P, rh = 85, noise = TRUE, n = 2e5) {
  q <- transition_logit(stage, P, rh, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
  if (!noise) return(c(nominal = 1 - (1 - q)^12))
  set.seed(1)
  eps <- matrix(rnorm(n * 12, 0, p$sigma), n, 12)  # engine: one eps per month, clipped [0,1]
  qm <- eps; qm[] <- pmin(1, pmax(0, q + eps))
  c(nominal = 1 - (1 - q)^12, realised = mean(1 - apply(1 - qm, 1, prod)), q_month = q,
    q_month_eff = mean(qm))
}
res$trans_Pref <- rbind(S = tr_annual("S", 2500), J = tr_annual("J", 2500))
res$trans_z0 <- c(S = 1 - (1 - plogis(p$psi0S))^12, J = 1 - (1 - plogis(p$psi0J))^12)

# ---- (b) realised size per stage in completed realistic runs -----------------
sites <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu")
size_rows <- list(); precip_rows <- list()
for (site in sites) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_realistic_h0.40.rds", site))
  if (!file.exists(f)) next
  r <- readRDS(f)
  runs <- Filter(Negate(is.null), r$runs)
  for (i in seq_along(runs)) {
    ru <- runs[[i]]
    for (st in c("S", "J", "A")) {
      ab <- ru[[paste0("abundance", st)]]; sz <- ru[[paste0("size_", st)]]
      nt <- dim(ab)[4]
      for (tt in list(all = seq_len(nt), final = nt)) {
        a <- ab[, , , tt, , drop = FALSE]; s <- sz[, , , tt, , drop = FALSE]
        occ <- which(a > 0)
        w <- as.numeric(a[occ]); z <- as.numeric(s[occ])
        if (!length(w)) next
        o <- order(z); cw <- cumsum(w[o]) / sum(w)
        q <- sapply(c(0.05, 0.5, 0.95), function(pp) z[o][which(cw >= pp)[1]])
        size_rows[[length(size_rows) + 1]] <- data.frame(
          site = site, rep = i, stage = st, window = if (length(tt) > 1) "all_t" else "final_t",
          n_ind = sum(w), mean_z = sum(w * z) / sum(w), q05 = q[1], q50 = q[2], q95 = q[3],
          min_z = min(z), max_z = max(z), regime = ru$regime %||% NA)
      }
    }
    # site annual precip as the engine computes it (mean hourly precip * 8760), per height
    pa <- vapply(ru$state$clim_by_height, function(cl) if (is.null(cl)) NA_real_ else mean(cl$precip, na.rm = TRUE) * 8760, numeric(1))
    precip_rows[[site]] <- data.frame(site = site, precip_min = min(pa, na.rm = TRUE), precip_mean = mean(pa, na.rm = TRUE), precip_max = max(pa, na.rm = TRUE))
  }
  rm(r, runs); gc()
}
res$size <- do.call(rbind, size_rows)
res$precip <- do.call(rbind, precip_rows)

# ---- (c) run_seconds by regime, v9 parts -------------------------------------
parts <- list.files(OUTPUT_DIR, pattern = "^part_.*\\.csv$", recursive = TRUE, full.names = TRUE)
parts <- parts[grepl("sensitivity_runs_v9", parts)]
v9 <- do.call(rbind, lapply(parts, function(f) read.csv(f, stringsAsFactors = FALSE)[, c("site", "regime", "run_seconds", "p_germ")]))
res$v9_n <- nrow(v9)
res$v9_secs <- do.call(rbind, lapply(split(v9$run_seconds, v9$regime), function(x)
  c(n = length(x), mean = mean(x), median = median(x), p90 = unname(quantile(x, 0.9)), max = max(x))))
res$v9_pgerm <- c(min = min(v9$p_germ), median = median(v9$p_germ), frac_le_0.001 = mean(v9$p_germ <= 0.001))
saveRDS(res, out_rds)
print(res)
