# debug_founder_crash2.R -- instrumented copy of run_pass3_survive_grow()
# to pinpoint exactly which intermediate quantity first goes non-finite,
# causing "NAs are not allowed in subscripted assignments" (see
# debug_founder_crash.R for the first reproduction + traceback).
options(error = function() { traceback(2); quit(save = "no", status = 1) })

site_name <- "Maquipucuna"
params_file <- "data/params/n_founders.rds"
height_step <- 0.4

library(parallel)
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")

manifest_suffix <- if (height_step != 0.1) sprintf("_h%.2f", height_step) else ""
microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s%s.rds", site_name, manifest_suffix))
microenv <- readRDS(microenv_path)
available_heights <- microenv_heights(microenv)

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
mean_canopy <- mean(niches$CanopyHeight_m[niches$Area_or_Site == site_name], na.rm = TRUE)
canopy_grid <- matrix(mean_canopy, nrow = 50, ncol = 50)
site <- list(Site = site_name)
forestparams <- default_forestparams()

params <- readRDS(params_file)
params$canopy_z <- mean_canopy
params$n_founders <- 30

log_msg <- function(msg) cat(msg, "\n")

clim_cache <- build_clim_cache(microenv)

# ── Instrumented override of run_pass3_survive_grow ──────────────────────────
run_pass3_survive_grow <- function(state, abundanceS, abundanceJ, abundanceA,
                                   size_S, size_J, size_A, fruited, t, stochastic = FALSE) {
  p <- state$params
  xDim <- state$xDim
  yDim <- state$yDim
  cat(sprintf("[pass3] entry: anyNA(abundanceS)=%s anyNA(abundanceJ)=%s anyNA(abundanceA)=%s anyNA(size_S)=%s anyNA(size_J)=%s anyNA(size_A)=%s\n",
    anyNA(abundanceS), anyNA(abundanceJ), anyNA(abundanceA), anyNA(size_S), anyNA(size_J), anyNA(size_A)))

  .surv_slice_vec <- function(stage_slice, prob_arr) {
    if (sum(stage_slice) == 0L) return(stage_slice * 0L)
    prob_arr <- pmin(pmax(prob_arr, 0), 1)
    prob_arr[!is.finite(prob_arr)] <- 0
    array(rbinom(length(stage_slice), as.integer(stage_slice), prob_arr), dim = dim(stage_slice))
  }

  for (sp in seq_len(state$n_species)) {
    for (zi in seq_len(state$zDim)) {
      if (sum(abundanceS[, , zi, t, sp]) + sum(abundanceJ[, , zi, t, sp]) +
        sum(abundanceA[, , zi, t, sp]) == 0L) next

      clim_year <- state$clim_by_height[[zi]]
      if (is.null(clim_year)) next
      precip_annual <- mean(clim_year$precip, na.rm = TRUE) * 8760
      if (!is.finite(precip_annual)) {
        cat(sprintf("*** sp=%d zi=%d (height=%.2f): precip_annual is NOT FINITE (%.4f). anyNA(clim_year$precip)=%s all-NA=%s n=%d\n",
          sp, zi, available_heights[zi], precip_annual, anyNA(clim_year$precip),
          all(is.na(clim_year$precip)), length(clim_year$precip)))
      }

      abundance_S_slice <- abundanceS[, , zi, t, sp]
      abundance_J_slice <- abundanceJ[, , zi, t, sp]
      abundance_A_slice <- abundanceA[, , zi, t, sp]
      size_S_slice <- size_S[, , zi, t, sp]
      size_J_slice <- size_J[, , zi, t, sp]
      size_A_slice <- size_A[, , zi, t, sp]
      delta_s_total_S <- array(0, dim = c(xDim, yDim))
      delta_s_total_J <- array(0, dim = c(xDim, yDim))
      delta_s_total_A <- array(0, dim = c(xDim, yDim))

      for (month in 1:12) {
        clim_month_table <- state$clim_month_by_height[[zi]][[month]]
        if (is.null(clim_month_table) || nrow(clim_month_table) == 0) next

        temp_mat <- .fill_na(.clim_voxel_slice(state, zi, month, "day", "temp", stochastic),
          mean(clim_month_table$temp, na.rm = TRUE), xDim, yDim)
        relhum_mat <- .fill_na(.clim_voxel_slice(state, zi, month, "day", "relhum", stochastic),
          mean(clim_month_table$relhum, na.rm = TRUE), xDim, yDim)
        swdown_mat <- .fill_na(.clim_voxel_slice(state, zi, month, "day", "swdown", stochastic),
          mean(clim_month_table$swdown[clim_month_table$swdown > 0], na.rm = TRUE), xDim, yDim)
        swdown_rel <- if (p$mean_swdown_site > 0) swdown_mat / p$mean_swdown_site else matrix(1.0, xDim, yDim)
        swdown_rel[!is.finite(swdown_rel)] <- 1.0

        if (anyNA(temp_mat) || anyNA(relhum_mat) || anyNA(swdown_rel)) {
          cat(sprintf("*** sp=%d zi=%d month=%d: anyNA temp=%s relhum=%s swdown_rel=%s\n",
            sp, zi, month, anyNA(temp_mat), anyNA(relhum_mat), anyNA(swdown_rel)))
        }

        survive_prob_S <- survival_logit("S", size_S_slice, temp_mat, relhum_mat, swdown_rel, p$beta0_S, p$beta0_J, p$beta0_A, p$beta1)
        survive_prob_J <- survival_logit("J", size_J_slice, temp_mat, relhum_mat, swdown_rel, p$beta0_S, p$beta0_J, p$beta0_A, p$beta1)
        survive_prob_A <- survival_logit("A", size_A_slice, temp_mat, relhum_mat, swdown_rel, p$beta0_S, p$beta0_J, p$beta0_A, p$beta1)
        if (anyNA(survive_prob_S) || anyNA(survive_prob_J) || anyNA(survive_prob_A)) {
          cat(sprintf("*** sp=%d zi=%d month=%d: anyNA survive_prob S=%s J=%s A=%s\n",
            sp, zi, month, anyNA(survive_prob_S), anyNA(survive_prob_J), anyNA(survive_prob_A)))
        }

        cat(sprintf("[pre-A-surv] sp=%d zi=%d month=%d: range(abundance_A_slice)=[%s,%s] range(survive_prob_A)=[%s,%s] sum(abundance_A_slice)=%s\n",
          sp, zi, month,
          suppressWarnings(min(abundance_A_slice)), suppressWarnings(max(abundance_A_slice)),
          suppressWarnings(min(survive_prob_A)), suppressWarnings(max(survive_prob_A)),
          sum(abundance_A_slice)))

        abundance_S_slice <- .surv_slice_vec(abundance_S_slice, survive_prob_S)
        abundance_J_slice <- .surv_slice_vec(abundance_J_slice, survive_prob_J)
        abundance_A_slice <- .surv_slice_vec(abundance_A_slice, survive_prob_A)
        if (anyNA(abundance_S_slice) || anyNA(abundance_J_slice) || anyNA(abundance_A_slice)) {
          cat(sprintf("*** sp=%d zi=%d month=%d: anyNA POST-SURVIVAL abundance S=%s J=%s A=%s\n",
            sp, zi, month, anyNA(abundance_S_slice), anyNA(abundance_J_slice), anyNA(abundance_A_slice)))
          if (anyNA(abundance_A_slice)) {
            cat("abundance_A_slice (post):\n"); print(abundance_A_slice)
          }
        }

        p_StoJ <- transition_logit("S", precip_annual, relhum_mat, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
        p_JtoA <- transition_logit("J", precip_annual, relhum_mat, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
        epsilon <- rnorm(1, 0, p$sigma)
        p_StoJ <- pmin(1, pmax(0, p_StoJ + epsilon))
        p_JtoA <- pmin(1, pmax(0, p_JtoA + epsilon))
        if (anyNA(p_StoJ) || anyNA(p_JtoA)) {
          cat(sprintf("*** sp=%d zi=%d month=%d: anyNA p_StoJ=%s p_JtoA=%s precip_annual=%.4f\n",
            sp, zi, month, anyNA(p_StoJ), anyNA(p_JtoA), precip_annual))
        }

        n_StoJ <- .surv_slice_vec(abundance_S_slice, p_StoJ)
        n_JtoA <- .surv_slice_vec(abundance_J_slice, p_JtoA)
        if (anyNA(n_StoJ) || anyNA(n_JtoA)) {
          cat(sprintf("*** sp=%d zi=%d month=%d: anyNA n_StoJ=%s n_JtoA=%s\n",
            sp, zi, month, anyNA(n_StoJ), anyNA(n_JtoA)))
          cat("n_JtoA values:\n"); print(n_JtoA)
          cat("abundance_J_slice (input to n_JtoA):\n"); print(abundance_J_slice)
          cat("p_JtoA:\n"); print(p_JtoA)
        }

        a_denom <- abundance_A_slice + n_JtoA
        if (anyNA(a_denom)) {
          cat(sprintf("*** sp=%d zi=%d month=%d: a_denom HAS NA -- STOPPING HERE\n", sp, zi, month))
          stop("Located the NA source -- see diagnostics above.")
        }
        a_mask <- a_denom > 0L
        a_blend <- (abundance_A_slice * size_A_slice + n_JtoA * pmax(p$s_A_min, size_J_slice)) / pmax(a_denom, 1)
        size_A_slice[a_mask] <- a_blend[a_mask]

        J_stayers <- abundance_J_slice - n_JtoA
        j_denom <- J_stayers + n_StoJ
        j_mask <- j_denom > 0L
        j_blend <- (J_stayers * size_J_slice + n_StoJ * size_S_slice) / pmax(j_denom, 1)
        size_J_slice[j_mask] <- j_blend[j_mask]

        abundance_S_slice <- abundance_S_slice - n_StoJ
        abundance_J_slice <- abundance_J_slice + n_StoJ - n_JtoA
        abundance_A_slice <- abundance_A_slice + n_JtoA

        delta_size_base <- (p$delta_s_base / 12) * (precip_annual / 2500) * (relhum_mat / 85)
        noise_sd <- p$sigma * 0.5 / sqrt(12)
        delta_size_S <- delta_size_base + rnorm(xDim * yDim, 0, noise_sd); dim(delta_size_S) <- c(xDim, yDim)
        delta_size_J <- delta_size_base + rnorm(xDim * yDim, 0, noise_sd); dim(delta_size_J) <- c(xDim, yDim)
        delta_size_A <- delta_size_base + rnorm(xDim * yDim, 0, noise_sd); dim(delta_size_A) <- c(xDim, yDim)
        delta_s_total_S <- delta_s_total_S + delta_size_S
        delta_s_total_J <- delta_s_total_J + delta_size_J
        delta_s_total_A <- delta_s_total_A + delta_size_A
      }

      abundanceS[, , zi, t, sp] <- abundance_S_slice
      abundanceJ[, , zi, t, sp] <- abundance_J_slice
      abundanceA[, , zi, t, sp] <- abundance_A_slice

      fruited_slice <- fruited[, , zi, sp]
      delta_s_total_A[fruited_slice] <- delta_s_total_A[fruited_slice] * p$cost_repro

      size_S_new <- pmin(p$s_S_max, pmax(p$s_S_min, size_S_slice + delta_s_total_S))
      size_J_new <- pmin(p$s_J_max, pmax(p$s_J_min, size_J_slice + delta_s_total_J))
      size_A_new <- pmin(p$s_A_max, pmax(p$s_A_min, size_A_slice + delta_s_total_A))

      has_S <- abundanceS[, , zi, t, sp] > 0L
      has_J <- abundanceJ[, , zi, t, sp] > 0L
      has_A <- abundanceA[, , zi, t, sp] > 0L
      size_S_slice[has_S] <- size_S_new[has_S]
      size_J_slice[has_J] <- size_J_new[has_J]
      size_A_slice[has_A] <- size_A_new[has_A]
      size_S[, , zi, t, sp] <- size_S_slice
      size_J[, , zi, t, sp] <- size_J_slice
      size_A[, , zi, t, sp] <- size_A_slice
    }
  }
  list(S = abundanceS, J = abundanceJ, A = abundanceA, size_S = size_S, size_J = size_J, size_A = size_A)
}

cat("\n=== Calling runcolonization() directly (no tryCatch), instrumented pass3 ===\n")
r <- runcolonization(
  site = site, niches = niches, canopy_grid = canopy_grid, microenv = microenv,
  timesteps = 50, resolution = 10, carCap = 5, maxDisp = 10, spinup = 5,
  Visualize = FALSE, parameters = params, forestparams = forestparams,
  clim_cache = clim_cache, seed = 1
)
cat("SUCCEEDED -- no crash reproduced with this seed/value.\n")
