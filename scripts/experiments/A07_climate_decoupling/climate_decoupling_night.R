# A07 night-time extension (2026-10-01) -- companion to climate_decoupling.R,
# which restricts everything to daylight (swdown > 0) because swdown = 0 at
# night has no variance: a swdown~X correlation is undefined there, not
# merely noisy. That result set is left untouched; this script asks the two
# questions night-time data CAN answer, so they can be set against it:
#
#   Option 1 -- temp~relhum covariation, computed identically for NIGHT and
#     DAY hours (overall / by month / by absolute 5m band / by relative
#     decile), pooled hourly x height exactly as climate_decoupling.R pools,
#     so the day/night contrast is like-for-like.
#
#   Option 2 -- carry-over: does the PREVIOUS day's insolation predict the
#     following night's microclimate (night mean temp/relhum, min temp, max
#     relhum)? The previous day has to come from ERA5 (microenv$.weather,
#     full 8760 h record, ABOVE-canopy swdown), not from the height files:
#     those are subsampled to every other day, and a UTC date holds the night
#     (UTC 00-11 = local 19:00-06:00, Ecuador UTC-5) BEFORE its own daylight
#     (UTC 12-23), so the day preceding any stored night is never stored.
#     Unit of analysis is the NIGHT (n ~ 183 per site), heights averaged
#     within a night -- not hours x heights, which would pseudo-replicate one
#     night many times. Reported raw AND as month-anomalies (each variable
#     minus its site-month mean), since dry-season sunny days / clear cold
#     nights would otherwise make the raw rho partly a seasonal artefact.
#
# Caveats carried into the outputs: these are microclimf model outputs
# forced by the same ERA5 record (not measurements), and within-canopy
# longwave -- the direct night-time driver -- is not stored in the height
# files (limitations: future version).
source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
source("scripts/02_model/engine/get_colonization.R")
library(ggplot2)

EXCLUDE_SALOYA <- Sys.getenv("EXCLUDE_SALOYA", unset = "0") == "1"
out_suffix <- if (EXCLUDE_SALOYA) "_noSaloya" else ""
SITES <- c("Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu", "LaElenita")
if (EXCLUDE_SALOYA) SITES <- setdiff(SITES, "Saloya")
LOCAL_OFFSET_H <- -5L   # Ecuador (no DST)
MIN_NIGHT_HOURS <- 8L   # a full night is ~11 h; drop truncated record edges
OUT <- file.path(OUTPUT_DIR, "climate_decoupling_night")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

niches <- load_observations()
niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                 !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))

spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10 || length(unique(x[ok])) < 3 || length(unique(y[ok])) < 3) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
rho_row <- function(x, y, ...) { r <- spearman_safe(x, y); data.frame(..., rho = r["rho"], p = r["p"], n = r["n"], row.names = NULL) }

# ── Hourly table WITH timestamps (climate_decoupling_table.rds keeps month only) ──
hourly <- list(); prevday <- list()
for (site in SITES) {
  microenv_path <- file.path(PROCESSED_DIR, sprintf("microenv_%s_h0.40.rds", site))
  if (!file.exists(microenv_path)) { message("Skipping ", site); next }
  message("Reading ", site, "...")
  microenv <- readRDS(microenv_path)
  heights <- microenv_heights(microenv)
  ceiling <- site_canopy_ceiling(site, niches, max(heights))
  px <- do.call(rbind, parallel::mclapply(heights, function(h) {
    hd <- load_height(microenv, h)
    if (is.null(hd)) return(NULL)
    data.frame(height = h, tme = hd$tme, temp = hd$temp_mean, relhum = hd$relhum_mean,
               swdown = hd$Rdirdown_mean + hd$Rdifdown_mean)
  }, mc.cores = n_cores))
  px$site <- site
  px$rel_height <- px$height / ceiling
  loc <- px$tme + LOCAL_OFFSET_H * 3600
  px$month <- as.integer(format(loc, "%m"))
  px$daypart <- ifelse(px$swdown > 0, "day", "night")
  # night_of = local calendar date on which the night STARTED (evening), so
  # 19:00 Jan 4 .. 06:00 Jan 5 are all night_of Jan 4 -- the day it follows.
  ld <- as.Date(format(loc, "%Y-%m-%d"))
  px$night_of <- ifelse(as.integer(format(loc, "%H")) >= 12L, ld, ld - 1L)
  px$height_band_abs <- cut(px$height, breaks = seq(0, ceiling(max(px$height) / 5) * 5, by = 5), include.lowest = TRUE)
  px$height_band_rel <- cut(px$rel_height, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  hourly[[site]] <- px

  # ERA5 above-canopy daily insolation per LOCAL date (kWh m-2 d-1; hourly W m-2 x 1 h)
  w <- microenv$.weather
  wl <- as.Date(format(w$obs_time + LOCAL_OFFSET_H * 3600, "%Y-%m-%d"))
  ins <- tapply(w$swdown, wl, function(v) if (sum(is.finite(v)) >= 20) sum(v, na.rm = TRUE) / 1000 else NA_real_)
  prevday[[site]] <- data.frame(site = site, night_of = as.numeric(as.Date(names(ins))), insolation_kwh = as.numeric(ins))
  rm(microenv); gc(verbose = FALSE)
}
all_df <- do.call(rbind, hourly); rownames(all_df) <- NULL
saveRDS(all_df, file.path(PROCESSED_DIR, sprintf("climate_decoupling_table_hourly_tme%s.rds", out_suffix)))
cat("Total rows:", nrow(all_df), " night share:", round(mean(all_df$daypart == "night"), 3), "\n\n")

# ── Option 1: temp~relhum, night vs day, like-for-like ──────────────────────
opt1 <- list()
for (s in names(hourly)) for (dp in c("night", "day")) {
  d <- all_df[all_df$site == s & all_df$daypart == dp, ]
  opt1[[length(opt1) + 1]] <- rho_row(d$temp, d$relhum, site = s, daypart = dp, stratum = "overall", level = "all")
  for (m in sort(unique(d$month))) { dm <- d[d$month == m, ]
    opt1[[length(opt1) + 1]] <- rho_row(dm$temp, dm$relhum, site = s, daypart = dp, stratum = "month", level = as.character(m)) }
  for (col in c("height_band_abs", "height_band_rel")) for (b in levels(d[[col]])) {
    db <- d[!is.na(d[[col]]) & d[[col]] == b, ]
    opt1[[length(opt1) + 1]] <- rho_row(db$temp, db$relhum, site = s, daypart = dp,
                                        stratum = if (col == "height_band_abs") "absolute_5m" else "relative_decile", level = b) }
}
opt1 <- do.call(rbind, opt1)
write.csv(opt1, file.path(OUT, sprintf("temp_relhum_night_vs_day%s.csv", out_suffix)), row.names = FALSE)
cat("=== Option 1: temp~relhum Spearman, overall, night vs day ===\n")
print(opt1[opt1$stratum == "overall", c("site", "daypart", "rho", "n")], row.names = FALSE)

# ── Option 2: previous-day insolation -> following night ────────────────────
nt <- all_df[all_df$daypart == "night", ]
agg_night <- function(d, keys) {
  g <- interaction(d[keys], drop = TRUE, sep = "\r")
  out <- data.frame(do.call(rbind, strsplit(levels(g), "\r", fixed = TRUE)), stringsAsFactors = FALSE)
  names(out) <- keys
  sp <- split(d, g)
  out$n_hours  <- vapply(sp, function(x) length(unique(x$tme)), 1)
  out$temp_mean <- vapply(sp, function(x) mean(x$temp), 1)
  out$relhum_mean <- vapply(sp, function(x) mean(x$relhum), 1)
  # min/max over the night's hourly HEIGHT-MEAN series, so they stay per-night
  out$temp_min <- vapply(sp, function(x) min(tapply(x$temp, x$tme, mean)), 1)
  out$relhum_max <- vapply(sp, function(x) max(tapply(x$relhum, x$tme, mean)), 1)
  out$night_of <- as.numeric(out$night_of)
  out
}
pd <- do.call(rbind, prevday)
resp <- c("temp_mean", "relhum_mean", "temp_min", "relhum_max")
opt2 <- list(); nights_all <- list()
for (strat in c("overall", "absolute_5m", "relative_decile")) {
  d <- nt
  d$stratum_level <- switch(strat, overall = "all",
                            absolute_5m = as.character(d$height_band_abs),
                            relative_decile = as.character(d$height_band_rel))
  d <- d[!is.na(d$stratum_level), ]
  a <- agg_night(d, c("site", "stratum_level", "night_of"))
  a <- a[a$n_hours >= MIN_NIGHT_HOURS, ]
  a <- merge(a, pd, by = c("site", "night_of"))
  a$month <- as.integer(format(as.Date(a$night_of, origin = "1970-01-01"), "%m"))
  # month anomalies within site x stratum level
  for (v in c("insolation_kwh", resp)) {
    a[[paste0(v, "_anom")]] <- a[[v]] - ave(a[[v]], a$site, a$stratum_level, a$month, FUN = function(z) mean(z, na.rm = TRUE))
  }
  a$stratum <- strat
  nights_all[[strat]] <- a
  for (key in unique(paste(a$site, a$stratum_level, sep = "\r"))) {
    s <- strsplit(key, "\r", fixed = TRUE)[[1]]
    da <- a[a$site == s[1] & a$stratum_level == s[2], ]
    for (v in resp) {
      opt2[[length(opt2) + 1]] <- rho_row(da$insolation_kwh, da[[v]], site = s[1], stratum = strat, level = s[2], response = v, version = "raw")
      opt2[[length(opt2) + 1]] <- rho_row(da$insolation_kwh_anom, da[[paste0(v, "_anom")]], site = s[1], stratum = strat, level = s[2], response = v, version = "month_anomaly")
    }
  }
}
opt2 <- do.call(rbind, opt2)
write.csv(opt2, file.path(OUT, sprintf("prevday_insolation_vs_night%s.csv", out_suffix)), row.names = FALSE)
saveRDS(do.call(rbind, nights_all), file.path(PROCESSED_DIR, sprintf("climate_decoupling_night_table%s.rds", out_suffix)))
cat("\n=== Option 2: ERA5 previous-day insolation vs following night (overall, n = nights) ===\n")
print(opt2[opt2$stratum == "overall", c("site", "response", "version", "rho", "p", "n")], row.names = FALSE)

# ── Figures ─────────────────────────────────────────────────────────────────
o1 <- opt1[opt1$stratum == "overall", ]
p1 <- ggplot(o1, aes(x = site, y = rho, fill = daypart)) +
  geom_col(position = "dodge") + geom_hline(yintercept = 0, linetype = "dashed") + coord_flip() +
  labs(title = "temp~relhum covariation: night vs day (Spearman rho, pooled hours x heights)", y = "Spearman rho", x = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT, sprintf("temp_relhum_night_vs_day_summary%s.png", out_suffix)), p1, width = 9, height = 6, dpi = 300, bg = "white")

o1r <- opt1[opt1$stratum == "relative_decile", ]
o1r$level <- factor(o1r$level, levels = levels(all_df$height_band_rel))
p1b <- ggplot(o1r, aes(x = level, y = rho, colour = daypart, group = daypart)) +
  geom_line() + geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + facet_wrap(~site) +
  labs(title = "temp~relhum by relative height decile, night vs day", x = "relative height decile", y = "Spearman rho") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, sprintf("temp_relhum_by_rel_height%s.png", out_suffix)), p1b, width = 12, height = 8, dpi = 300, bg = "white")

o2 <- opt2[opt2$stratum == "overall", ]
p2 <- ggplot(o2, aes(x = site, y = rho, fill = version)) +
  geom_col(position = "dodge") + geom_hline(yintercept = 0, linetype = "dashed") + coord_flip() + facet_wrap(~response) +
  labs(title = "Previous-day ERA5 insolation vs following night (Spearman rho, n = nights)", y = "Spearman rho", x = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(OUT, sprintf("prevday_insolation_summary%s.png", out_suffix)), p2, width = 11, height = 7, dpi = 300, bg = "white")

na <- nights_all[["overall"]]
p2b <- ggplot(na, aes(x = insolation_kwh, y = relhum_mean)) +
  geom_point(alpha = 0.4, size = 0.8) + geom_smooth(method = "loess", se = FALSE) + facet_wrap(~site, scales = "free_y") +
  labs(title = "Previous-day insolation vs following-night mean relhum", x = "ERA5 insolation (kWh m-2 d-1)", y = "night mean relhum (%)") +
  theme_minimal(base_size = 10)
ggsave(file.path(OUT, sprintf("prevday_insolation_vs_night_relhum%s.png", out_suffix)), p2b, width = 12, height = 8, dpi = 300, bg = "white")
p2c <- p2b + aes(y = temp_min) + labs(title = "Previous-day insolation vs following-night minimum temp", y = "night min temp (C)")
ggsave(file.path(OUT, sprintf("prevday_insolation_vs_night_tempmin%s.png", out_suffix)), p2c, width = 12, height = 8, dpi = 300, bg = "white")
cat("\nOutputs in", OUT, "\nDone.\n")
