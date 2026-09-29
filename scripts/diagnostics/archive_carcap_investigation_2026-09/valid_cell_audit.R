# valid_cell_audit.R -- per-site, old vs new: valid ERA5 grid cells and the
# finite fraction of the resulting per-pixel voxel_quantiles array, per
# variable, at one representative height tier (the lowest, h0.10 -- ground
# level, present at every site).
source("scripts/02_model/config/patches.R")
source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

sites <- c("Maquipucuna","Mashpi","Yanayacu","MindoMirador","MindoTarabita","Saloya","LaElenita")
vars  <- c("temp","relhum","windspeed","swdown","difrad")

finite_frac <- function(manifest_path, var) {
  m <- readRDS(manifest_path)
  h0 <- file.path(m$.height_dir, sprintf("h%.2f.rds", m$.heights[1]))
  if (!file.exists(h0)) return(NA_real_)
  hd <- readRDS(h0)
  key <- if (var %in% c("temp","relhum")) sprintf("annual_both_%s", var) else
         if (var == "windspeed") "annual_both_windspeed" else
         sprintf("annual_day_%s", var)
  q <- hd$voxel_quantiles$quantiles[[key]]
  if (is.null(q)) return(NA_real_)
  mean(is.finite(q))
}

valid_cells_log_line <- function(site, tag) {
  pat <- if (tag == "new") sprintf("logs/microclim_%s_20260907_11", site) else NULL
  files <- Sys.glob(sprintf("logs/microclim_%s_*.log", site))
  files <- files[!grepl("20260907_11", files)]  # exclude this week's new-run logs for "old"
  if (tag == "new") {
    f <- Sys.glob(sprintf("logs/microclim_%s_20260907_11*.log", site))
    f <- f[1]
  } else {
    f <- tail(sort(files), 1)  # most recent OLD run's own log, i.e. the one behind the archived manifest
  }
  if (length(f) == 0 || is.na(f)) return(NA_character_)
  ln <- grep("Valid grid cells", readLines(f), value = TRUE)
  if (length(ln) == 0) return(NA_character_)
  sub(".*Valid grid cells: ", "", ln[1])
}

rows <- list()
for (s in sites) {
  new_path <- sprintf("data/processed/microenv_%s_h0.40.rds", s)
  old_path <- sprintf("data/processed/archive_pre_v7pix/microenv_%s_h0.40.rds", s)
  vc_new <- valid_cells_log_line(s, "new")
  vc_old <- valid_cells_log_line(s, "old")
  for (v in vars) {
    rows[[length(rows) + 1]] <- data.frame(
      site = s, variable = v,
      valid_cells_old = vc_old, valid_cells_new = vc_new,
      finite_frac_old = finite_frac(old_path, v),
      finite_frac_new = finite_frac(new_path, v)
    )
  }
}
out <- do.call(rbind, rows)
write.csv(out, "output/valid_cell_audit.csv", row.names = FALSE)
print(out, row.names = FALSE)
