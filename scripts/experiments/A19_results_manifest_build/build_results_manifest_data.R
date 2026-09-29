# build_results_manifest_data.R -- extracts headline numbers from every
# completed Phase F output file into one CSV, for results_manifest.md to
# be built from (no interpretation here, just extraction).
source("scripts/02_model/config/paths.R")

SITES <- c("Maquipucuna", "Mashpi", "Yanayacu", "MindoMirador", "MindoTarabita", "Saloya")
EXP_TAGS <- c("best_case", "realistic", "founder_number",
             "horizontal_resolution_pooled", "horizontal_resolution_voxel",
             "realistic_273founders")

extract_stats <- function(result) {
  if (is.null(result) || is.null(result$summary)) return(NULL)
  s <- result$summary
  final_t <- max(s$t)
  final <- s[s$t == final_t, ]
  n_reps <- length(unique(s$rep))
  data.frame(
    n_reps = n_reps,
    n_recruited = sum(final$recruited, na.rm = TRUE),
    frac_recruited = mean(final$recruited, na.rm = TRUE),
    n_not_extinct = sum(!final$extinct, na.rm = TRUE),
    frac_not_extinct = mean(!final$extinct, na.rm = TRUE),
    mean_final_total = mean(final$total, na.rm = TRUE),
    sd_final_total = sd(final$total, na.rm = TRUE)
  )
}

rows <- list()
for (site in SITES) {
  for (tag in EXP_TAGS) {
    f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_%s_h0.40.rds", site, tag))
    if (!file.exists(f)) next
    result <- tryCatch(readRDS(f), error = function(e) NULL)
    st <- extract_stats(result)
    if (is.null(st)) next
    rows[[length(rows) + 1]] <- data.frame(site = site, experiment = tag, st, file = f)
  }
  # founder_number is a sweep -- also break out by n_founders level
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_founder_number_h0.40.rds", site))
  if (file.exists(f)) {
    result <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(result) && !is.null(result$summary) && "n_founders" %in% names(result$summary)) {
      s <- result$summary
      final_t <- max(s$t)
      final <- s[s$t == final_t, ]
      by_nf <- do.call(rbind, lapply(sort(unique(final$n_founders)), function(nf) {
        d <- final[final$n_founders == nf, ]
        data.frame(site = site, experiment = "founder_number_sweep", n_founders = nf,
                  n_reps = nrow(d), frac_recruited = mean(d$recruited, na.rm = TRUE),
                  mean_final_total = mean(d$total, na.rm = TRUE), file = f)
      }))
      rows[[length(rows) + 1]] <- by_nf
    }
  }
}
out <- do.call(rbind, lapply(rows, function(r) {
  # normalize columns across the two row shapes (overall vs by-founder-level)
  cols <- c("site", "experiment", "n_founders", "n_reps", "n_recruited", "frac_recruited",
           "n_not_extinct", "frac_not_extinct", "mean_final_total", "sd_final_total", "file")
  for (c in cols) if (!c %in% names(r)) r[[c]] <- NA
  r[cols]
}))
write.csv(out, "output/results_manifest_data.csv", row.names = FALSE)
cat("Wrote output/results_manifest_data.csv --", nrow(out), "rows\n")
print(out, row.names = FALSE)
