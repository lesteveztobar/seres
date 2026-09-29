source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")

sites <- c("Maquipucuna","Mashpi","MindoMirador","MindoTarabita","Yanayacu")

main_effect_range <- function(d, param, t_max) {
  final <- d[d$t == t_max, ]
  agg <- aggregate(as.formula(paste("total ~", param)), data = final, FUN = mean)
  rng <- diff(range(agg$total))
  rel_rng <- rng / max(agg$total)  # relative range as fraction of max level's abundance
  list(abs_range = rng, rel_range = rel_rng, levels = nrow(agg))
}

cat("=== v4 (reproduction_factorial_v4_v6): per-parameter main-effect range ===\n")
v4_params <- c("p_poll","p_germ","p_s1","n_founders")
v4_summary <- list()
for (s in sites) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_reproduction_factorial_v4_v6_h0.40.rds", s))
  d <- readRDS(f)
  t_max <- max(d$t)
  cat(sprintf("--- %s (t_max=%d) ---\n", s, t_max))
  for (p in v4_params) {
    r <- main_effect_range(d, p, t_max)
    cat(sprintf("  %-12s abs_range=%8.1f  rel_range=%.3f  (%d levels)\n", p, r$abs_range, r$rel_range, r$levels))
    v4_summary[[p]] <- c(v4_summary[[p]], r$rel_range)
  }
}
cat("\n=== v4 mean relative-range across sites (ranking) ===\n")
v4_means <- sort(sapply(v4_summary, mean), decreasing = TRUE)
print(round(v4_means, 4))

cat("\n=== v5 (survival_factorial_v5_v6): per-parameter main-effect range ===\n")
v5_params <- c("beta0S","beta0J","beta0A")
v5_summary <- list()
for (s in sites) {
  f <- file.path(PROCESSED_DIR, sprintf("colonization_%s_survival_factorial_v5_v6_h0.40.rds", s))
  d <- readRDS(f)
  t_max <- max(d$t)
  cat(sprintf("--- %s (t_max=%d) ---\n", s, t_max))
  for (p in v5_params) {
    r <- main_effect_range(d, p, t_max)
    cat(sprintf("  %-12s abs_range=%8.1f  rel_range=%.3f  (%d levels)\n", p, r$abs_range, r$rel_range, r$levels))
    v5_summary[[p]] <- c(v5_summary[[p]], r$rel_range)
  }
}
cat("\n=== v5 mean relative-range across sites (ranking) ===\n")
v5_means <- sort(sapply(v5_summary, mean), decreasing = TRUE)
print(round(v5_means, 4))

cat("\n=== Actual tested levels (for building v6's grid) ===\n")
d4 <- readRDS(file.path(PROCESSED_DIR, "colonization_Maquipucuna_reproduction_factorial_v4_v6_h0.40.rds"))
for (p in v4_params) cat(p, ":", paste(sort(unique(d4[[p]])), collapse=", "), "\n")
d5 <- readRDS(file.path(PROCESSED_DIR, "colonization_Maquipucuna_survival_factorial_v5_v6_h0.40.rds"))
for (p in v5_params) cat(p, ":", paste(sort(unique(d5[[p]])), collapse=", "), "\n")
