# compare_niche_caches.R -- old (pre-rebuild) vs rebuilt niche cache. Raw
# per-site climate values are not stored in the cache (only the derived
# density-ratio score curves), and the pre-rebuild per-site manifests were
# overwritten in place, so the comparison is on what downstream consumes:
# per-species axis score curves and the pooled background density.
args <- commandArgs(TRUE); old_dir <- args[1]
source("scripts/02_model/config/paths.R")
on <- readRDS(file.path(old_dir, "species_niches.rds")); nn <- readRDS(NICHE_CACHE_PATH)
ob <- readRDS(file.path(old_dir, "niche_background_density.rds")); nb <- readRDS(NICHE_BACKGROUND_PATH)
cat(sprintf("species: old %d, new %d, identical names: %s\n", length(on), length(nn), identical(sort(names(on)), sort(names(nn)))))
d <- do.call(rbind, lapply(intersect(names(on), names(nn)), function(sp) do.call(rbind, lapply(names(on[[sp]]$axes), function(ax) {
  a <- on[[sp]]$axes[[ax]]$score; b <- nn[[sp]]$axes[[ax]]$score
  data.frame(species = sp, axis = ax, max_abs_diff = max(abs(a - b), na.rm = TRUE), mean_abs_diff = mean(abs(a - b), na.rm = TRUE),
             argmax_shift = nn[[sp]]$axes[[ax]]$x[which.max(b)] - on[[sp]]$axes[[ax]]$x[which.max(a)])
}))))
cat("\nper-axis score-curve change (0-100 scale), summary over species:\n"); print(aggregate(cbind(max_abs_diff, mean_abs_diff) ~ axis, d, function(x) round(c(median = median(x), max = max(x)), 2)))
cat("\nspecies with largest change:\n"); print(head(d[order(-d$max_abs_diff), ], 8), row.names = FALSE)
cat("\npooled background density, quantile of x weighted by y (old vs new):\n")
for (v in names(ob)) { q <- function(b) { o <- b[[v]]; cs <- cumsum(o$y) / sum(o$y); approx(cs, o$x, c(.05, .25, .5, .75, .95), rule = 2, ties = "ordered")$y }
  cat(sprintf("  %-8s old %s | new %s\n", v, paste(round(q(ob), 2), collapse = " "), paste(round(q(nb), 2), collapse = " "))) }
out_path <- resolve_output_path("niche_cache_compare", "niche_cache_old_vs_new.csv",
                                 legacy_path = file.path(OUTPUT_DIR, "niche_cache_old_vs_new.csv"))
write.csv(d, out_path, row.names = FALSE)
