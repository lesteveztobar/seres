# build_lhs_design.R -- generic Latin hypercube design builder for the
# sensitivity run harness (sensitivity_run.R).
#
# Usage: Rscript build_lhs_design.R <tag> <sites> <n_points> <n_reps> <ranges_csv> [seed]
#   tag         label for this design (used in file names and run ids)
#   sites       comma-separated site names; every site gets the SAME parameter
#               sets and seeds, so results are comparable across sites
#   n_points    number of parameter sets
#   n_reps      replicates per parameter set (seeds 1..n_reps)
#   ranges_csv  columns: param, min, max, integer (TRUE = round to integer)
#   seed        RNG seed for the hypercube (default 1)
# Writes data/params/sensitivity_design_<tag>_<site>.rds, one per site.
suppressPackageStartupMessages({
  .libPaths(c("~/R/marvin_libs", .libPaths()))
  library(lhs)
})
source("scripts/02_model/config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: build_lhs_design.R <tag> <sites> <n_points> <n_reps> <ranges_csv> [seed]")
tag <- args[1]
sites <- strsplit(args[2], ",")[[1]]
n_points <- as.integer(args[3])
n_reps <- as.integer(args[4])
ranges <- read.csv(args[5], stringsAsFactors = FALSE)
seed <- if (length(args) >= 6) as.integer(args[6]) else 1L
stopifnot(all(c("param", "min", "max", "integer") %in% names(ranges)), all(ranges$max >= ranges$min))

set.seed(seed)
u <- randomLHS(n_points, nrow(ranges))
points <- as.data.frame(lapply(seq_len(nrow(ranges)), function(j) {
  v <- ranges$min[j] + u[, j] * (ranges$max[j] - ranges$min[j])
  if (isTRUE(as.logical(ranges$integer[j]))) as.integer(round(v)) else v
}))
names(points) <- ranges$param

for (site in sites) {
  p <- points
  p$point_id <- sprintf("%s_%04d", site, seq_len(n_points))
  p$site <- site
  design <- do.call(rbind, lapply(seq_len(n_reps), function(r) { x <- p; x$rep <- r; x$seed <- r; x }))
  design <- design[order(design$point_id, design$rep), ]
  design$run_id <- sprintf("%s_%s_r%d", tag, design$point_id, design$rep)
  design <- design[, c("run_id", "point_id", "site", "rep", "seed", ranges$param)]
  out <- file.path(PARAMS_DIR, sprintf("sensitivity_design_%s_%s.rds", tag, site))
  saveRDS(design, out)
  cat(sprintf("%-14s %d runs (%d parameter sets x %d seeds) -> %s\n", site, nrow(design), n_points, n_reps, out))
}
for (j in seq_len(nrow(ranges))) cat(sprintf("  %-11s [%.6g, %.6g]\n", ranges$param[j], min(points[[j]]), max(points[[j]])))
