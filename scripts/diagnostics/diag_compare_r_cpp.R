# diag_compare_r_cpp.R — Step 2 correctness check: compares runmicro()'s
# method="R" vs method="Cpp" full-field output for the same site/height,
# produced by diag_wrap_method_test.R (jobs 26906059 Cpp, 26906060 R).
# Not part of the pipeline; submitted manually, read-only against the two
# .rds files -- writes nothing back except a small summary.
setwd("/home/s38leste_hpc/seres")

cat(sprintf("[%s] Loading R-path output...\n", format(Sys.time(), "%H:%M:%S")))
r_out <- readRDS("data/processed/diag_wrap_heights/r_serial_h5.10.rds")
cat(sprintf("[%s] Loading Cpp-path output...\n", format(Sys.time(), "%H:%M:%S")))
cpp_out <- readRDS("data/processed/diag_wrap_heights/cpp_serial_h5.10.rds")

cat("\n=== structure ===\n")
cat("R Tz dim:  ", paste(dim(r_out$Tz), collapse = " x "), "\n")
cat("Cpp Tz dim:", paste(dim(cpp_out$Tz), collapse = " x "), "\n")
cat("R vars:  ", paste(names(r_out), collapse = ", "), "\n")
cat("Cpp vars:", paste(names(cpp_out), collapse = ", "), "\n")
cat("R n timesteps:", r_out$n, " Cpp n timesteps:", cpp_out$n, "\n")
cat("tme identical:", isTRUE(all.equal(r_out$tme, cpp_out$tme)), "\n")

set.seed(1)
nr <- dim(r_out$Tz)[1]; nc <- dim(r_out$Tz)[2]; nt <- dim(r_out$Tz)[3]
idx <- data.frame(row = sample(nr, 8, TRUE), col = sample(nc, 8, TRUE), t = sample(nt, 8, TRUE))

vars <- c("Tz", "relhum", "windspeed", "Rdirdown", "Rdifdown")

cat("\n=== sample pixel/hour comparisons ===\n")
for (v in vars) {
  cat(sprintf("--- %s ---\n", v))
  for (k in seq_len(nrow(idx))) {
    i <- idx$row[k]; j <- idx$col[k]; t <- idx$t[k]
    rv <- r_out[[v]][i, j, t]; cv <- cpp_out[[v]][i, j, t]
    cat(sprintf("  [%d,%d,t=%d]  R=%.6f  Cpp=%.6f  diff=%.6g\n", i, j, t, rv, cv, rv - cv))
  }
}

cat("\n=== overall summary diffs (full arrays, one var at a time to limit peak memory) ===\n")
for (v in vars) {
  rv <- r_out[[v]]; cv <- cpp_out[[v]]
  n_na_r <- sum(is.na(rv)); n_na_cpp <- sum(is.na(cv))
  d <- as.vector(rv) - as.vector(cv)
  cat(sprintf("%s: max|diff|=%.6g  mean|diff|=%.6g  n_NA_R=%d  n_NA_Cpp=%d  n_total=%d\n",
              v, max(abs(d), na.rm = TRUE), mean(abs(d), na.rm = TRUE), n_na_r, n_na_cpp, length(d)))
  rm(rv, cv, d); gc(FALSE)
}
cat(sprintf("[%s] Done.\n", format(Sys.time(), "%H:%M:%S")))
