source("scripts/02_model/config/patches.R")
reticulate::use_python(Sys.getenv("CANOPY_PYTHON",
  unset = "/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python"), required = TRUE)
library(microclimf)
cat("=== runpointmodela ===\n")
print(microclimf::runpointmodela)
cat("\n=== .pointmodel (if it exists, the likely internal worker) ===\n")
tryCatch(print(microclimf:::.pointmodel), error = function(e) cat("not found:", conditionMessage(e), "\n"))
