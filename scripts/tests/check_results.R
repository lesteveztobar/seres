source("scripts/02_model/config/paths.R")

state_avg  <- test$avg$state
state_spat <- test$spatial$state

clim_avg  <- state_avg$clim_by_height          # list, one flat scalar/data.frame per height
clim_spat <- state_spat$clim_spatial_by_height # list, one x,y matrix per height

heights <- state_avg$heights
n_h <- length(heights)

# For each height tier: compare the single flat averaged value against
# EVERY individual voxel value from the spatial version, for temp, relhum,
# and swdown -- the three variables that feed survival and establishment
# (see survival_logit()/run_pass2_establish(), get_colonization.R).
compare_one_height <- function(zi) {
  cl_avg  <- clim_avg[[zi]]
  cl_spat <- clim_spat[[zi]]
  if (is.null(cl_avg) || is.null(cl_spat)) return(NULL)

  avg_temp   <- mean(cl_avg$temp, na.rm = TRUE)
  avg_relhum <- mean(cl_avg$relhum, na.rm = TRUE)
  avg_swdown <- mean(cl_avg$swdown, na.rm = TRUE)

  spat_temp   <- as.numeric(cl_spat$year$temp)
  spat_relhum <- as.numeric(cl_spat$year$relhum)
  spat_swdown <- as.numeric(cl_spat$year$swdown)

  data.frame(
    height            = heights[zi],
    avg_temp          = avg_temp,
    spat_temp_mean    = mean(spat_temp, na.rm = TRUE),
    spat_temp_sd      = sd(spat_temp, na.rm = TRUE),
    spat_temp_range   = diff(range(spat_temp, na.rm = TRUE)),
    avg_relhum        = avg_relhum,
    spat_relhum_mean  = mean(spat_relhum, na.rm = TRUE),
    spat_relhum_sd    = sd(spat_relhum, na.rm = TRUE),
    spat_relhum_range = diff(range(spat_relhum, na.rm = TRUE)),
    avg_swdown        = avg_swdown,
    spat_swdown_mean  = mean(spat_swdown, na.rm = TRUE),
    spat_swdown_sd    = sd(spat_swdown, na.rm = TRUE),
    spat_swdown_range = diff(range(spat_swdown, na.rm = TRUE))
  )
}

comparison <- do.call(rbind, lapply(seq_len(n_h), compare_one_height))

cat("=== Per-height summary (first and last rows) ===\n")
print(head(comparison, 10))
print(tail(comparison, 10))

# QUESTION 1: How much does a "typical" voxel deviate from the single flat
# value lookup_climate_by_height() currently assigns to the whole height tier? A large SD/
# range relative to the biologically meaningful scale of each variable
# (e.g. survival_logit()'s ~1.5x swdown_rel thtesthold, or a few degrees C
# mattering for the heat-sttests term) would mean the flat-average
# simplification is discarding real, decision-relevant heterogeneity.
cat("\n=== How much does the typical voxel deviate from the flat value? ===\n")
cat("Temp   -- mean SD across voxels:", round(mean(comparison$spat_temp_sd, na.rm = TRUE), 3),
    "C | mean range:", round(mean(comparison$spat_temp_range, na.rm = TRUE), 3), "C\n")
cat("RelHum -- mean SD across voxels:", round(mean(comparison$spat_relhum_sd, na.rm = TRUE), 3),
    "% | mean range:", round(mean(comparison$spat_relhum_range, na.rm = TRUE), 3), "%\n")
cat("SWdown -- mean SD across voxels:", round(mean(comparison$spat_swdown_sd, na.rm = TRUE), 3),
    "W/m2 | mean range:", round(mean(comparison$spat_swdown_range, na.rm = TRUE), 3), "W/m2\n")

# QUESTION 2: Is the flat average value itself systematically biased
# relative to the true spatial mean, or does it at least track the spatial
# mean correctly even though it discards the spread? A paired t-test across
# height tiers checks for a consistent offset (e.g. the flat value always
# running warmer/cooler than the real spatial mean), which would matter
# even independently of whether the spread itself is biologically relevant.
cat("\n=== Is there a systematic bias (flat average vs. true spatial mean)? ===\n")
t_temp   <- t.test(comparison$avg_temp,   comparison$spat_temp_mean,   paired = TRUE)
t_relhum <- t.test(comparison$avg_relhum, comparison$spat_relhum_mean, paired = TRUE)
t_swdown <- t.test(comparison$avg_swdown, comparison$spat_swdown_mean, paired = TRUE)
print(t_temp)
print(t_relhum)
print(t_swdown)

saveRDS(comparison, file.path("scripts", "02_model", "tests", "test_testults",
                               "clim_avg_vs_spatial_comparison.rds"))
