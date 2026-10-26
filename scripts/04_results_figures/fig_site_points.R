# fig_site_points.R -- Methods Figure (fig:sitemap): field observation points
# at all seven sites, Mindo-cluster inset, Yanayacu labelled as the
# eastern-slope outgroup. Wraps plot_site_map() (plot_functions.R) and checks
# that every site, Saloya included, made it onto the map.
#   sbatch scripts/04_results_figures/run_fig.sh fig_site_points.R
source("scripts/02_model/config/paths.R")
suppressMessages(source("scripts/02_model/plots/plot_functions.R"))
p <- plot_site_map(out_dir = file.path(OUTPUT_DIR, "figures", "sites"))
cat("plot_site_map() checks internally that all seven sites, Saloya included, have points.\n")
