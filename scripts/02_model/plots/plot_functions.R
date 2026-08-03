# plot_functions.R
# Reusable plotting functions for the whole canopymicroenv pipeline. Each
# function loads its own inputs from data/processed (or geojson_to_csv) and
# saves its output(s) into OUTPUT_DIR — nothing here plots-and-forgets.
# If an input file doesn't exist yet, the function skips with a message
# instead of erroring, so plot_all.R can be re-run at any point in the
# pipeline and it just draws whatever is available so far.
#
# Refactored from the standalone plotmap.R / plot_temperatures.R /
# plot_bestfit_3d.R, plus new export wrappers around the plot_abundance() /
# plot_3d_abundance() / plot_experiment() functions already defined in
# get_colonization.R.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2)
library(patchwork)
library(plotly)
library(htmlwidgets)
library(abind)
library(scatterplot3d)

source("scripts/02_model/config/paths.R")
source("scripts/02_model/engine/get_colonization.R")

# ── Field site map ─────────────────────────────────────────────────────────────
# NW Ecuador Maxillariinae field sites, from GeoJSON transects/points.
# sf/dplyr/rnaturalearth(data)/ggspatial/scico are loaded here rather than at
# file level: this is the only function in the file that needs them (all are
# geospatial/mapping-specific, several with heavy system-library dependencies
# -- GDAL/PROJ/GEOS/UDUNITS -- that every other function in this file has no
# reason to require just to be sourced).
plot_site_map <- function(geojson_dir = file.path(BASE_DIR, "geojson_to_csv", "raw"),
                           out_dir = OUTPUT_DIR) {
  library(sf)
  library(dplyr)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggspatial)
  library(scico)

  # MiradorMindo.geojson retired 2026-07-22: it actually bundled observations
  # from 3 separate locations (LaElenita, MindoMirador, Saloya), now split
  # into their own GeoJSON exports (see scripts/00_data_conversion/rebuild_combined_csv.py).
  # Left on disk as an audit trail but no longer read here.
  geojson_files <- c(
    Maquipucuna   = "Maquipucuna.geojson",
    Mashpi        = "Mashpi.geojson",
    MindoTarabita = "MindoTarabita.geojson",   # merged with TarabitaMindo below
    TarabitaMindo = "TarabitaMindo.geojson",
    LaElenita     = "LaElenita.geojson",
    MindoMirador  = "MindoMirador.geojson",
    Saloya        = "Saloya.geojson",
    Yanayacu      = "Yanayacu.geojson"
  )

  load_geojson <- function(filename, site_name, dir = geojson_dir) {
    sf::st_read(file.path(dir, filename), quiet = TRUE) %>%
      transmute(site = site_name,
                description = if ("description" %in% names(.)) description else NA_character_)
  }

  site_list    <- mapply(load_geojson, geojson_files, names(geojson_files), SIMPLIFY = FALSE)
  all_features <- do.call(rbind, site_list)
  all_features <- all_features %>%
    mutate(site = if_else(site == "TarabitaMindo", "MindoTarabita", site))

  obs_points <- all_features %>% filter(sf::st_geometry_type(geometry) == "POINT")
  transects  <- all_features %>% filter(sf::st_geometry_type(geometry) == "LINESTRING")

  transect_labels <- transects %>%
    group_by(site) %>%
    summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    mutate(geometry = sf::st_centroid(geometry))

  fallback_labels <- obs_points %>%
    filter(!site %in% transect_labels$site) %>%
    group_by(site) %>%
    summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    mutate(geometry = sf::st_centroid(geometry))

  site_labels <- rbind(transect_labels, fallback_labels)

  # Per-site label offsets (degrees) — keeps long names off the dots.
  # LaElenita/MindoMirador/Saloya added 2026-07-24 with nudge (0,0)
  # placeholders (untuned) -- adjust once you've seen how their labels sit
  # on the actual map.
  label_nudges <- data.frame(
    site    = c("Mashpi", "Maquipucuna", "MindoTarabita", "Yanayacu",
                "LaElenita", "MindoMirador", "Saloya"),
    nudge_x = c( 0.00,     0.08,          0.15,            0.10,
                 0.00,      0.00,          0.00),
    nudge_y = c(-0.04,     0.07,         -0.05,            0.06,
                 0.00,      0.00,          0.00)
  )
  label_pos <- site_labels |>
    dplyr::mutate(X = sf::st_coordinates(geometry)[, 1],
                  Y = sf::st_coordinates(geometry)[, 2]) |>
    sf::st_drop_geometry() |>
    dplyr::left_join(label_nudges, by = "site") |>
    dplyr::mutate(X = X + nudge_x, Y = Y + nudge_y)

  ecuador  <- ne_states(country = "Ecuador", returnclass = "sf")
  nw_provs <- c("Pichincha", "Esmeraldas", "Imbabura",
                "Santo Domingo de los Tsáchilas", "Cotopaxi",
                "Napo", "Sucumbios", "Tungurahua")
  map_area <- ecuador %>% filter(name %in% nw_provs)

  bbox <- sf::st_bbox(obs_points)
  xlim <- c(bbox["xmin"] - 0.30, bbox["xmax"] + 0.55)  # extra east for Yanayacu
  ylim <- c(bbox["ymin"] - 0.40, bbox["ymax"] + 0.30)

  site_colours <- setNames(
    scico::scico(7, palette = "lipari", begin = 0.10, end = 0.88),
    c("Maquipucuna", "Mashpi", "MindoTarabita", "LaElenita", "MindoMirador", "Saloya", "Yanayacu")
  )

  p <- ggplot() +
    geom_sf(data = map_area, fill = "#f5f0e8", colour = "grey55", linewidth = 0.35) +
    geom_sf(data = obs_points, aes(colour = site), size = 1.8, alpha = 0.55, shape = 16) +
    geom_sf(data = transects, aes(colour = site), linewidth = 1.3, alpha = 0.85) +
    geom_text(data = label_pos, aes(x = X, y = Y, label = site, colour = site),
              size = 3, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = site_colours, name = "Site") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tr", style = north_arrow_fancy_orienteering(),
      height = unit(1.2, "cm"), width = unit(1.2, "cm")
    ) +
    labs(title = "Field sites", subtitle = "NW Ecuador · Chocó Andino ",
         x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey40")
    )

  out_path <- file.path(out_dir, "site_map.png")
  ggsave(out_path, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Modelled tree diagram ────────────────────────────────────────────────────
# Static illustrative figure for report/methods.tex (Section "Landscape and
# forest structure"): renders ONE tree exactly as build_forest()
# (get_colonization.R) classifies it -- Johansson zones 1-5, trunk restricted
# to a single column, crown following the bell-shaped sin() profile peaking
# at 65% of tree height. Not a simulation output -- draws one tree at its
# forestparams mean height/crown radius (deterministic, not a stochastic
# draw), at fine continuous resolution, purely to visualize the geometry the
# methods text describes. Defaults match the forestparams list used
# everywhere else (run_colonization.R etc.); pass a different list to
# match a different site/config if that ever stops being shared across sites
# (see methods.tex, Landscape and forest structure, on why it's shared now).
plot_tree_diagram <- function(forestparams = list(mean_hgt = 8.4, mean_crown_r = 2.0,
                                                    trunk_r = 0.114),
                               out_dir = OUTPUT_DIR, res = 0.25) {
  th <- forestparams$mean_hgt
  cr <- forestparams$mean_crown_r

  heights      <- seq(0.5, th, by = res)
  half_extent  <- ceiling(cr / res) * res
  xy           <- seq(-half_extent, half_extent, by = res)

  rows <- list()
  for (h in heights) {
    rel_h <- h / th
    # Same thresholds as build_forest() (get_colonization.R) -- keep in sync.
    jzone <- if      (rel_h < 0.10) 1L
             else if (rel_h < 0.30) 2L
             else if (rel_h < 0.50) 3L
             else if (rel_h < 0.80) 4L
             else                   5L
    if (jzone <= 2) {
      rows[[length(rows) + 1]] <- data.frame(x = 0, y = 0, z = h, zone = jzone)
    } else {
      crown_fraction <- (rel_h - 0.30) / 0.70
      effective_r    <- cr * sin(crown_fraction * pi)
      if (effective_r <= 0) next
      grid      <- expand.grid(x = xy, y = xy)
      grid$dist <- sqrt(grid$x^2 + grid$y^2)
      grid      <- grid[grid$dist <= effective_r, ]
      if (nrow(grid) == 0) next
      rows[[length(rows) + 1]] <- data.frame(x = grid$x, y = grid$y, z = h, zone = jzone)
    }
  }
  pts <- do.call(rbind, rows)

  zone_cols   <- c(`1` = "#6b4423", `2` = "#8a6d3b", `3` = "#a8d18d",
                    `4` = "#4f9a4f", `5` = "#2d6a2d")
  zone_labels <- c(`1` = "Zone 1 (trunk base)",  `2` = "Zone 2 (trunk)",
                    `3` = "Zone 3 (lower crown)", `4` = "Zone 4 (mid crown)",
                    `5` = "Zone 5 (upper crown)")
  present <- as.character(sort(unique(pts$zone)))

  out_path <- file.path(out_dir, "tree_diagram.png")
  # width/height/mar tuned so neither the title nor the legend run past the
  # canvas edge (base R does not wrap or auto-shrink `main=` text -- a title
  # wider than the device just gets clipped symmetrically at both edges,
  # which is what happened before: the single-line title was ~8in wide
  # rendered into a 6.67in-wide canvas). Splitting the title onto two lines
  # and keeping the legend INSIDE the plot box (rather than pushed out into
  # the margin with a negative inset) fixes both at once.
  png(out_path, width = 2200, height = 1800, res = 300, type = "cairo", bg = "white")
  on.exit(dev.off(), add = TRUE)
  par(mar = c(2, 2, 4, 1))
  scatterplot3d::scatterplot3d(
    pts$x, pts$y, pts$z,
    color = zone_cols[as.character(pts$zone)],
    pch = 16, cex.symbols = 0.5,
    xlab = "x (m)", ylab = "y (m)", zlab = "Height (m)",
    main = sprintf("Modelled tree geometry\n(height = %.1f m, crown radius = %.1f m)", th, cr),
    cex.main = 1.1,
    angle = 55, scale.y = 0.7, grid = TRUE, box = FALSE,
    col.axis = "grey40", col.grid = "grey88", col.lab = "grey20",
    cex.axis = 0.8, cex.lab = 0.9
  )
  legend("topleft", inset = 0.02, legend = zone_labels[present], col = zone_cols[present],
         pch = 16, cex = 0.7, bty = "n")
  message("Saved: ", out_path)
  invisible(pts)
}

# ── Microclimate temperature profile ────────────────────────────────────────────
# 3D temperature volume (interactive, plotly) plus a side-view heatmap +
# boxplot with Kruskal-Wallis / Dunn post-hoc height comparison, for one
# site's microenv RDS.

plot_temperature_profile <- function(site_name, out_dir = OUTPUT_DIR,
                                      processed_dir = PROCESSED_DIR,
                                      height_step = 0.25) {
  # height_step defaults to 0.25 (the production resolution used everywhere
  # else in this pipeline -- run_colonization.R, characterize_niches.R
  # -- via the same manifest_suffix convention), NOT the unsuffixed 0.1m file.
  # 2026-07-17: this function used to hardcode the unsuffixed name, which only
  # Maquipucuna happens to have (and at ~2.5x more height tiers than the
  # 0.25m file, ~168 vs ~67) -- every other site was silently skipped with
  # "no such file" on every run_plots.sh pass.
  manifest_suffix <- if (height_step != 0.1) sprintf("_h%.2f", height_step) else ""
  microenv_path <- file.path(processed_dir, sprintf("microenv_%s%s.rds", site_name, manifest_suffix))
  if (!file.exists(microenv_path)) {
    message("Skipping ", site_name, " temperature profile -- no ", microenv_path)
    return(invisible(NULL))
  }
  env       <- readRDS(microenv_path)
  heights_m <- microenv_heights(env)
  # load_height()/microenv_heights() (get_colonization.R) handle both the new
  # manifest-only format (.heights + .height_dir, per-height RDS on scratch)
  # and the old format where each height is embedded as its own list element.
  #
  # Stream one height at a time instead of loading every height's full raw
  # per-timestep spatial array into memory at once -- each height's raw
  # tmax/tmin object is ~6 GB (48 hourly rasters), and with ~170 heights
  # that's >1 TB, reliably OOM-killing the job regardless of how much memory
  # is requested. Everything downstream (volume plot, cross-section heatmap,
  # boxplot) only ever needs the time-averaged 2D Tz layer per height, so
  # reduce immediately on load and keep only that (tiny) result.
  tz_layers <- lapply(heights_m, function(hgt) {
    h <- load_height(env, hgt)
    # New format (h$tme present, full-year rewrite): Tz sits at the top
    # level, no tmax/tmin day-type split. Old format: nested under $tmax.
    r <- if (!is.null(h$tme)) h$Tz else h$tmax$Tz
    if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE) else r
  })
  vol_tmax <- abind::abind(tz_layers, along = 3)

  nrow_r  <- dim(vol_tmax)[1]
  ncol_r  <- dim(vol_tmax)[2]
  nheight <- dim(vol_tmax)[3]

  if (!is.null(env$.spatial)) {
    # .spatial$ext is a plain named numeric vector (xmin/xmax/ymin/ymax), not
    # a terra::ext() S4 object -- see run_microclimate_site.R for why.
    e        <- env$.spatial$ext
    x_coords <- seq(e[["xmin"]], e[["xmax"]], length.out = ncol_r)
    y_coords <- seq(e[["ymax"]], e[["ymin"]], length.out = nrow_r)  # north→south
    x_lab <- "Longitude"; y_lab <- "Latitude"
  } else {
    x_coords <- seq_len(ncol_r)
    y_coords <- seq_len(nrow_r)
    x_lab <- "col (W→E)"; y_lab <- "row (S→N)"
  }

  # ── 3D volume (plotly) ────────────────────────────────────────────────────
  grid <- expand.grid(xi = x_coords, yi = y_coords, zi = seq_len(nheight))
  grid$z  <- heights_m[grid$zi]
  grid$Tz <- vol_tmax[cbind(nrow_r + 1 - grid$yi, grid$xi, grid$zi)]
  grid    <- grid[!is.na(grid$Tz), ]

  fig <- plot_ly(
    data = grid, x = ~xi, y = ~yi, z = ~z, value = ~Tz, type = "volume",
    isomin = quantile(grid$Tz, 0.05), isomax = quantile(grid$Tz, 0.95),
    opacity = 0.15, surface = list(count = 10),
    colorscale = list(
      c(0,    "#313695"), c(0.25, "#74add1"), c(0.5, "#ffffbf"),
      c(0.75, "#f46d43"), c(1,    "#a50026")
    ),
    colorbar = list(title = "T (°C)")
  ) |>
    layout(
      title = sprintf("Air temperature 3D volume — %s (warmest representative day)", site_name),
      scene = list(
        xaxis = list(title = x_lab), yaxis = list(title = y_lab),
        zaxis = list(title = "Height (m)")
      )
    )

  volume_path <- file.path(out_dir, sprintf("temp_volume_%s.html", site_name))
  # selfcontained=TRUE needs pandoc (not installed on the cluster); FALSE
  # writes a small "<name>_files/" dependency folder alongside the HTML
  # instead -- keep the two together when copying/viewing elsewhere.
  htmlwidgets::saveWidget(fig, volume_path, selfcontained = FALSE)
  message("Saved: ", volume_path)

  # ── Side-view heatmap + per-height boxplot ────────────────────────────────
  side_df <- do.call(rbind, lapply(seq_along(tz_layers), function(i) {
    col_means <- colMeans(tz_layers[[i]], na.rm = TRUE)
    data.frame(col = seq_along(col_means), height = heights_m[i], Tz = col_means)
  }))
  side_df <- side_df[!is.na(side_df$Tz), ]

  # All pixel values per height — more power for the test than col means alone
  all_px <- do.call(rbind, lapply(seq_along(tz_layers), function(i) {
    data.frame(height = factor(heights_m[i]), Tz = as.vector(tz_layers[[i]]))
  }))
  all_px <- all_px[!is.na(all_px$Tz), ]

  # Kruskal-Wallis: are any height levels different?
  kw <- kruskal.test(Tz ~ height, data = all_px)
  kw_label <- sprintf("Kruskal-Wallis: χ²(%.0f) = %.2f, p %s",
                      kw$parameter, kw$statistic,
                      ifelse(kw$p.value < 0.001, "< 0.001",
                             sprintf("= %.3f", kw$p.value)))
  message(kw_label)

  # Dunn post-hoc pairwise comparisons (Holm correction)
  if (!requireNamespace("dunn.test", quietly = TRUE)) install.packages("dunn.test")
  dunn_res <- dunn.test::dunn.test(all_px$Tz, all_px$height,
                                    method = "holm", kw = FALSE, label = FALSE)
  dunn_df <- data.frame(
    comparison = dunn_res$comparisons,
    p_adj      = dunn_res$P.adjusted,
    sig        = ifelse(dunn_res$P.adjusted < 0.001, "***",
                 ifelse(dunn_res$P.adjusted < 0.01,  "**",
                 ifelse(dunn_res$P.adjusted < 0.05,  "*", "ns")))
  )
  message("Pairwise Dunn tests (Holm-adjusted):")
  print(dunn_df[order(dunn_df$p_adj), ], row.names = FALSE)

  # Compact letter display for the boxplot annotation
  if (!requireNamespace("multcompView", quietly = TRUE)) install.packages("multcompView")
  p_mat  <- setNames(dunn_res$P.adjusted, dunn_res$comparisons)
  cld    <- multcompView::multcompLetters(p_mat, threshold = 0.05)$Letters
  cld_df <- data.frame(height = as.numeric(names(cld)), letter = cld)

  p_heat <- ggplot(side_df, aes(x = col, y = height, fill = Tz)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"),
      name = "T (°C)", na.value = "grey90"
    ) +
    labs(x = "Horizontal position (col, W→E)", y = "Height (m)",
         title = sprintf("Vertical cross-section — %s (warmest representative day)", site_name),
         subtitle = kw_label) +
    theme_minimal(base_size = 11)

  p_box <- ggplot(all_px, aes(x = Tz, y = height, group = height)) +
    geom_boxplot(aes(fill = after_stat(middle)), width = 0.06, outlier.size = 0.5) +
    geom_text(data = cld_df,
              aes(x = max(all_px$Tz, na.rm = TRUE), y = height, label = letter),
              hjust = -0.2, size = 3.5, inherit.aes = FALSE) +
    scale_fill_gradientn(
      colours = c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"), guide = "none"
    ) +
    labs(x = "T (°C)", y = NULL,
         caption = "Letters = Dunn post-hoc (Holm); shared letter → no significant difference") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_blank())

  p_combined  <- p_heat + p_box + plot_layout(widths = c(3, 1))
  profile_path <- file.path(out_dir, sprintf("temp_profile_%s.png", site_name))
  ggsave(profile_path, plot = p_combined, width = 11, height = 5, dpi = 300, bg = "white")
  message("Saved: ", profile_path)

  invisible(list(volume = fig, profile = p_combined))
}

# ── Best-fit vs. realistic 3D comparison (simple standalone model) ─────────────
plot_bestfit_3d_comparison <- function(out_dir = OUTPUT_DIR) {
  if (!exists("best_run", inherits = TRUE)) {
    source("scripts/simple_model/simple_colonization.R")
    best_params <- params
    best_params$establishment_prob <- 0.10
    best_params$repro_rate         <- 500
    best_params$survival_A         <- 0.95
    best_params$maxDisp            <- 4
    best_params$n_founders         <- 100
    best_params$carCap             <- 40
    best_params$timesteps          <- 30
    message("Re-running best fit...")
    best_run <- run_simple_colonization(best_params, seed = 42)
  }

  p      <- best_run$params
  T_last <- p$timesteps

  abundA_last <- best_run$abundA[, , , T_last]
  zone_arr    <- best_run$zone

  idx <- which(abundA_last > 0, arr.ind = TRUE)
  df  <- data.frame(
    x = idx[, 1], y = idx[, 2], z = idx[, 3],
    n = abundA_last[idx], zone = zone_arr[idx]
  )
  df$x_m <- (df$x - 0.5) * p$resolution
  df$y_m <- (df$y - 0.5) * p$resolution
  df$z_m <- (df$z - 0.5) * (p$max_height / p$zDim)

  zone_cols <- c(
    "1" = "#8B4513", "2" = "#A0522D", "3" = "#6B8E23",
    "4" = "#228B22", "5" = "#32CD32"
  )
  df$col     <- zone_cols[as.character(df$zone)]
  df$pt_size <- pmin(df$n / max(df$n) * 6 + 1, 7)

  # Realistic-params run for comparison (above extinction threshold, not saturating)
  source("scripts/simple_model/simple_colonization.R")
  base_p <- params
  base_p$establishment_prob <- 0.05
  base_p$repro_rate         <- 150
  base_p$carCap             <- 8
  base_p$n_founders         <- 30
  base_p$timesteps          <- 30
  message("Running baseline simulation...")
  invisible(capture.output(
    base_out <- run_simple_colonization(base_p, seed = 7),
    type = "output"
  ))

  make_df <- function(run_out, timestep) {
    abA  <- run_out$abundA[, , , timestep]
    zArr <- run_out$zone
    idx  <- which(abA > 0, arr.ind = TRUE)
    if (nrow(idx) == 0) return(NULL)
    data.frame(
      x_m = (idx[, 1] - 0.5) * run_out$params$resolution,
      y_m = (idx[, 2] - 0.5) * run_out$params$resolution,
      z_m = (idx[, 3] - 0.5) * (run_out$params$max_height / run_out$params$zDim),
      n = abA[idx], zone = zArr[idx]
    )
  }

  df_base <- make_df(base_out, base_p$timesteps)
  df_best <- make_df(best_run, T_last)

  set.seed(42)
  df_best_plot <- df_best[sample(nrow(df_best), min(nrow(df_best), 5000)), ]
  df_base_plot <- df_base   # usually sparse enough to keep all

  add_cols <- function(df) { df$col <- zone_cols[as.character(df$zone)]; df }
  df_base_plot <- add_cols(df_base_plot)
  df_best_plot <- add_cols(df_best_plot)

  out_path <- file.path(out_dir, "bestfit_3d.png")
  png(out_path, width = 3200, height = 1400, res = 180, type = "cairo")
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(10, 10, 3))

  plot_panel <- function(df_p, title_str, cex_sym = 0.55) {
    par(mar = c(2, 2, 3, 1))
    scatterplot3d(
      x = df_p$x_m, y = df_p$y_m, z = df_p$z_m,
      color = df_p$col, pch = 16, cex.symbols = cex_sym,
      xlab = "East–West (m)", ylab = "South–North (m)", zlab = "Height (m)",
      main = title_str, angle = 35, scale.y = 0.6, grid = TRUE, box = FALSE,
      col.axis = "grey40", col.grid = "grey88", col.lab = "grey20",
      cex.axis = 0.8, cex.lab = 0.9
    )
  }

  plot_panel(df_base_plot,
    sprintf("Realistic params  (p_e=0.05, λ=150, carCap=8)\nyear %d — %d adults",
            base_p$timesteps, sum(df_base$n)))

  plot_panel(df_best_plot,
    sprintf("Best-fit params  (p_e=0.10, λ=500, carCap=40)\nyear %d — %d adults  [5k sample]",
            T_last, sum(df_best$n)))

  par(mar = c(2, 0, 3, 1))
  plot.new()
  legend(
    "center",
    legend = c("Zone 1 – trunk base", "Zone 2 – lower trunk",
               "Zone 3 – upper trunk", "Zone 4 – inner crown",
               "Zone 5 – outer crown"),
    col = unname(zone_cols), pch = 16, pt.cex = 1.5, bty = "n", cex = 0.95,
    title = "Johansson zone", title.col = "grey20", title.font = 2
  )

  message("Saved: ", out_path)
  invisible(out_path)
}

# ── Colonization sensitivity-experiment sweeps ──────────────────────────────────
# One RDS per site x experiment tag, produced by run_colonization.R
# via batch_exp.sh. This mapping must stay in sync with the EXP/PARAMS
# pairing in scripts/batch_exp.sh.
EXP_PARAM_MAP <- c(
  pollination_success       = "p_poll",
  adult_survival_intercept  = "beta0A",
  germination_probability   = "p_germ",
  reproduction_cost         = "cost_repro",
  climate_sensitivity_rh    = "beta_rh",
  precipitation_sensitivity = "beta_precip",
  founder_number            = "n_founders"
)

plot_colonization_experiment <- function(site_name, exp_tag, out_dir = OUTPUT_DIR,
                                          processed_dir = PROCESSED_DIR) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (!is.data.frame(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- single-run result, use plot_default_colonization_run() instead")
    return(invisible(NULL))
  }
  if (!"param_value" %in% names(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- factorial result (no param_value column), use plot_factorial_experiment() instead")
    return(invisible(NULL))
  }

  param_name <- EXP_PARAM_MAP[[exp_tag]] %||% exp_tag
  p <- plot_experiment(result, param_name,
                       title = sprintf("Effect of %s — %s", param_name, site_name))

  out_path <- file.path(out_dir, sprintf("experiment_%s_%s.png", site_name, exp_tag))
  ggsave(out_path, plot = p, width = 12, height = 5, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Factorial experiment (e.g. n_founders x p_poll x p_germ x p_s1) ────────────
# Extinction-rate heatmap across two swept parameters, faceted by however many
# more are present — generalizes to any factorial from run_factorial_experiment(),
# not just a specific 3- or 4-parameter design. Swept columns are detected as
# whatever's left after removing the fixed output columns, with n_founders
# (if present) placed on the primary x-axis since it's usually the parameter
# of most direct interest. Mirrors the simple model's Fig. 5 (report.pdf
# sec. 3.3), generalized past a single facet variable via facet_grid.
plot_factorial_experiment <- function(site_name, exp_tag = "reproduction_factorial",
                                      out_dir = OUTPUT_DIR, processed_dir = PROCESSED_DIR) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (!is.data.frame(result) || "param_value" %in% names(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- not a factorial result, use plot_colonization_experiment() instead")
    return(invisible(NULL))
  }
  fixed_cols <- c("rep", "t", "totalS", "totalJ", "totalA", "total", "extinct")
  swept <- setdiff(names(result), fixed_cols)
  if ("n_founders" %in% swept) swept <- c("n_founders", setdiff(swept, "n_founders"))
  if (length(swept) < 2) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- fewer than 2 swept columns found (", paste(swept, collapse = ", "), ")")
    return(invisible(NULL))
  }

  t_max <- max(result$t)
  final <- result[result$t == t_max, ]
  form  <- as.formula(paste("cbind(extinct, total) ~", paste(swept, collapse = " + ")))
  combo_summary <- aggregate(form, data = final, FUN = mean)

  x_var <- swept[1]; y_var <- swept[2]
  facet_vars <- swept[-(1:2)]

  p <- ggplot(combo_summary,
             aes(x = factor(.data[[x_var]]), y = factor(.data[[y_var]]), fill = extinct)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("#08519c", "#f1a340", "#a50026"),
                         limits = c(0, 1), name = "Extinction\nrate") +
    labs(x = x_var, y = y_var,
         title = sprintf("Factorial — %s", site_name),
         subtitle = if (length(facet_vars) > 0)
           sprintf("Extinction rate at year %d, faceted by %s",
                   t_max, paste(facet_vars, collapse = " x "))
         else
           sprintf("Extinction rate at year %d", t_max)) +
    theme_minimal(base_size = 11)

  if (length(facet_vars) == 1) {
    p <- p + facet_wrap(as.formula(paste("~", facet_vars[1])), labeller = label_both)
  } else if (length(facet_vars) >= 2) {
    p <- p + facet_grid(as.formula(paste(facet_vars[1], "~", facet_vars[2])), labeller = label_both)
    if (length(facet_vars) > 2)
      message("Note: faceting by ", facet_vars[1], " and ", facet_vars[2],
              " only; ", paste(facet_vars[-(1:2)], collapse = ", "),
              " collapsed via aggregation.")
  }

  out_path <- file.path(out_dir, sprintf("factorial_%s_%s.png", site_name, exp_tag))
  ggsave(out_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Shared setup for the niche-suitability plotting functions below ────────
# Building the climate cache is the expensive part (one Lustre read per
# height tier) -- computing it once and sharing it between
# plot_niche_suitability() and plot_niche_profile_curves() (as plot_all.R
# does, both called back to back for the same site) avoids paying that cost
# twice per site. height_step defaults to 0.25, the production resolution
# used everywhere else in this pipeline, NOT the unsuffixed 0.1m file --
# 2026-07-17: the unsuffixed default previously meant every niche-suitability
# plot either skipped with "no such file" (4/5 sites only have the _h0.25
# microenv) or, for Maquipucuna (which has both), rebuilt an unnecessarily
# large 168-tier cache instead of the 67-tier 0.25m one, slow enough to blow
# through run_plots.sh's 2-hour budget on its own.
.niche_plot_context <- function(site_name, processed_dir = PROCESSED_DIR, height_step = 0.25) {
  manifest_suffix  <- if (height_step != 0.1) sprintf("_h%.2f", height_step) else ""
  niche_cache_path <- file.path(processed_dir, "species_niches.rds")
  microenv_path    <- file.path(processed_dir, sprintf("microenv_%s%s.rds", site_name, manifest_suffix))
  if (!file.exists(niche_cache_path) || !file.exists(microenv_path)) {
    message("No niche context for ", site_name, " -- need both ",
            niche_cache_path, " and ", microenv_path)
    return(NULL)
  }
  niche_cache <- readRDS(niche_cache_path)
  microenv    <- readRDS(microenv_path)
  heights     <- microenv_heights(microenv)

  message("Building climate cache for ", site_name, " (reads all ", length(heights), " height files once)...")
  cc <- build_clim_cache(microenv)

  height_scalars      <- height_clim_scalars(cc$clim_by_height)
  landscape_clim_vals <- height_scalars[stats::complete.cases(height_scalars), , drop = FALSE]

  niches <- load_observations()
  niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                   !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
  site_obs <- niches[niches$Area_or_Site == site_name, ]

  list(niche_cache = niche_cache, heights = heights, height_scalars = height_scalars,
       landscape_clim_vals = landscape_clim_vals, site_obs = site_obs)
}

# ── Niche suitability: per-axis scores, and the combined score before/after
# the per-site ceiling rescale ──────────────────────────────────────────────
# Pools every (species x height tier) suitability value at a site into three
# views: the three per-axis 0-100 scores (temp/relhum/swdown), the combined
# score BEFORE the per-site rescale (geometric mean of the three axis scores
# — "multiply, then normalize back to 100", prof feedback 2026-07-15), and
# the same combined score AFTER (rescaled so this species' best-scoring OWN
# observed presence height at this site reads 100 — prof feedback 2026-07-15
# follow-up; see niche_ceiling()/niche_overall_score() in get_colonization.R).
# A large gap between BEFORE's max and 100 means this site's climate profile
# never offers a height where all three axes are simultaneously ideal for
# that species, so the ceiling rescale is doing real work.
#
# extra_species: species to include even if never observed at this site
# (e.g. to preview a species_subset transplant — see run_colonization.R's
# species_file arg). Falls back to this landscape's own best-available height
# as the ceiling instead of a local observed presence — see niche_ceiling().
# context: pass a pre-built .niche_plot_context() to skip rebuilding the
# climate cache (see plot_all.R, which shares one context with
# plot_niche_profile_curves()); left NULL to build it standalone.
# See check_niche_suitability.R for the text-only per-species version.
plot_niche_suitability <- function(site_name, out_dir = OUTPUT_DIR,
                                    processed_dir = PROCESSED_DIR,
                                    height_step = 0.25,
                                    extra_species = character(0),
                                    context = NULL) {
  ctx <- context %||% .niche_plot_context(site_name, processed_dir, height_step)
  if (is.null(ctx)) {
    message("Skipping ", site_name, " niche suitability -- no context available")
    return(invisible(NULL))
  }
  niche_cache <- ctx$niche_cache; heights <- ctx$heights
  height_scalars <- ctx$height_scalars; landscape_clim_vals <- ctx$landscape_clim_vals
  site_obs <- ctx$site_obs

  clim_by_height <- lapply(seq_len(nrow(height_scalars)), function(i) as.list(height_scalars[i, ]))
  clim_by_height <- Filter(function(cl) !anyNA(unlist(cl)), clim_by_height)

  site_species <- sort(unique(c(site_obs$FinalID, extra_species)))
  site_species <- site_species[!vapply(niche_cache[site_species], is.null, logical(1))]
  if (length(site_species) == 0) {
    message("Skipping ", site_name, " niche suitability -- no cached niches for this site's species")
    return(invisible(NULL))
  }

  rows <- do.call(rbind, lapply(site_species, function(sp) {
    niche_sp <- niche_cache[[sp]]
    obs_sp   <- site_obs[site_obs$FinalID == sp & !is.na(site_obs$Height_m), ]
    obs_clim_vals <- if (nrow(obs_sp) > 0) {
      do.call(rbind, lapply(obs_sp$Height_m, function(h)
        height_scalars[which.min(abs(heights - h)), , drop = TRUE]))
    } else NULL
    ceiling <- niche_ceiling(niche_sp, obs_clim_vals, landscape_clim_vals)

    axis_scores <- sapply(clim_by_height, function(clim) niche_axis_scores(clim, niche_sp))
    before <- apply(axis_scores, 2, function(s) .geomean(s))
    after  <- pmin(100, 100 * before / ceiling)
    data.frame(species = sp, temp = axis_scores["temp", ], relhum = axis_scores["relhum", ],
               swdown = axis_scores["swdown", ], before = before, after = after)
  }))

  # Categorical slots 1/2/3 (blue/green/magenta) for the three niche axes --
  # first four slots validate all-pairs, fixed order per palette convention.
  axis_cols <- c(temp = "#2a78d6", relhum = "#008300", swdown = "#e87ba4")
  axis_df <- data.frame(
    axis  = factor(rep(c("temp", "relhum", "swdown"), each = nrow(rows)),
                   levels = c("temp", "relhum", "swdown")),
    score = c(rows$temp, rows$relhum, rows$swdown)
  )
  p_axes <- ggplot(axis_df, aes(x = score, fill = axis)) +
    geom_histogram(binwidth = 5, boundary = 0, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = axis_cols, guide = "none") +
    facet_wrap(~axis, nrow = 1) +
    labs(x = "Per-axis suitability (0-100)", y = "Height tiers x species",
         title = sprintf("Niche suitability by axis — %s", site_name)) +
    theme_minimal(base_size = 11)

  seq_blue <- "#3987e5"
  p_before <- ggplot(rows, aes(x = before)) +
    geom_histogram(binwidth = 5, boundary = 0, fill = seq_blue, color = "white", linewidth = 0.2) +
    labs(x = "Geometric mean (0-100)", y = NULL,
         title = "Combined score — before per-site ceiling rescale") +
    theme_minimal(base_size = 11)

  p_after <- ggplot(rows, aes(x = after)) +
    geom_histogram(binwidth = 5, boundary = 0, fill = seq_blue, color = "white", linewidth = 0.2) +
    labs(x = "Rescaled (0-100)", y = NULL,
         title = "Combined score — after per-site ceiling rescale") +
    theme_minimal(base_size = 11)

  p <- p_axes / (p_before + p_after)

  out_path <- file.path(out_dir, sprintf("niche_suitability_%s.png", site_name))
  ggsave(out_path, plot = p, width = 11, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Niche suitability CURVES: score vs. height, before/after the per-site
# ceiling rescale ──────────────────────────────────────────────────────────
# Companion to plot_niche_suitability() (which pools every height tier into
# histograms, losing the height axis): this draws the actual continuous
# curve of suitability against height for each species observed at a site --
# height being the natural axis for this thesis' central question (vertical
# stratification), and the "curve" framing prof feedback asked for. One
# panel BEFORE the per-site ceiling rescale (geometric mean of the three
# axis scores) and one AFTER (rescaled so the species' own best-observed
# height reads 100) -- see niche_raw_score()/niche_ceiling() in
# get_colonization.R. species_niches.rds must already exist
# (characterize_niches.R). context: see plot_niche_suitability() -- pass a
# pre-built .niche_plot_context() to skip rebuilding the climate cache.
plot_niche_profile_curves <- function(site_name, out_dir = OUTPUT_DIR,
                                       processed_dir = PROCESSED_DIR,
                                       height_step = 0.25,
                                       extra_species = character(0),
                                       context = NULL) {
  ctx <- context %||% .niche_plot_context(site_name, processed_dir, height_step)
  if (is.null(ctx)) {
    message("Skipping ", site_name, " niche profile curves -- no context available")
    return(invisible(NULL))
  }
  niche_cache <- ctx$niche_cache; heights <- ctx$heights
  height_scalars <- ctx$height_scalars; landscape_clim_vals <- ctx$landscape_clim_vals
  site_obs <- ctx$site_obs
  valid_h  <- which(stats::complete.cases(height_scalars))

  site_species <- sort(unique(c(site_obs$FinalID, extra_species)))
  site_species <- site_species[!vapply(niche_cache[site_species], is.null, logical(1))]
  if (length(site_species) == 0) {
    message("Skipping ", site_name, " niche profile curves -- no cached niches for this site's species")
    return(invisible(NULL))
  }

  rows <- do.call(rbind, lapply(site_species, function(sp) {
    niche_sp <- niche_cache[[sp]]
    obs_sp   <- site_obs[site_obs$FinalID == sp & !is.na(site_obs$Height_m), ]
    obs_clim_vals <- if (nrow(obs_sp) > 0) {
      do.call(rbind, lapply(obs_sp$Height_m, function(h)
        height_scalars[which.min(abs(heights - h)), , drop = TRUE]))
    } else NULL
    ceiling <- niche_ceiling(niche_sp, obs_clim_vals, landscape_clim_vals)

    before <- vapply(valid_h, function(i) {
      .geomean(niche_axis_scores(as.list(height_scalars[i, ]), niche_sp))
    }, numeric(1))
    after <- pmin(100, 100 * before / ceiling)

    data.frame(species = sp, height = heights[valid_h], before = before, after = after)
  }))

  long_df <- rbind(
    data.frame(species = rows$species, height = rows$height, score = rows$before,
               stage = "Before ceiling rescale"),
    data.frame(species = rows$species, height = rows$height, score = rows$after,
               stage = "After ceiling rescale")
  )
  long_df$stage <- factor(long_df$stage, levels = c("Before ceiling rescale", "After ceiling rescale"))

  # Categorical palette (fixed slot order), cycled if a site has more species
  # than slots -- most sites here have well under 8.
  cat_slots <- c("#2a78d6", "#008300", "#e87ba4", "#eda100",
                "#1baf7a", "#eb6834", "#4a3aa7", "#e34948")
  sp_cols <- setNames(cat_slots[((seq_along(site_species) - 1) %% length(cat_slots)) + 1], site_species)

  p <- ggplot(long_df, aes(x = height, y = score, colour = species)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = sp_cols, name = "Species") +
    facet_wrap(~stage, nrow = 1) +
    labs(x = "Height (m)", y = "Combined suitability (0-100)",
         title = sprintf("Niche suitability profile — %s", site_name)) +
    theme_minimal(base_size = 11)

  out_path <- file.path(out_dir, sprintf("niche_profile_curves_%s.png", site_name))
  ggsave(out_path, plot = p, width = 12, height = 5.5, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Single default colonization run (no swept parameter) ───────────────────────
plot_default_colonization_run <- function(site_name, exp_tag = "default",
                                           out_dir = OUTPUT_DIR,
                                           processed_dir = PROCESSED_DIR,
                                           animate = TRUE) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (is.data.frame(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- sweep result, use plot_colonization_experiment() instead")
    return(invisible(NULL))
  }
  # run_replicated() output: list(runs = <one runcolonization() per
  # replicate>, summary = <tidy data frame>). Falls back to treating `result`
  # itself as a single run for any older RDS saved before that change.
  runs <- if (is.list(result) && !is.null(result$runs)) result$runs else list(result)
  runs <- Filter(Negate(is.null), runs)
  if (length(runs) == 0) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no successful replicates")
    return(invisible(NULL))
  }
  message(length(runs), " replicate(s) found for ", site_name, " / ", exp_tag)

  saved <- lapply(seq_along(runs), function(i) {
    suffix <- if (length(runs) > 1) sprintf("_rep%d", i) else ""
    abundance_path <- file.path(out_dir, sprintf("abundance_%s_%s%s.png", site_name, exp_tag, suffix))
    png(abundance_path, width = 1800, height = 900, res = 150, type = "cairo")
    plot_abundance(runs[[i]])
    dev.off()
    message("Saved: ", abundance_path)
    abundance_path
  })

  # 3D: static snapshot (final year) + animated (all years) for the first
  # replicate only, to avoid generating one large HTML per replicate by default.
  fig_3d <- plot_3d_abundance(runs[[1]])
  volume_path <- file.path(out_dir, sprintf("abundance_3d_%s_%s.html", site_name, exp_tag))
  if (!is.null(fig_3d)) {
    # selfcontained=TRUE needs pandoc (not installed on the cluster); FALSE
    # writes a small "<name>_files/" dependency folder alongside the HTML
    # instead -- keep the two together when copying/viewing elsewhere.
    htmlwidgets::saveWidget(fig_3d, volume_path, selfcontained = FALSE)
    message("Saved: ", volume_path)
  }

  anim_path <- NULL
  if (animate) {
    anim_path <- file.path(out_dir, sprintf("abundance_3d_animated_%s_%s.html", site_name, exp_tag))
    plot_3d_abundance_animated(runs[[1]], out_path = anim_path)
  }

  invisible(list(abundance_png = saved, abundance_3d = fig_3d, abundance_3d_animated = anim_path))
}
