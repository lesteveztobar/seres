# plot_functions.R
# Reusable plotting functions for the whole seres pipeline. Each
# function loads its own inputs from data/processed (or geojson_to_csv) and
# saves its output(s) into OUTPUT_DIR — nothing here plots-and-forgets.
# If an input file doesn't exist yet, the function skips with a message
# instead of erroring, so plot_all.R can be re-run at any point in the
# pipeline and it just draws whatever is available so far.
#
# Refactored from the standalone plotmap.R / plot_temperatures.R /
# plot_bestfit_3d.R, plus .plot_live()/plot_abundance()/.canopy_context_
# trace()/plot_3d_abundance()/plot_3d_abundance_animated() (moved here from
# get_colonization.R -- pure visualization, no climate/niche logic) and
# export wrappers around them, plus plot_experiment() which is still defined
# in get_colonization.R.
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2)
library(patchwork)
library(plotly)
library(htmlwidgets)
library(abind)
library(scatterplot3d)
library(scico)

source("scripts/02_model/config/paths.R")
source("scripts/02_model/config/shared_helpers.R")
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
  # into their own GeoJSON exports (see scripts/data_prep/rebuild_combined_csv.py).
  # Left on disk as an audit trail but no longer read here.
  # Saloya excluded from this map (2026-08-29, matches batch_exp.sh's
  # 2026-08-28 exclusion from the colonization sensitivity experiments) --
  # this figure is meant to show the sites the model was actually run at.
  geojson_files <- c(
    Maquipucuna   = "Maquipucuna.geojson",
    Mashpi        = "Mashpi.geojson",
    MindoTarabita = "MindoTarabita.geojson",   # merged with TarabitaMindo below
    TarabitaMindo = "TarabitaMindo.geojson",
    LaElenita     = "LaElenita.geojson",
    MindoMirador  = "MindoMirador.geojson",
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

  # 2026-08-29: this figure now shows observation points only -- no transect
  # lines (removed per editorial review: too small/short at this figure's
  # spatial scope to be legible, and only added visual clutter). Every
  # site's label is therefore positioned from its points' own centroid;
  # the transect-vs-point label-source split this used to need is gone.
  obs_points <- all_features %>% filter(sf::st_geometry_type(geometry) == "POINT")

  site_labels <- obs_points %>%
    group_by(site) %>%
    summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    mutate(geometry = sf::st_centroid(geometry))

  # Per-site label offsets (degrees) — keeps long names off the dots.
  # LaElenita/MindoMirador added 2026-07-24 with nudge (0,0) placeholders
  # (untuned) -- adjust once you've seen how their labels sit on the actual
  # map.
  label_nudges <- data.frame(
    site    = c("Mashpi", "Maquipucuna", "MindoTarabita", "Yanayacu",
                "LaElenita", "MindoMirador"),
    nudge_x = c( 0.00,     0.08,          0.15,            0.10,
                 0.00,      0.00),
    nudge_y = c(-0.04,     0.07,         -0.05,            0.06,
                 0.00,      0.00)
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
    scico::scico(6, palette = "lipari", begin = 0.10, end = 0.88),
    c("Maquipucuna", "Mashpi", "MindoTarabita", "LaElenita", "MindoMirador", "Yanayacu")
  )

  # ── Mindo cluster inset ───────────────────────────────────────────────────
  # MindoTarabita / LaElenita / MindoMirador fall within ~7 km of each other
  # and collapse to a single smudge at the regional scale. Draw a zoomed panel
  # over the empty SW corner of the main map (no sites there) so the three
  # read separately, and outline the same extent on the main map so the reader
  # can locate it. Inset labels are placed from each site's own point centroid
  # (the regional `label_nudges` are far too large for this ~3 km window).
  mindo_sites  <- c("MindoTarabita", "LaElenita", "MindoMirador")
  mindo_points <- obs_points %>% filter(site %in% mindo_sites)
  mindo_bbox   <- sf::st_bbox(mindo_points)
  # Zoom extent: pad in y, then set x from the box centre to whatever half-width
  # gives the inset panel a roughly square aspect (the raw cluster is a tall
  # sliver -- ~2 km E-W vs ~7 km N-S -- and would render as an unreadable
  # ribbon at equal scale).
  pad_y      <- max(as.numeric(mindo_bbox["ymax"] - mindo_bbox["ymin"]) * 0.28, 0.006)
  mindo_cx   <- mean(c(mindo_bbox["xmin"], mindo_bbox["xmax"]))
  half_y     <- as.numeric(mindo_bbox["ymax"] - mindo_bbox["ymin"]) / 2 + pad_y
  half_x     <- half_y * 1.4
  mindo_xlim <- c(mindo_cx - half_x, mindo_cx + half_x)
  mindo_ylim <- c(mindo_bbox["ymin"] - pad_y, mindo_bbox["ymax"] + pad_y)

  mindo_labels <- site_labels %>%
    dplyr::filter(site %in% mindo_sites) %>%
    dplyr::mutate(X = sf::st_coordinates(geometry)[, 1],
                  Y = sf::st_coordinates(geometry)[, 2]) %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(
      lab_hjust = dplyr::if_else(site == "LaElenita", 1, 0),
      X_lab     = X + dplyr::if_else(site == "LaElenita", -half_x * 0.05, half_x * 0.05)
    )

  p_inset <- ggplot() +
    geom_sf(data = map_area, fill = "#f5f0e8", colour = "grey55", linewidth = 0.3) +
    geom_sf(data = mindo_points, aes(colour = site), size = 2.6, alpha = 0.9,
            shape = 16, show.legend = FALSE) +
    geom_text(data = mindo_labels,
              aes(x = X_lab, y = Y, label = site, colour = site, hjust = lab_hjust),
              size = 2.3, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = site_colours) +
    coord_sf(xlim = mindo_xlim, ylim = mindo_ylim, expand = FALSE) +
    annotation_scale(location = "br", width_hint = 0.3,
                     height = unit(0.1, "cm"), text_cex = 0.55,
                     pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm")) +
    labs(title = "Mindo sites (zoom)") +
    theme_bw(base_size = 8) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      plot.title       = element_text(size = 7.5, face = "bold"),
      plot.background  = element_rect(fill = "white", colour = "grey30", linewidth = 0.6),
      plot.margin      = margin(3, 3, 3, 3),
      panel.grid       = element_blank()
    )

  # Main-map labels: LaElenita / MindoMirador are unreadably stacked on
  # MindoTarabita at this scale -- the inset carries them, so the main map
  # just points to it.
  main_labels <- label_pos %>%
    dplyr::filter(!site %in% c("LaElenita", "MindoMirador")) %>%
    dplyr::mutate(lab = dplyr::if_else(site == "MindoTarabita",
                                       "Mindo sites (see inset)", site))

  p <- ggplot() +
    geom_sf(data = map_area, fill = "#f5f0e8", colour = "grey55", linewidth = 0.35) +
    geom_sf(data = obs_points, aes(colour = site), size = 1.8, alpha = 0.55, shape = 16) +
    annotate("rect", xmin = mindo_xlim[1], xmax = mindo_xlim[2],
             ymin = mindo_ylim[1], ymax = mindo_ylim[2],
             fill = NA, colour = "grey30", linewidth = 0.5) +
    geom_text(data = main_labels, aes(x = X, y = Y, label = lab, colour = site),
              size = 3, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = site_colours, name = "Site") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "br", width_hint = 0.25) +
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

  p <- p + patchwork::inset_element(
    p_inset, left = 0.015, bottom = 0.02, right = 0.42, top = 0.44,
    align_to = "panel"
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
  tree_height <- forestparams$mean_hgt
  crown_r <- forestparams$mean_crown_r

  heights      <- seq(0.5, tree_height, by = res)
  half_extent  <- ceiling(crown_r / res) * res
  xy           <- seq(-half_extent, half_extent, by = res)

  rows <- list()
  for (h in heights) {
    rel_height <- h / tree_height
    # Same geometry as build_forest() (get_colonization.R) -- calls its
    # .classify_tree_zone()/.crown_fraction()/.effective_crown_radius()
    # helpers directly rather than re-deriving the thresholds/formula here,
    # so this stays in sync automatically.
    zone_id <- .classify_tree_zone(rel_height)
    if (zone_id <= 2) {
      rows[[length(rows) + 1]] <- data.frame(x = 0, y = 0, z = h, zone = zone_id)
    } else {
      effective_r <- .effective_crown_radius(crown_r, rel_height)
      if (effective_r <= 0) next
      grid      <- expand.grid(x = xy, y = xy)
      grid$dist <- sqrt(grid$x^2 + grid$y^2)
      grid      <- grid[.voxel_in_tree(zone_id, grid$dist, effective_r), ]
      if (nrow(grid) == 0) next
      rows[[length(rows) + 1]] <- data.frame(x = grid$x, y = grid$y, z = h, zone = zone_id)
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
    main = sprintf("Modelled tree geometry\n(height = %.1f m, crown radius = %.1f m)", tree_height, crown_r),
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
                                      height_step = 0.4) {
  # height_step defaults to 0.4 (the production resolution as of
  # 2026-08-17, was 0.25 -- used everywhere else in this pipeline --
  # run_colonization.R, characterize_niches.R -- via the same
  # manifest_suffix convention), NOT the unsuffixed 0.1m file.
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
  # per-timestep spatial array into memory at once -- raw hourly per-pixel
  # arrays no longer exist on disk at all as of the 2026-08-05 write-time
  # quantile-reduction move (get_colonization.R's .compute_voxel_quantiles()/
  # get_clim_voxel()); each height file now carries per-pixel quantiles
  # instead. The per-pixel annual median temperature (the "annual"/"both"
  # daypart quantile's 50th percentile, one value per raster pixel) is this
  # function's closest available equivalent to the old per-pixel annual
  # MEAN -- both are a single robust central-tendency value per pixel, and
  # this still needs (and gets) full per-pixel spatial resolution, unlike
  # the separately-flattened lookup_climate_by_height() pathway.
  tz_layers <- lapply(heights_m, function(hgt) {
    h <- load_height(env, hgt)
    vq <- h$voxel_quantiles
    matrix(vq$quantiles[["annual_both_temp"]][, 3], nrow = vq$nr, ncol = vq$nc)
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
    sig        = .sig_stars(dunn_res$P.adjusted)  # shared_helpers.R
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
                                      out_dir = OUTPUT_DIR, processed_dir = PROCESSED_DIR,
                                      metric = c("extinction", "survival")) {
  metric <- match.arg(metric)
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
  shape <- .classify_result_shape(result)  # shared_helpers.R
  swept <- shape$swept
  if ("n_founders" %in% swept) swept <- c("n_founders", setdiff(swept, "n_founders"))
  if (length(swept) < 2) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- fewer than 2 swept columns found (", paste(swept, collapse = ", "), ")")
    return(invisible(NULL))
  }

  t_max <- shape$t_max
  combo_summary <- shape$combo
  # combo$extinct is already the per-combo MEAN extinction rate (extinct is a
  # 0/1 flag per replicate, aggregated with FUN=mean in .classify_result_
  # shape()) -- so 1-extinct is directly the survival rate, not just a
  # majority-vote boolean like combo$persisted.
  combo_summary$survival <- 1 - combo_summary$extinct

  x_var <- swept[1]; y_var <- swept[2]
  facet_vars <- swept[-(1:2)]

  fill_var   <- if (metric == "survival") "survival" else "extinct"
  legend_lab <- if (metric == "survival") "Survival\nrate" else "Extinction\nrate"
  subtitle_lab <- if (metric == "survival") "Survival rate" else "Extinction rate"
  # Reversed ramp for survival so high (good) still reads as the "safe" blue
  # end and low (bad) as red, matching the extinction ramp's color sense
  # rather than just flipping the number and keeping red-for-high.
  ramp <- if (metric == "survival")
    c("#a50026", "#f1a340", "#08519c") else c("#08519c", "#f1a340", "#a50026")

  p <- ggplot(combo_summary,
             aes(x = factor(.data[[x_var]]), y = factor(.data[[y_var]]), fill = .data[[fill_var]])) +
    geom_tile() +
    scale_fill_gradientn(colours = ramp, limits = c(0, 1), name = legend_lab) +
    labs(x = x_var, y = y_var,
         title = sprintf("Factorial — %s", site_name),
         subtitle = if (length(facet_vars) > 0)
           sprintf("%s at year %d, faceted by %s",
                   subtitle_lab, t_max, paste(facet_vars, collapse = " x "))
         else
           sprintf("%s at year %d", subtitle_lab, t_max)) +
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

  suffix <- if (metric == "survival") "_survival" else ""
  out_path <- file.path(out_dir, sprintf("factorial_%s_%s%s.png", site_name, exp_tag, suffix))
  ggsave(out_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# Every pairwise combination of the 4 reproduction_factorial_v3 swept
# params, each as its own small heatmap, filled by % change in final
# abundance relative to the literature-realistic baseline -- a sibling to
# plot_factorial_experiment() above (which stays as-is: 2 params as tile
# axes + the rest faceted, filled by extinction/survival), not a
# replacement. Answers "which parameter(s) actually move abundance",
# distinct from plot_factorial_experiment()'s "where does this specific
# 2-axis slice go extinct". 2026-08-25.
plot_factorial_pairwise_heatmap <- function(site_name, exp_tag = "reproduction_factorial_v3",
                                            out_dir = OUTPUT_DIR, processed_dir = PROCESSED_DIR) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (!is.data.frame(result) || "param_value" %in% names(result)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- not a factorial result")
    return(invisible(NULL))
  }
  shape <- .classify_result_shape(result)  # shared_helpers.R
  swept <- shape$swept
  if (length(swept) < 2) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- fewer than 2 swept columns found (", paste(swept, collapse = ", "), ")")
    return(invisible(NULL))
  }
  # Subtitle below is hardcoded to reproduction_factorial_v3's own 4 params
  # (this function's documented target) -- guard against a differently-
  # shaped factorial (e.g. only 2-3 swept cols) hitting an undefined field.
  if (!all(c("p_poll", "p_germ", "p_s1", "n_founders") %in% swept)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- plot_factorial_pairwise_heatmap() expects exactly p_poll/p_germ/p_s1/n_founders, found: ",
            paste(swept, collapse = ", "))
    return(invisible(NULL))
  }
  t_max <- shape$t_max
  final <- result[result$t == t_max, ]

  # Baseline = the literature-realistic combo -- confirmed (make_params.R,
  # params_reprofactorial_v3) to be exactly the FIRST level of every swept
  # param (p_poll=0.30, p_germ=0.00100, p_s1=0.450, n_founders=30), i.e. a
  # literal row already in this data -- no separate "realistic" run needed.
  baseline_vals <- setNames(lapply(swept, function(v) sort(unique(final[[v]]))[1]), swept)
  baseline_mask <- Reduce(`&`, lapply(swept, function(v) final[[v]] == baseline_vals[[v]]))
  if (!any(baseline_mask)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- literature-realistic baseline combo not found in this data")
    return(invisible(NULL))
  }
  baseline_abundance <- mean(final$total[baseline_mask])

  # Best combo -- same aggregate()+which.max logic as pick_best_combo.R,
  # reused rather than reimplemented.
  combo_form <- as.formula(paste("total ~", paste(swept, collapse = " + ")))
  combo_agg  <- aggregate(combo_form, data = final, FUN = mean)
  best_combo <- combo_agg[which.max(combo_agg$total), ]

  # effect_size (2026-08-28, replaces raw pct_change): 100 * mean_total /
  # best_combo$total. Naturally bounded to [0, 100] with no clamping needed
  # -- abundance is never negative, and best_combo$total is by construction
  # the maximum mean-total among all 625 tested combos for this site, so
  # every other combo's total is <= it. This fixes two things at once: (1)
  # pct_change was unbounded and each of the 6 panels trained its own
  # independent color scale (plot_layout(guides="collect") only visually
  # merges legends that already match, so it was silently misrepresenting
  # 6 different scales as one) -- fixed limits=c(0,100) make the scale
  # trivially, truthfully shared across every panel; (2) reusing the same
  # scico "lipari" 0-100 ramp as the niche-suitability figures unifies the
  # report's visual language for "0-100 suitability-style" scores.
  pairs <- combn(swept, 2, simplify = FALSE)
  panels <- lapply(pairs, function(pr) {
    A <- pr[1]; B <- pr[2]
    form <- as.formula(paste("total ~", A, "+", B))
    agg  <- aggregate(form, data = final, FUN = mean)
    agg$effect_size <- 100 * agg$total / best_combo$total
    agg$is_best <- agg[[A]] == best_combo[[A]] & agg[[B]] == best_combo[[B]]
    agg$is_baseline <- agg[[A]] == baseline_vals[[A]] & agg[[B]] == baseline_vals[[B]]

    ggplot(agg, aes(x = factor(.data[[A]]), y = factor(.data[[B]]), fill = effect_size)) +
      geom_tile() +
      geom_tile(data = agg[agg$is_baseline, ], fill = NA, colour = "white",
               linewidth = 1, linetype = "dashed") +
      geom_tile(data = agg[agg$is_best, ], fill = NA, colour = "black", linewidth = 1) +
      scale_fill_gradientn(colours = scico::scico(100, palette = "lipari"),
                           limits = c(0, 100), name = "Effect size\n(0-100)") +
      labs(x = A, y = B) +
      theme_minimal(base_size = 10)
  })

  p <- patchwork::wrap_plots(panels, nrow = 2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Factorial parameter effects — %s", site_name),
      subtitle = sprintf(
        "Effect size: 0 = extinct, 100 = this site's best-tested combo (black outline) | Baseline (literature-realistic, dashed white outline): p_poll=%.2f, p_germ=%.4f, p_s1=%.2f, n_founders=%d -> N=%.1f (effect size %.0f) | Best: p_poll=%.2f, p_germ=%.4f, p_s1=%.2f, n_founders=%d -> N=%.1f",
        baseline_vals$p_poll, baseline_vals$p_germ, baseline_vals$p_s1, baseline_vals$n_founders,
        baseline_abundance, 100 * baseline_abundance / best_combo$total,
        best_combo$p_poll, best_combo$p_germ, best_combo$p_s1, best_combo$n_founders, best_combo$total))

  out_path <- file.path(out_dir, sprintf("factorial_pairwise_%s_%s.png", site_name, exp_tag))
  ggsave(out_path, plot = p, width = 16, height = 9, dpi = 300, bg = "white")
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
.niche_plot_context <- function(site_name, processed_dir = PROCESSED_DIR, height_step = 0.4,
                                niche_cache_path = NICHE_CACHE_PATH) {
  manifest_suffix  <- if (height_step != 0.1) sprintf("_h%.2f", height_step) else ""
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
  cc <- build_clim_cache_voxel(microenv)

  # No landscape/footprint defined for this diagnostic (no simulated grid
  # the way init_colonization() has) -- footprint = NULL pools per-height
  # quantiles across the whole raster, same one-row-per-height shape
  # height_clim_scalars() used to produce (see get_colonization.R).
  height_scalars      <- voxel_background_table(cc, microenv, footprint = NULL, months = "annual")
  landscape_clim_vals <- height_scalars[stats::complete.cases(height_scalars), , drop = FALSE]

  niches <- load_observations()
  niches <- niches[!is.na(niches$lat) & !is.na(niches$lon) &
                   !is.na(niches$Height_m) & !is.na(niches$FinalID), ]
  site_obs <- niches[niches$Area_or_Site == site_name, ]

  # `niches` (the full, un-site-filtered table) is carried through so callers
  # can check genuine identification via .confirmed_species_sites()
  # (shared_helpers.R), which needs the raw Identification column across
  # every site, not just this one -- see .niche_score_rows()'s species
  # filter fix, 2026-08-28.
  list(niche_cache = niche_cache, heights = heights, height_scalars = height_scalars,
       landscape_clim_vals = landscape_clim_vals, site_obs = site_obs, niches = niches)
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
                                    height_step = 0.4,
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

  # See .niche_score_rows() for why this filters on genuine identification
  # (.confirmed_species_sites(), shared_helpers.R) rather than the
  # post-default FinalID.
  candidate_species <- sort(unique(c(site_obs$FinalID, extra_species)))
  site_species <- candidate_species[
    candidate_species %in% extra_species |
    vapply(candidate_species, function(sp) site_name %in% .confirmed_species_sites(sp, ctx$niches), logical(1))
  ]
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
# ── Shared per-species x per-height suitability score matrix ───────────────
# Factored out of plot_niche_profile_curves() (2026-08-24) so the new
# heatmap version (plot_niche_suitability_heatmap()) computes the exact same
# numbers instead of duplicating the niche-scoring math -- both just re-
# encode this same (species, height, before, after) table differently
# (geom_line vs geom_tile). Returns NULL (with a message()) for the same
# "no context" / "no cached niches" cases the line-plot already handles, so
# both callers can share one skip-and-continue check.
.niche_score_rows <- function(ctx, site_name, extra_species = character(0)) {
  niche_cache <- ctx$niche_cache; heights <- ctx$heights
  height_scalars <- ctx$height_scalars; landscape_clim_vals <- ctx$landscape_clim_vals
  site_obs <- ctx$site_obs
  valid_h  <- which(stats::complete.cases(height_scalars))

  # 2026-08-28: filter site_obs$FinalID down to species this site actually,
  # genuinely had identified -- not the post-default FinalID, which
  # silently relabels every blank/unverified-photo row as "Maxillaria
  # acutifolia" (load_observations(), paths.R). Without this, a site with
  # mostly-unidentified observations (e.g. LaElenita/MindoMirador/Saloya)
  # borrows a well-supported niche model built from OTHER sites' genuine
  # acutifolia sightings, making its suitability figure look convincing for
  # reasons that have nothing to do with that site's own plants. Reuses
  # .confirmed_species_sites() (shared_helpers.R) rather than duplicating
  # its raw-Identification-column logic. extra_species (explicitly
  # requested by the caller) is exempt from this filter.
  candidate_species <- sort(unique(c(site_obs$FinalID, extra_species)))
  site_species <- candidate_species[
    candidate_species %in% extra_species |
    vapply(candidate_species, function(sp) site_name %in% .confirmed_species_sites(sp, ctx$niches), logical(1))
  ]
  site_species <- site_species[!vapply(niche_cache[site_species], is.null, logical(1))]
  if (length(site_species) == 0) {
    message("Skipping ", site_name, " -- no cached niches for this site's species")
    return(NULL)
  }
  if (length(valid_h) == 0) {
    # 2026-08-25: confirmed pre-existing, NOT specific to this site or
    # introduced by this refactor -- voxel_background_table(footprint=NULL)
    # (get_colonization.R) synthesizes NA lon/lat "observations" intending
    # to pool across every pixel per height, but .voxel_point_quantile()'s
    # pixel-mode branch (the mode every real production microenv uses) has
    # no such pooling path for NA coordinates -- .lonlat_to_pixel(NA, NA)
    # never matches any pixel, so height_scalars comes back entirely NA for
    # every site tested (Maquipucuna, MindoMirador). Without this check,
    # `data.frame(species=sp, height=heights[valid_h], before=before,
    # after=after)` recycles the length-1 `species=sp` against length-0
    # height/before/after and errors ("differing number of rows: 1, 0").
    message("Skipping ", site_name,
      " -- no height tier has complete landscape climate data (see voxel_background_table() NULL-footprint note)")
    return(NULL)
  }

  do.call(rbind, lapply(site_species, function(sp) {
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
}

plot_niche_profile_curves <- function(site_name, out_dir = OUTPUT_DIR,
                                       processed_dir = PROCESSED_DIR,
                                       height_step = 0.4,
                                       extra_species = character(0),
                                       confirmed_only = TRUE,
                                       context = NULL) {
  ctx <- context %||% .niche_plot_context(site_name, processed_dir, height_step)
  if (is.null(ctx)) {
    message("Skipping ", site_name, " niche profile curves -- no context available")
    return(invisible(NULL))
  }
  rows <- .niche_score_rows(ctx, site_name, extra_species)
  if (is.null(rows)) {
    message("Skipping ", site_name, " niche profile curves -- no cached niches for this site's species")
    return(invisible(NULL))
  }
  if (confirmed_only) {
    rows <- rows[.niche_is_confirmed(rows$species) | rows$species %in% extra_species, , drop = FALSE]
  }
  if (nrow(rows) == 0) {
    message("Skipping ", site_name, " niche profile curves -- no confirmed species")
    return(invisible(NULL))
  }

  # One small panel per species (niche-height order), two lines: combined
  # suitability along canopy height before vs. after the per-species ceiling
  # rescale. Faceting by species instead of overplotting every species in
  # two shared panels -- the old form was an unreadable hairball once the
  # after-rescale curves saturate near 100. Heights binned for display, same
  # as the all-sites niche figures.
  sp_order <- .niche_species_order(rows)
  long_df <- .niche_long_df(rows, c("before", "after"))
  long_df <- .niche_bin_heights(long_df, "score", .NICHE_DISP_BIN)
  long_df$species <- factor(.niche_sp_abbr(long_df$species), levels = .niche_sp_abbr(sp_order))
  long_df$stage <- factor(long_df$stage, levels = c("Before ceiling rescale", "After ceiling rescale"))

  stage_cols <- setNames(scico::scico(2, palette = "lipari", begin = 0.30, end = 0.72),
                         c("Before ceiling rescale", "After ceiling rescale"))

  n_sp <- nlevels(long_df$species)
  ncol_f <- min(4, n_sp)

  p <- ggplot(long_df, aes(x = height, y = score, colour = stage)) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = stage_cols, name = NULL) +
    facet_wrap(~species, ncol = ncol_f) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = "Height above ground (m)", y = "Combined suitability (0–100)",
         title = sprintf("Niche suitability profile — %s", site_name)) +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "italic"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")

  w <- max(7, 1.4 + ncol_f * 2.3)
  h <- max(4, 1.6 + ceiling(n_sp / ncol_f) * 2.1)
  out_path <- file.path(out_dir, sprintf("niche_profile_curves_%s.png", site_name))
  ggsave(out_path, plot = p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  message("Saved: ", out_path)
  invisible(p)
}

# Species (y) x height (x) suitability heatmap -- the "hotbox" version of
# plot_niche_profile_curves() above: same (species, height, before/after)
# score table via .niche_score_rows(), re-encoded as geom_tile instead of
# geom_line so every species x height cell's suitability is readable at a
# glance instead of needing to trace overlapping lines. No continuous
# suitability fill convention existed anywhere in this codebase before this
# (confirmed 2026-08-24) -- scico "lipari" chosen to match the discrete
# species/stage palette already used everywhere else (.species_colors()/
# .stage_colors(), shared_helpers.R), rather than introducing an unrelated
# palette family just for this one plot.
plot_niche_suitability_heatmap <- function(site_name, out_dir = OUTPUT_DIR,
                                           processed_dir = PROCESSED_DIR,
                                           height_step = 0.4,
                                           extra_species = character(0),
                                           context = NULL) {
  ctx <- context %||% .niche_plot_context(site_name, processed_dir, height_step)
  if (is.null(ctx)) {
    message("Skipping ", site_name, " niche suitability heatmap -- no context available")
    return(invisible(NULL))
  }
  rows <- .niche_score_rows(ctx, site_name, extra_species)
  if (is.null(rows)) {
    message("Skipping ", site_name, " niche suitability heatmap -- no cached niches for this site's species")
    return(invisible(NULL))
  }

  long_df <- rbind(
    data.frame(species = rows$species, height = rows$height, score = rows$before,
               stage = "Before ceiling rescale"),
    data.frame(species = rows$species, height = rows$height, score = rows$after,
               stage = "After ceiling rescale")
  )
  long_df$stage <- factor(long_df$stage, levels = c("Before ceiling rescale", "After ceiling rescale"))

  p <- ggplot(long_df, aes(x = height, y = species, fill = score)) +
    geom_tile() +
    scale_fill_gradientn(colours = scico::scico(100, palette = "lipari"),
                         limits = c(0, 100), name = "Suitability\n(0-100)") +
    facet_wrap(~stage, nrow = 1) +
    labs(x = "Height (m)", y = "Species",
         title = sprintf("Niche suitability hotbox — %s", site_name)) +
    theme_minimal(base_size = 11)

  out_path <- file.path(out_dir, sprintf("niche_suitability_heatmap_%s.png", site_name))
  ggsave(out_path, plot = p, width = 12, height = 5.5, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── Combined, all-sites niche-suitability data prep ─────────────────────────
# Shared by every figure in the "Combined all-sites niche-suitability
# figures" block below (plot_niche_suitability_by_site(), plot_niche_
# forest(), plot_niche_multisite_species()) -- all need the exact same
# (site, species, height, before/after) long table, just re-encoded with
# different aesthetics. Callers normally reach it through
# .niche_score_rows_all_sites_cached(). Builds each site's context
# fresh (.niche_plot_context()) and calls the now-fixed .niche_score_rows()
# (genuine-identification filter), tags rows with `site`, and drops/reports
# any site left with zero confirmed species after that filter -- same
# skip-and-message convention as everywhere else in this file.
.niche_score_rows_all_sites <- function(sites = NULL, processed_dir = PROCESSED_DIR,
                                         height_step = 0.4) {
  if (is.null(sites)) {
    files <- list.files(processed_dir, pattern = sprintf("^microenv_.*_h%.2f\\.rds$", height_step))
    sites <- sub(sprintf("^microenv_(.*)_h%.2f\\.rds$", height_step), "\\1", files)
  }
  excluded <- character(0)
  rows_by_site <- lapply(sites, function(site_name) {
    ctx <- .niche_plot_context(site_name, processed_dir, height_step)
    if (is.null(ctx)) {
      excluded[[length(excluded) + 1]] <<- site_name
      return(NULL)
    }
    rows <- .niche_score_rows(ctx, site_name)
    if (is.null(rows)) {
      excluded[[length(excluded) + 1]] <<- site_name
      return(NULL)
    }
    cbind(site = site_name, rows)
  })
  all_rows <- do.call(rbind, Filter(Negate(is.null), rows_by_site))
  if (length(excluded) > 0) {
    message("Excluded (no confirmed species): ", paste(excluded, collapse = ", "))
  }
  list(rows = all_rows, excluded = excluded)
}

# Cached wrapper: .niche_score_rows_all_sites() is a full per-site
# climate-cache pass (~30 min for 7 sites) and every figure in the family
# below needs the exact same table. Cache it to data/processed/ keyed by the
# active observations CSV + niche cache (so the default and v6 datasets get
# separate caches) and invalidate whenever any input file (CSV, niche cache,
# background, or a microenv_*_h<step>.rds) is newer than the cache. Callers
# that already hold a `built` list pass it straight through.
.niche_score_rows_all_sites_cached <- function(processed_dir = PROCESSED_DIR,
                                               height_step = 0.4) {
  microenv <- list.files(processed_dir,
    pattern = sprintf("^microenv_.*_h%.2f\\.rds$", height_step), full.names = TRUE)
  inputs <- c(OBSERVATIONS_CSV, NICHE_CACHE_PATH, NICHE_BACKGROUND_PATH, microenv)
  inputs <- inputs[file.exists(inputs)]
  tag <- paste0(tools::file_path_sans_ext(basename(OBSERVATIONS_CSV)), "__",
                tools::file_path_sans_ext(basename(NICHE_CACHE_PATH)))
  cache <- file.path(processed_dir,
    sprintf("niche_score_rows__%s__h%.2f.rds", tag, height_step))
  if (file.exists(cache) && length(inputs) > 0 &&
      file.mtime(cache) >= max(file.mtime(inputs))) {
    message("Reusing cached niche score table: ", cache)
    return(readRDS(cache))
  }
  built <- .niche_score_rows_all_sites(processed_dir = processed_dir, height_step = height_step)
  saveRDS(built, cache)
  message("Wrote niche score table cache: ", cache)
  built
}

# ── Combined all-sites niche-suitability figures ──────────────────────────────
# 2026-08-29 rewrite. The previous single figure (site x height per species,
# facet_grid(species ~ stage)) was ~35 species rows x 12 columns and mostly
# empty -- each species occurs at only 1-2 sites -- so it rendered ~40 in
# tall and unreadable. Replaced with a small family of focused figures that
# share one builder (.niche_score_rows_all_sites) and helper set:
#
#   plot_niche_suitability_by_site()  Fig A -- per-site small multiples, the
#                                     quantitative comparison figure.
#   plot_niche_forest()               Fig B -- per-site bars on a fixed
#                                     global species axis (gaps = absence).
#   plot_niche_multisite_species()    Supp -- only species confirmed at 2+
#                                     sites, horizontal bars faceted per species.
#   plot_niche_suitability_by_site(stages = c("before","after"))
#                                     Supp -- the ceiling-rescale effect.
#   plot_niche_suitability_by_site(confirmed_only = FALSE)
#                                     Supp -- unidentified morphospecies.

# A species is "confirmed" when its Identification carries a real epithet --
# anything of the form "<Genus> sp."/"sp"/"sp. A"/"sp1"/"sp2" is an
# unidentified morphospecies. Kept deliberately simple (token test on the
# epithet) rather than reusing .confirmed_species_sites(), which answers a
# different question (was THIS row genuinely identified, vs. a blank-ID
# default) -- here we're binning a species name, not a site membership.
.niche_is_confirmed <- function(species) {
  epithet <- sub("^\\S+\\s*", "", trimws(species))
  nzchar(epithet) & !grepl("^sp(\\.|[0-9]|\\s|$)", epithet, ignore.case = TRUE)
}

# "Maxillaria acutifolia" -> "M. acutifolia" (axis labels; drawn italic).
.niche_sp_abbr <- function(x) sub("^([A-Z])[a-z]+\\s+", "\\1. ", x)

# Long (site, species, height, score, stage) table from a .niche_score_rows_
# all_sites() result, for whichever ceiling-rescale stage(s) are asked for.
.niche_long_df <- function(rows, stages = c("before", "after")) {
  parts <- list()
  if ("before" %in% stages)
    parts$b <- data.frame(site = rows$site, species = rows$species, height = rows$height,
                          score = rows$before, stage = "Before ceiling rescale")
  if ("after" %in% stages)
    parts$a <- data.frame(site = rows$site, species = rows$species, height = rows$height,
                          score = rows$after, stage = "After ceiling rescale")
  df <- do.call(rbind, parts)
  df$stage <- factor(df$stage, levels = c("Before ceiling rescale", "After ceiling rescale"))
  df
}

# Order species by where their niche actually sits in the canopy: the
# suitability-weighted mean height (after-rescale, pooled across sites), so
# every figure lists species low-canopy -> high-canopy and the height
# gradient reads as a diagonal. Species with no positive score anywhere fall
# back to their median modelled height.
.niche_species_order <- function(rows) {
  agg <- tapply(seq_len(nrow(rows)), rows$species, function(ix) {
    w <- pmax(rows$after[ix], 0); h <- rows$height[ix]
    if (sum(w) > 0) sum(w * h) / sum(w) else stats::median(h)
  })
  names(sort(unlist(agg)))
}

.NICHE_FILL <- function()
  scale_fill_gradientn(colours = scico::scico(100, palette = "lipari"),
                       limits = c(0, 100), name = "Climate-niche suitability (0–100)")

# Collapse the per-tier suitability to coarser height bins (mean) for
# display. The raw 0.4 m tiers are spiky enough that the tiled figures read
# as noise; ~1.5 m bins keep the vertical structure without the strobing.
# `by` names the non-height grouping columns present in `df`.
.niche_bin_heights <- function(df, value = "score", bin = 1.5,
                               by = c("site", "species", "stage")) {
  by <- intersect(by, names(df))
  df$height <- (floor(df$height / bin) + 0.5) * bin
  form <- stats::as.formula(sprintf("%s ~ %s", value, paste(c(by, "height"), collapse = " + ")))
  stats::aggregate(form, data = df, FUN = mean)
}
.NICHE_DISP_BIN <- 1.5

# ── Fig A: per-site suitability-by-height small multiples ─────────────────────
# One panel per site; x = species (ordered by niche height, italic), y =
# height above ground, fill = suitability. Defaults to confirmed species and
# the after-ceiling-rescale stage only (the headline number). Pass
# stages = c("before","after") for the rescale-effect supplement, or
# confirmed_only = FALSE for the morphospecies supplement.
plot_niche_suitability_by_site <- function(sites = NULL, out_dir = OUTPUT_DIR,
                                           processed_dir = PROCESSED_DIR,
                                           height_step = 0.4, built = NULL,
                                           confirmed_only = TRUE,
                                           stages = "after", file_tag = NULL) {
  built <- built %||% .niche_score_rows_all_sites(sites, processed_dir, height_step)
  rows  <- built$rows
  if (is.null(rows) || nrow(rows) == 0) {
    message("Skipping niche suitability by site -- no sites had confirmed species")
    return(invisible(NULL))
  }
  is_conf <- .niche_is_confirmed(rows$species)
  rows <- rows[if (confirmed_only) is_conf else !is_conf, , drop = FALSE]
  if (nrow(rows) == 0) {
    message("Skipping niche suitability by site -- no ",
            if (confirmed_only) "confirmed" else "morphospecies", " rows")
    return(invisible(NULL))
  }

  sp_order <- .niche_species_order(rows)
  df <- .niche_long_df(rows, stages)
  df <- .niche_bin_heights(df, "score", .NICHE_DISP_BIN)
  df$species <- factor(.niche_sp_abbr(df$species), levels = .niche_sp_abbr(sp_order))
  df$site    <- factor(df$site)

  omitted <- setdiff(sort(unique(as.character(built$rows$site))),
                     levels(df$site))

  multi_stage <- length(stages) > 1
  facet <- if (multi_stage) facet_grid(stage ~ site) else facet_wrap(~site, nrow = 2)

  p <- ggplot(df, aes(x = species, y = height, fill = score)) +
    geom_tile(width = 0.9, height = .NICHE_DISP_BIN * 1.02) +
    .NICHE_FILL() +
    facet +
    labs(x = NULL, y = "Height above ground (m)",
         title = if (confirmed_only)
           "Modelled climate-niche suitability by canopy height, per site"
         else
           "Climate-niche suitability by canopy height — unidentified morphospecies",
         subtitle = if (length(omitted) > 0)
           sprintf("Not shown (no %s species): %s",
                   if (confirmed_only) "identified" else "morphospecies",
                   paste(omitted, collapse = ", ")) else NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(2.4, "cm"))

  n_site <- nlevels(df$site); n_sp <- nlevels(df$species)
  ncol_f <- if (multi_stage) n_site else ceiling(n_site / 2)
  nrow_f <- 2
  w <- max(8, 1.8 + ncol_f * (0.30 * n_sp + 0.9))
  h <- max(5.5, 1.8 + nrow_f * 3.4)

  tag <- file_tag %||% (if (!confirmed_only) "_morphospecies"
                        else if (multi_stage) "_rescale" else "")
  out_path <- file.path(out_dir, sprintf("niche_suitability_by_site%s.png", tag))
  ggsave(out_path, plot = p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  message("Saved: ", out_path)
  invisible(p)
}

# ── Fig B: the "forest" ─────────────────────────────────────────────────────
# One panel per confirmed species (niche-height order); within a panel a
# horizontal bar per site -- height above ground on x, colour = after-rescale
# suitability, same tiled encoding as Fig A. Every site row is drawn in every
# panel, so a blank row reads directly as "this species is not confirmed at
# that site". This is the "for each species, where in the canopy -- and at
# which sites -- is the climate suitable" view; the stylised-tree-glyph
# version was dropped 2026-08-30 (the shaping obscured more than it
# conveyed), and the per-site orientation on 2026-09-01.
plot_niche_forest <- function(sites = NULL, out_dir = OUTPUT_DIR,
                              processed_dir = PROCESSED_DIR,
                              height_step = 0.4, built = NULL) {
  built <- built %||% .niche_score_rows_all_sites(sites, processed_dir, height_step)
  rows  <- built$rows
  if (is.null(rows) || nrow(rows) == 0) {
    message("Skipping niche forest -- no sites had confirmed species")
    return(invisible(NULL))
  }
  rows <- rows[.niche_is_confirmed(rows$species), , drop = FALSE]
  if (nrow(rows) == 0) {
    message("Skipping niche forest -- no confirmed species")
    return(invisible(NULL))
  }

  rows <- .niche_bin_heights(rows, "after", .NICHE_DISP_BIN, by = c("site", "species"))

  sp_order   <- .niche_species_order(rows)
  # Sites ordered by their modelled canopy height (shortest at the bottom of
  # each panel), so every panel's y-axis reads the same way.
  site_hmax  <- tapply(rows$height, rows$site, max)
  rows$site  <- factor(rows$site, levels = names(sort(site_hmax)))
  rows$species <- factor(.niche_sp_abbr(rows$species), levels = .niche_sp_abbr(sp_order))

  p <- ggplot(rows, aes(x = height, y = site, fill = after)) +
    geom_tile(height = 0.78, width = .NICHE_DISP_BIN * 1.02) +
    .NICHE_FILL() +
    scale_y_discrete(drop = FALSE) +
    facet_wrap(~species, ncol = 4) +
    labs(x = "Height above ground (m)", y = NULL,
         title = "Confirmed-species climate niches, by species and site",
         subtitle = paste("Colour = climate-niche suitability along each site's canopy;",
                          "a blank row = species not confirmed at that site")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "italic"),
          legend.position = "bottom",
          legend.key.width = unit(2.4, "cm"))

  n_sp <- length(sp_order); n_site <- nlevels(rows$site)
  ncol_f <- min(4, n_sp)
  w <- max(10, 1.6 + ncol_f * 2.7)
  h <- max(4.5, 1.4 + ceiling(n_sp / ncol_f) * (0.5 + 0.28 * n_site))
  out_path <- file.path(out_dir, "niche_forest_all_sites.png")
  ggsave(out_path, plot = p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  message("Saved: ", out_path)
  invisible(p)
}

# ── Supp: species confirmed at 2+ sites, cross-site niche comparison ─────────
# The subset where a site-to-site comparison is even meaningful. One panel
# per species; horizontal bars -- y = site, x = height, fill = after-rescale
# suitability -- so each site's niche profile reads as a strip and the sites
# stack for direct comparison.
plot_niche_multisite_species <- function(sites = NULL, out_dir = OUTPUT_DIR,
                                         processed_dir = PROCESSED_DIR,
                                         height_step = 0.4, built = NULL) {
  built <- built %||% .niche_score_rows_all_sites(sites, processed_dir, height_step)
  rows  <- built$rows
  if (is.null(rows) || nrow(rows) == 0) {
    message("Skipping niche multi-site species -- no sites had confirmed species")
    return(invisible(NULL))
  }
  rows <- rows[.niche_is_confirmed(rows$species), , drop = FALSE]
  site_counts <- table(unique(rows[, c("site", "species")])$species)
  multi <- names(site_counts)[site_counts >= 2]
  rows <- rows[rows$species %in% multi, , drop = FALSE]
  if (nrow(rows) == 0) {
    message("Skipping niche multi-site species -- none confirmed at 2+ sites")
    return(invisible(NULL))
  }

  sp_order <- .niche_species_order(rows)
  df <- .niche_long_df(rows, "after")
  df <- .niche_bin_heights(df, "score", .NICHE_DISP_BIN)
  df$species <- factor(.niche_sp_abbr(df$species), levels = .niche_sp_abbr(sp_order))

  p <- ggplot(df, aes(x = height, y = site, fill = score)) +
    geom_tile(height = 0.82, width = .NICHE_DISP_BIN * 1.02) +
    .NICHE_FILL() +
    facet_wrap(~species, scales = "free_y", ncol = 2) +
    labs(x = "Height above ground (m)", y = NULL,
         title = "Species confirmed at 2+ sites — cross-site niche comparison") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(2.4, "cm"),
          strip.text = element_text(face = "italic"))

  n_sp <- length(multi); n_site <- length(unique(df$site))
  ncol_f <- min(2, n_sp)
  w <- max(8, 1.5 + ncol_f * 4.2)
  h <- max(3.5, 1.2 + ceiling(n_sp / ncol_f) * (0.5 + 0.42 * n_site))
  out_path <- file.path(out_dir, "niche_multisite_species.png")
  ggsave(out_path, plot = p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  message("Saved: ", out_path)
  invisible(p)
}

# ── Internal live visualisation ───────────────────────────────────────────────

.plot_live <- function(state, abundanceS, abundanceJ, abundanceA,
                       totalS, totalJ, totalA, t, carCap, sleeptime = 0.2) {
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols <- .species_colors(state$n_species)  # shared_helpers.R
  total <- totalS + totalJ + totalA
  n_sp <- state$n_species
  par(mfrow = c(1, n_sp + 1), mar = c(4, 4, 3, 2))
  for (sp in seq_len(n_sp)) {
    ts_S <- sapply(1:t, function(i) sum(abundanceS[, , , i, sp]))
    ts_J <- sapply(1:t, function(i) sum(abundanceJ[, , , i, sp]))
    ts_A <- sapply(1:t, function(i) sum(abundanceA[, , , i, sp]))
    ts_total <- ts_S + ts_J + ts_A
    plot(ts_total,
      type = "b", col = sp_cols[sp], lwd = 2,
      ylim = c(0, max(ts_total, 1)), xlab = "Year", ylab = "Abundance",
      main = paste0(state$species_ids[sp], " (t=", t, ")"), las = 1
    )
    lines(ts_S, type = "b", col = stage_cols[1], pch = 16, lty = 2)
    lines(ts_J, type = "b", col = stage_cols[2], pch = 17, lty = 2)
    lines(ts_A, type = "b", col = stage_cols[3], pch = 15, lty = 2)
    legend("topleft",
      legend = c("Total", "S", "J", "A"),
      col = c(sp_cols[sp], stage_cols), lty = c(1, 2, 2, 2),
      pch = c(NA, 16, 17, 15), cex = 0.6
    )
  }
  plot(total[1:t],
    type = "b", col = "black", lwd = 2,
    ylim = c(0, max(total, 1)), xlab = "Year", ylab = "Abundance",
    main = paste0("All species (t=", t, ")"), las = 1
  )
  lines(totalS[1:t], type = "b", col = stage_cols[1], pch = 16)
  lines(totalJ[1:t], type = "b", col = stage_cols[2], pch = 17)
  lines(totalA[1:t], type = "b", col = stage_cols[3], pch = 15)
  abline(h = carCap * state$xDim * state$yDim * state$zDim, col = "red", lty = 2)
  legend("topleft",
    legend = c("Total", "S", "J", "A"),
    col = c("black", stage_cols), lty = 1, pch = c(NA, 16, 17, 15), cex = 0.6
  )
  dev.flush()
  Sys.sleep(sleeptime)
}

# ── Post-hoc abundance plot ───────────────────────────────────────────────────

plot_abundance <- function(result, t = NULL, species_specific = TRUE) {
  state <- result$state
  t_max <- if (is.null(t)) length(result$totalabundanceA) else t
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols <- .species_colors(state$n_species)  # shared_helpers.R
  totalS <- result$totalabundanceS
  totalJ <- result$totalabundanceJ
  totalA <- result$totalabundanceA
  total <- totalS + totalJ + totalA
  n_panels <- if (species_specific) state$n_species + 1L else 1L
  par(mfrow = c(1, n_panels), mar = c(4, 4, 3, 2))
  if (species_specific) {
    for (sp in seq_len(state$n_species)) {
      ts_S <- sapply(1:t_max, function(i) sum(result$abundanceS[, , , i, sp]))
      ts_J <- sapply(1:t_max, function(i) sum(result$abundanceJ[, , , i, sp]))
      ts_A <- sapply(1:t_max, function(i) sum(result$abundanceA[, , , i, sp]))
      ts_total <- ts_S + ts_J + ts_A
      plot(ts_total,
        type = "b", col = sp_cols[sp], lwd = 2,
        ylim = c(0, max(ts_total, 1)), xlab = "Year", ylab = "Abundance",
        main = paste0(state$species_ids[sp], " (t=", t_max, ")"), las = 1
      )
      lines(ts_S, type = "b", col = stage_cols[1], pch = 16, lty = 2)
      lines(ts_J, type = "b", col = stage_cols[2], pch = 17, lty = 2)
      lines(ts_A, type = "b", col = stage_cols[3], pch = 15, lty = 2)
      legend("topleft",
        legend = c("Total", "S", "J", "A"),
        col = c(sp_cols[sp], stage_cols), lty = c(1, 2, 2, 2),
        pch = c(NA, 16, 17, 15), cex = 0.6
      )
    }
  }
  plot(total[1:t_max],
    type = "b", col = "black", lwd = 2,
    ylim = c(0, max(total[1:t_max], 1)), xlab = "Year", ylab = "Abundance",
    main = paste0(state$site_name, " — All species (t=", t_max, ")"), las = 1
  )
  lines(totalS[1:t_max], type = "b", col = stage_cols[1], pch = 16)
  lines(totalJ[1:t_max], type = "b", col = stage_cols[2], pch = 17)
  lines(totalA[1:t_max], type = "b", col = stage_cols[3], pch = 15)
  abline(h = state$carCap * state$xDim * state$yDim * state$zDim, col = "red", lty = 2)
  legend("topleft",
    legend = c("Total", "S", "J", "A"),
    col = c("black", stage_cols), lty = 1, pch = c(NA, 16, 17, 15), cex = 0.6
  )
}

# ggplot version of plot_abundance()'s per-species-per-stage breakdown
# (2026-08-24) -- reuses the exact same summing logic
# (sum(result$abundanceS[,,,i,sp]) per timestep/species, see plot_abundance()
# above) but as a proper multi-replicate ggplot figure instead of a base-R
# device plot, so it can show replicate spread the same way
# plot_experiment() (get_colonization.R) already does for OAT sweeps: thin
# per-replicate lines at low alpha, thick mean line on top. Reads its own
# input file (house convention -- every plotting function in this file loads
# what it needs rather than taking an already-loaded object), same
# result$runs/list(result) fallback plot_default_colonization_run() already
# uses for an RDS saved before run_replicated() existed.
plot_species_stage_curves <- function(site_name, exp_tag = "best_combo",
                                      out_dir = OUTPUT_DIR,
                                      processed_dir = PROCESSED_DIR) {
  in_path <- file.path(processed_dir, sprintf("colonization_%s_%s_h0.40.rds", site_name, exp_tag))
  if (!file.exists(in_path)) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no results at ", in_path)
    return(invisible(NULL))
  }
  result <- readRDS(in_path)
  if (is.data.frame(result)) {
    message("Skipping ", site_name, " / ", exp_tag,
            " -- sweep result, has no per-species abundanceS/J/A arrays")
    return(invisible(NULL))
  }
  runs <- if (is.list(result) && !is.null(result$runs)) result$runs else list(result)
  runs <- Filter(Negate(is.null), runs)
  if (length(runs) == 0) {
    message("Skipping ", site_name, " / ", exp_tag, " -- no successful replicates")
    return(invisible(NULL))
  }

  # 2026-08-25: fixed a real crash ("subscript out of bounds") -- different
  # replicates of the SAME site+params run can have DIFFERENT species_ids/
  # n_species (confirmed on Mashpi's best_combo run: 10/8/9 species across
  # 3 reps). A species with very few observations apparently can drop out
  # of a replicate's modeled set depending on that replicate's random
  # train/val split (runcolonization()'s per-species niche availability
  # check). Assuming every replicate shares rep 1's species_ids/n_species
  # (indexing every run's abundanceS[,,,i,sp] the same way) breaks the
  # moment a later replicate has fewer species than rep 1. Fixed by looking
  # up each replicate's OWN species_ids by NAME, skipping a species
  # entirely for whichever replicate(s) didn't model it, rather than
  # assuming a shared, fixed species list/index across replicates.
  all_species <- sort(unique(unlist(lapply(runs, function(r) r$species_ids))))
  sp_cols <- .species_colors(length(all_species))  # shared_helpers.R
  names(sp_cols) <- all_species

  long_df <- do.call(rbind, lapply(seq_along(runs), function(rep_i) {
    run <- runs[[rep_i]]
    t_max <- length(run$totalabundanceA)
    do.call(rbind, lapply(seq_along(run$species_ids), function(sp) {
      ts_S <- sapply(1:t_max, function(i) sum(run$abundanceS[, , , i, sp]))
      ts_J <- sapply(1:t_max, function(i) sum(run$abundanceJ[, , , i, sp]))
      ts_A <- sapply(1:t_max, function(i) sum(run$abundanceA[, , , i, sp]))
      rbind(
        data.frame(rep = rep_i, species = run$species_ids[sp], stage = "Seedling (S)", t = 1:t_max, abundance = ts_S),
        data.frame(rep = rep_i, species = run$species_ids[sp], stage = "Juvenile (J)", t = 1:t_max, abundance = ts_J),
        data.frame(rep = rep_i, species = run$species_ids[sp], stage = "Adult (A)",    t = 1:t_max, abundance = ts_A)
      )
    }))
  }))
  long_df$stage <- factor(long_df$stage, levels = c("Seedling (S)", "Juvenile (J)", "Adult (A)"))

  mean_df <- aggregate(abundance ~ species + stage + t, data = long_df, FUN = mean)

  p <- ggplot(long_df, aes(x = t, y = abundance, colour = species)) +
    geom_line(aes(group = interaction(species, rep)), alpha = 0.25, linewidth = 0.4) +
    geom_line(data = mean_df, linewidth = 1.0) +
    scale_colour_manual(values = sp_cols, name = "Species") +
    facet_wrap(~stage, nrow = 1, scales = "free_y") +
    labs(x = "Year", y = "Abundance",
         title = sprintf("Per-species abundance by life stage — %s", site_name),
         subtitle = sprintf("%s (%d replicate%s); thick = mean, thin = individual replicates",
                            exp_tag, length(runs), if (length(runs) > 1) "s" else "")) +
    theme_minimal(base_size = 11)

  out_path <- file.path(out_dir, sprintf("species_stage_curves_%s_%s.png", site_name, exp_tag))
  ggsave(out_path, plot = p, width = 14, height = 5.5, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# Same per-species/per-stage summing logic as plot_species_stage_curves()
# above, but for ONE species across MULTIPLE sites' best_combo runs --
# color = site instead of species, so it directly compares how the same
# taxon fares under each site's own best-performing parameter combo and
# landscape. 2026-08-25. `sites` defaults to .confirmed_species_sites()
# (shared_helpers.R) filtered to sites with an actual best_combo result on
# disk -- see that helper's header for why it checks the raw Identification
# column rather than FinalID (blank IDs default to "Maxillaria acutifolia").
plot_species_across_sites <- function(species_name, sites = NULL, exp_tag = "best_combo",
                                      out_dir = OUTPUT_DIR, processed_dir = PROCESSED_DIR) {
  requested_sites <- sites %||% .confirmed_species_sites(species_name)
  have_data <- vapply(requested_sites, function(s)
    file.exists(file.path(processed_dir, sprintf("colonization_%s_%s_h0.40.rds", s, exp_tag))),
    logical(1))
  excluded <- requested_sites[!have_data]
  use_sites <- requested_sites[have_data]
  if (length(use_sites) < 2) {
    message("Skipping ", species_name, " -- fewer than 2 sites with a ", exp_tag,
            " result available (requested: ", paste(requested_sites, collapse = ", "), ")")
    return(invisible(NULL))
  }

  long_df <- do.call(rbind, lapply(use_sites, function(site_name) {
    in_path <- file.path(processed_dir, sprintf("colonization_%s_%s_h0.40.rds", site_name, exp_tag))
    result  <- readRDS(in_path)
    runs <- if (is.list(result) && !is.null(result$runs)) result$runs else list(result)
    runs <- Filter(Negate(is.null), runs)
    if (length(runs) == 0) {
      message("Skipping ", site_name, " for ", species_name, " -- no successful replicates")
      return(NULL)
    }
    # 2026-08-25: same per-replicate species_ids fix as
    # plot_species_stage_curves() -- a species can be present in some of a
    # site's replicates and absent from others (see that function's header
    # note), so `sp` must be looked up fresh per replicate rather than once
    # from runs[[1]]$state, and a replicate lacking the species contributes
    # no rows instead of indexing the wrong column or crashing.
    rows_per_rep <- lapply(seq_along(runs), function(rep_i) {
      run <- runs[[rep_i]]
      sp <- match(species_name, run$species_ids)
      if (is.na(sp)) {
        return(NULL)
      }
      t_max <- length(run$totalabundanceA)
      ts_S <- sapply(1:t_max, function(i) sum(run$abundanceS[, , , i, sp]))
      ts_J <- sapply(1:t_max, function(i) sum(run$abundanceJ[, , , i, sp]))
      ts_A <- sapply(1:t_max, function(i) sum(run$abundanceA[, , , i, sp]))
      rbind(
        data.frame(site = site_name, rep = rep_i, stage = "Seedling (S)", t = 1:t_max, abundance = ts_S),
        data.frame(site = site_name, rep = rep_i, stage = "Juvenile (J)", t = 1:t_max, abundance = ts_J),
        data.frame(site = site_name, rep = rep_i, stage = "Adult (A)",    t = 1:t_max, abundance = ts_A)
      )
    })
    rows_per_rep <- Filter(Negate(is.null), rows_per_rep)
    if (length(rows_per_rep) == 0) {
      message("Skipping ", site_name, " for ", species_name,
              " -- not modeled in any replicate (", paste(runs[[1]]$species_ids, collapse = ", "), ")")
      return(NULL)
    }
    do.call(rbind, rows_per_rep)
  }))
  if (is.null(long_df) || length(unique(long_df$site)) < 2) {
    message("Skipping ", species_name, " -- fewer than 2 sites actually modeled this species")
    return(invisible(NULL))
  }
  long_df$stage <- factor(long_df$stage, levels = c("Seedling (S)", "Juvenile (J)", "Adult (A)"))

  site_names <- sort(unique(long_df$site))
  site_cols  <- setNames(.species_colors(length(site_names)), site_names)  # shared_helpers.R

  mean_df <- aggregate(abundance ~ site + stage + t, data = long_df, FUN = mean)

  excluded_note <- if (length(excluded) > 0)
    sprintf(" | Excluded: %s (no %s result, or unconfirmed ID -- see .confirmed_species_sites())",
            paste(excluded, collapse = ", "), exp_tag)
  else ""

  p <- ggplot(long_df, aes(x = t, y = abundance, colour = site)) +
    geom_line(aes(group = interaction(site, rep)), alpha = 0.25, linewidth = 0.4) +
    geom_line(data = mean_df, linewidth = 1.0) +
    scale_colour_manual(values = site_cols, name = "Site") +
    facet_wrap(~stage, nrow = 1, scales = "free_y") +
    labs(x = "Year", y = "Abundance",
         title = species_name,
         subtitle = sprintf("Compared across: %s%s", paste(use_sites, collapse = ", "), excluded_note)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "italic"))

  slug <- gsub("[^A-Za-z0-9]+", "_", species_name)
  out_path <- file.path(out_dir, sprintf("species_across_sites_%s.png", slug))
  ggsave(out_path, plot = p, width = 14, height = 5.5, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

# ── 3D post-hoc visualisation ─────────────────────────────────────────────────

# Translucent point-cloud underlay shared by plot_3d_abundance() and
# plot_3d_abundance_animated(): state$landscape is the full boolean canopy-
# occupancy array already computed for the run, reused here purely as a
# visual backdrop so abundance markers read as embedded in the canopy rather
# than floating in empty space. Subsampled (default cap 20,000 voxels) --
# this is a shape cue, not a faithful full-resolution render, and plotly gets
# slow/heavy (especially the animated HTML export) well before every valid
# voxel is actually needed to convey "there is canopy here."
.canopy_context_trace <- function(state, max_points = 20000, seed = 1) {
  idx <- which(state$landscape, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(NULL)
  }
  if (nrow(idx) > max_points) {
    set.seed(seed)
    idx <- idx[sample.int(nrow(idx), max_points), , drop = FALSE]
  }
  data.frame(x = idx[, 1], y = idx[, 2], z = idx[, 3])
}

plot_3d_abundance <- function(result, t = NULL, show_canopy = TRUE,
                              canopy_opacity = 0.05, canopy_max_points = 20000) {
  state <- result$state
  if (is.null(t)) t <- dim(result$abundanceA)[4]
  sp_cols <- .species_colors(state$n_species)  # shared_helpers.R
  rows <- list()
  for (sp in seq_len(state$n_species)) {
    sp_name <- state$species_ids[sp]
    col <- sp_cols[sp]
    for (stg in list(
      list(arr = result$abundanceS, nm = "S", sym = "circle"),
      list(arr = result$abundanceJ, nm = "J", sym = "diamond"),
      list(arr = result$abundanceA, nm = "A", sym = "square")
    )) {
      idx <- which(stg$arr[, , , t, sp] > 0, arr.ind = TRUE)
      if (nrow(idx) > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          x = idx[, 1], y = idx[, 2], z = idx[, 3],
          species = sp_name, stage = stg$nm, color = col, symbol = stg$sym,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) {
    message("No individuals to plot at t=", t)
    return(invisible(NULL))
  }
  df <- do.call(rbind, rows)
  traces <- split(df, paste0(df$species, "_", df$stage))
  fig <- plotly::plot_ly()
  if (show_canopy) {
    canopy_df <- .canopy_context_trace(state, max_points = canopy_max_points)
    if (!is.null(canopy_df)) {
      fig <- plotly::add_trace(fig,
        data = canopy_df, x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        name = "Canopy", showlegend = TRUE,
        marker = list(
          color = "#6b4423", size = 2,
          opacity = canopy_opacity
        )
      )
    }
  }
  for (tr in traces) {
    fig <- plotly::add_trace(fig,
      data = tr, x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      name = paste0(tr$species[1], " ", tr$stage[1]),
      marker = list(
        symbol = tr$symbol[1], color = tr$color[1],
        size = 6, opacity = 0.85
      )
    )
  }
  fig <- plotly::layout(fig,
    title = paste0("Abundance (t=", t, ")"),
    scene = list(
      xaxis = list(title = "x"), yaxis = list(title = "y"),
      zaxis = list(title = "height tier")
    )
  )
  print(fig)
  invisible(fig)
}

# Same idea as plot_3d_abundance() but across every timestep, using plotly's
# built-in frame/animation support (play button + slider) instead of a
# single static scatter. Saved as a self-contained HTML if out_path is
# given. Not yet run against real output — verify once you have a result
# worth animating (e.g. from a best_case replicate that actually persists).
# The canopy-context trace (show_canopy=TRUE default, added 2026-07-24) in
# particular needs a visual check: mixing an unframed static trace with a
# framed animated one in the same figure is standard plotly behavior, but
# wasn't exercised against a real result in this session -- open the saved
# HTML once and confirm the canopy points stay put while the slider moves.
plot_3d_abundance_animated <- function(result, out_path = NULL, show_canopy = TRUE,
                                       canopy_opacity = 0.05, canopy_max_points = 20000) {
  state <- result$state
  n_t <- dim(result$abundanceA)[4]
  sp_cols <- .species_colors(state$n_species)  # shared_helpers.R

  rows <- list()
  for (t in seq_len(n_t)) {
    for (sp in seq_len(state$n_species)) {
      sp_name <- state$species_ids[sp]
      col <- sp_cols[sp]
      for (stg in list(
        list(arr = result$abundanceS, nm = "S", sym = "circle"),
        list(arr = result$abundanceJ, nm = "J", sym = "diamond"),
        list(arr = result$abundanceA, nm = "A", sym = "square")
      )) {
        idx <- which(stg$arr[, , , t, sp] > 0, arr.ind = TRUE)
        if (nrow(idx) > 0) {
          rows[[length(rows) + 1]] <- data.frame(
            x = idx[, 1], y = idx[, 2], z = idx[, 3], t = t,
            species = sp_name, stage = stg$nm, color = col,
            symbol = stg$sym, trace = paste0(sp_name, " ", stg$nm),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(rows) == 0) {
    message("No individuals to plot across any timestep")
    return(invisible(NULL))
  }
  df <- do.call(rbind, rows)

  # trace -> color lookup (one row per unique trace, in matching order —
  # safer than pairing two independently-deduplicated vectors)
  trace_lu <- df[!duplicated(df$trace), c("trace", "color")]
  colors_named <- setNames(trace_lu$color, trace_lu$trace)

  # Canopy trace is added first, with no `frame` mapping, so it renders as a
  # static backdrop that persists unchanged across every animation frame
  # (plotly supports mixing framed and unframed traces in one figure) --
  # only the abundance trace below actually animates by year.
  fig <- plotly::plot_ly()
  if (show_canopy) {
    canopy_df <- .canopy_context_trace(state, max_points = canopy_max_points)
    if (!is.null(canopy_df)) {
      fig <- plotly::add_trace(fig,
        data = canopy_df, x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        name = "Canopy", showlegend = TRUE,
        marker = list(
          color = "#6b4423", size = 2,
          opacity = canopy_opacity
        )
      )
    }
  }
  fig <- plotly::add_trace(
    fig,
    data = df, x = ~x, y = ~y, z = ~z, frame = ~t, color = ~trace,
    colors = colors_named,
    symbol = ~symbol, symbols = c(circle = "circle", diamond = "diamond", square = "square"),
    type = "scatter3d", mode = "markers",
    marker = list(size = 6, opacity = 0.85)
  )
  fig <- fig |>
    plotly::layout(
      title = "Abundance over time",
      scene = list(
        xaxis = list(title = "x"), yaxis = list(title = "y"),
        zaxis = list(title = "height tier")
      )
    ) |>
    plotly::animation_opts(frame = 400, transition = 200, redraw = TRUE) |>
    plotly::animation_slider(currentvalue = list(prefix = "Year: "))

  if (!is.null(out_path)) {
    # selfcontained=TRUE needs pandoc (not installed on the cluster); FALSE
    # writes a small "<name>_files/" dependency folder alongside the HTML
    # instead -- keep the two together when copying/viewing elsewhere.
    htmlwidgets::saveWidget(fig, out_path, selfcontained = FALSE)
    message("Saved: ", out_path)
  }
  invisible(fig)
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
