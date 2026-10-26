# fig_site_map.R -- Results figure 2.3: all seven field sites with the ERA5
# 0.25 deg cells that force them. Replaces output/figures/sites/site_map.png
# (plot_site_map(), plot_functions.R), which omits Saloya and predates the
# palette; that file is left in place.
#   (a) region: every ERA5 cell in view, the forcing cells outlined, Ecuador inset
#   (b) the shared forcing cell: five sites, Mindo-cluster hull, DTM hillshade
# Inputs:
#   site coordinates   data/csv/combinedv6.csv (the runs' CANOPY_OBS_CSV):
#                      mean lon/lat of every record with coordinates, per site
#   forcing cells      output/tables/audits/home_cell_validity_lookup.csv
#                      (home_cell_validity_lookup.R, 2026-09-08: the 0.25 deg cell
#                      containing each site footprint; hour-1 validity checked)
#   terrain            data/raw/<site>/dtm/dtm.tif (already on disk; nothing downloaded)
#   Rscript scripts/04_results_figures/fig_site_map.R
source("scripts/04_results_figures/fig_common.R")
suppressPackageStartupMessages({ library(sf); library(patchwork); library(ggrepel); library(ggspatial) })
sf_use_s2(FALSE)

# ── sites ────────────────────────────────────────────────────────────────────
obs <- read.csv("data/csv/combinedv6.csv", stringsAsFactors = FALSE)
obs$lon <- suppressWarnings(as.numeric(obs$lon)); obs$lat <- suppressWarnings(as.numeric(obs$lat))
obs <- obs[is.finite(obs$lon) & is.finite(obs$lat), ]
sites <- do.call(rbind, lapply(split(obs, obs$Area_or_Site), function(d)
  data.frame(site = d$Area_or_Site[1], lon = mean(d$lon), lat = mean(d$lat), n_records = nrow(d))))
sites <- sites[match(intersect(SITE_LEVELS, sites$site), sites$site), ]
cat("Site coordinates (data/csv/combinedv6.csv, mean of records with lon/lat):\n"); print(sites, row.names = FALSE)
stopifnot(setequal(sites$site, SITE_LEVELS))
sites$fill <- SITE_COL[sites$site]

# ── ERA5 cells ───────────────────────────────────────────────────────────────
CELLS_CSV <- "output/tables/audits/home_cell_validity_lookup.csv"
cells <- read.csv(CELLS_CSV)
print(cells)
stopifnot(setequal(cells$site, SITE_LEVELS), !anyNA(cells$home_cell_lon), all(cells$hour1_valid))
fc <- unique(cells[, c("home_cell_lon", "home_cell_lat")])
names(fc) <- c("era5_lon", "era5_lat")
fc$members <- mapply(function(x, y) sum(cells$home_cell_lon == x & cells$home_cell_lat == y), fc$era5_lon, fc$era5_lat)
fc <- fc[order(-fc$members), ]
fc$id <- LETTERS[seq_len(nrow(fc))]
cell_poly <- function(lon, lat, h = 0.125)
  st_polygon(list(rbind(c(lon - h, lat - h), c(lon + h, lat - h), c(lon + h, lat + h), c(lon - h, lat + h), c(lon - h, lat - h))))
fc_sf <- st_sf(fc, geometry = st_sfc(mapply(cell_poly, fc$era5_lon, fc$era5_lat, SIMPLIFY = FALSE), crs = 4326))
fc_sf$lab <- sprintf("Cell %s: %d site%s", fc$id, fc$members, ifelse(fc$members > 1, "s", ""))

bb <- st_bbox(fc_sf)
xlim_a <- c(min(bb["xmin"], sites$lon) - 0.05, max(bb["xmax"], sites$lon) + 0.05)
ylim_a <- c(min(bb["ymin"], sites$lat) - 0.05, max(bb["ymax"], sites$lat) + 0.05)
gx <- seq(floor(xlim_a[1] * 4) / 4, ceiling(xlim_a[2] * 4) / 4, 0.25)
gy <- seq(floor(ylim_a[1] * 4) / 4, ceiling(ylim_a[2] * 4) / 4, 0.25)
grid_pts <- expand.grid(lon = gx, lat = gy)
grid_sf <- st_sf(grid_pts, geometry = st_sfc(mapply(cell_poly, grid_pts$lon, grid_pts$lat, SIMPLIFY = FALSE), crs = 4326))

# ── Mindo cluster hull (the three pink sites, 1.5 km buffer) ─────────────────
mindo <- c("Saloya", "MindoTarabita", "MindoMirador")
mindo_pts <- st_as_sf(obs[obs$Area_or_Site %in% mindo, ], coords = c("lon", "lat"), crs = 4326)
hull <- st_transform(st_buffer(st_convex_hull(st_union(st_transform(mindo_pts, 32717))), 1500), 4326)
d_km <- as.matrix(dist(st_coordinates(st_transform(st_as_sf(sites[sites$site %in% mindo, ], coords = c("lon", "lat"), crs = 4326), 32717)))) / 1000
dimnames(d_km) <- list(sites$site[sites$site %in% mindo], sites$site[sites$site %in% mindo])
cat("Mindo cluster centroid distances (km):\n"); print(round(d_km, 1))

# ── hillshade for panel b from the per-site DTMs ─────────────────────────────
shared <- fc_sf[1, ]
xlim_b <- c(st_bbox(shared)["xmin"] - 0.02, st_bbox(shared)["xmax"] + 0.03)
ylim_b <- c(st_bbox(shared)["ymin"] - 0.02, st_bbox(shared)["ymax"] + 0.03)
hs_df <- NULL
dtm_f <- file.path("data/raw", SITE_LEVELS, "dtm", "dtm.tif"); dtm_f <- dtm_f[file.exists(dtm_f)]
cat("DTMs found:", length(dtm_f), "\n")
if (length(dtm_f) && requireNamespace("terra", quietly = TRUE)) {
  hs_df <- tryCatch({
    r <- lapply(dtm_f, terra::rast)
    r <- lapply(r, function(x) if (!terra::same.crs(x, r[[1]])) terra::project(x, r[[1]]) else x)
    # DTM resolutions differ: project each onto one lon/lat grid, then average
    tmpl <- terra::rast(terra::ext(xlim_b[1], xlim_b[2], ylim_b[1], ylim_b[2]), resolution = 0.0008, crs = "EPSG:4326")
    m <- terra::app(terra::rast(lapply(r, function(x) terra::project(x, tmpl))), mean, na.rm = TRUE)
    sl <- terra::terrain(m, "slope", unit = "radians"); as <- terra::terrain(m, "aspect", unit = "radians")
    hs <- terra::shade(sl, as, angle = 40, direction = 315)
    df <- as.data.frame(hs, xy = TRUE, na.rm = TRUE); names(df)[3] <- "shade"
    cat(sprintf("Hillshade: %d cells, resolution %.5f deg\n", nrow(df), terra::res(hs)[1]))
    df
  }, error = function(e) { cat("Hillshade skipped:", conditionMessage(e), "\n"); NULL })
}

# ── panels ───────────────────────────────────────────────────────────────────
site_pts <- function(d, size = 3)
  geom_point(data = d, aes(lon, lat, fill = fill), shape = 21, colour = STROKE_DARK, stroke = STROKE_W, size = size)

pa <- ggplot() +
  geom_sf(data = grid_sf, fill = NA, colour = "grey80", linewidth = 0.25) +
  geom_sf(data = fc_sf, fill = "grey93", colour = "grey30", linewidth = 0.5) +
  geom_text(data = grid_pts, aes(lon, lat - 0.095, label = sprintf("%.2f, %.2f", lat, lon)),
            size = 1.5 + TXT_UP, colour = "grey60") +
  geom_text(data = st_drop_geometry(fc_sf), aes(era5_lon, era5_lat + 0.095, label = lab), size = 2.1 + TXT_UP, colour = "grey15", fontface = "bold") +
  geom_sf(data = hull, fill = "grey85", colour = "grey35", linewidth = 0.3, alpha = 0.6) +
  site_pts(sites, 2.6) +
  geom_text_repel(data = sites[!sites$site %in% c(mindo, "LaElenita", "Maquipucuna"), ],
                  aes(lon, lat, label = site), size = 2.5 + TXT_UP, colour = "grey10", box.padding = 0.5, seed = 1) +
  annotate("text", x = st_bbox(shared)["xmax"] + 0.02, y = st_bbox(shared)["ymin"] - 0.03,
           label = "five sites: see (b)", fontface = "italic", size = 2.2 + TXT_UP, colour = "grey25", hjust = 0) +
  scale_fill_identity() +
  coord_sf(xlim = xlim_a, ylim = ylim_a, expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.55 + 2 / 12, height = unit(0.12, "cm"),
                   bar_cols = c("grey30", "white"), line_col = "grey30", text_col = "grey20") +
  labs(title = "a", x = NULL, y = NULL) +
  theme_seres(9) + theme(panel.grid.major = element_blank(), axis.text = element_text(size = 7.5))

ecu <- rnaturalearth::ne_countries(country = "Ecuador", scale = "medium", returnclass = "sf")
inset <- ggplot() +
  geom_sf(data = ecu, fill = "grey92", colour = "grey45", linewidth = 0.2) +
  annotate("rect", xmin = xlim_a[1], xmax = xlim_a[2], ymin = ylim_a[1], ymax = ylim_a[2],
           fill = NA, colour = "grey10", linewidth = 0.4) +
  coord_sf(xlim = c(-81.2, -75.1), ylim = c(-5.1, 1.6), expand = FALSE) +
  theme_void() + theme(panel.background = element_rect(fill = "white", colour = "grey50", linewidth = 0.3))
pa <- pa + inset_element(inset, left = 0.70, bottom = 0.66, right = 0.99, top = 0.99, align_to = "panel")

sites_b <- sites[sites$lon >= xlim_b[1] & sites$lon <= xlim_b[2] & sites$lat >= ylim_b[1] & sites$lat <= ylim_b[2], ]
cat("Panel b sites:", paste(sites_b$site, collapse = ", "), "\n")
hull_c <- st_coordinates(st_centroid(hull))
pb <- ggplot() +
  { if (!is.null(hs_df)) geom_raster(data = hs_df, aes(x, y, alpha = 1 - shade), fill = "grey20") } +
  scale_alpha_continuous(range = c(0, 0.45), guide = "none") +
  geom_sf(data = shared, fill = NA, colour = "grey30", linewidth = 0.5) +
  geom_sf(data = hull, fill = "grey90", colour = "grey35", linewidth = 0.35, alpha = 0.5) +
  annotate("text", x = st_bbox(hull)["xmin"], y = st_bbox(hull)["ymin"] - 0.008, label = "Mindo cluster: one locality",
           size = 2.3 + TXT_UP, colour = "grey20", hjust = 0, fontface = "italic") +
  site_pts(sites_b, 3.2) +
  geom_text_repel(data = sites_b, aes(lon, lat, label = site), size = 2.6 + TXT_UP, colour = "grey10",
                  box.padding = 0.6, point.padding = 0.3, min.segment.length = 0, segment.colour = "grey40",
                  segment.size = 0.25, seed = 2) +
  annotate("text", x = xlim_b[1] + 0.005, y = ylim_b[2] - 0.008, hjust = 0, vjust = 1, size = 2.1 + TXT_UP, colour = "grey20",
           label = sprintf("ERA5 cell %s: one forcing series", shared$id)) +
  scale_fill_identity() +
  coord_sf(xlim = xlim_b, ylim = ylim_b, expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.3, text_cex = 0.55 + 2 / 12, height = unit(0.12, "cm"),
                   bar_cols = c("grey30", "white"), line_col = "grey30", text_col = "grey20") +
  labs(title = "b", x = NULL, y = NULL) +
  theme_seres(9) + theme(panel.grid.major = element_blank(), axis.text = element_text(size = 7.5))

cap <- paste0("Site positions: mean of observation coordinates (combinedv6.csv). Grey squares: ERA5 0.25° cells ",
              "(centre lat, lon); outlined: the home cell of at least one site (home_cell_validity_lookup.csv). ",
              if (!is.null(hs_df)) "Hillshade in (b) from the per-site DTMs." else "No terrain underlay.")
p <- (pa | pb) + plot_layout(widths = c(1.15, 1))  # caption text lives in the LaTeX caption (printed below)
save_fig("fig_site_map", p, width = FULL_W, height = 100)
cat(sprintf("Plotted %d sites in (a), %d in (b); dropped 0\n", nrow(sites), nrow(sites_b)))
cat("Caption text for LaTeX:", cap, "\n")
