# fig_fig1_landscape.R -- Build the main 4-panel SERES overview figure
# (report/figures/fig1_landscape.png) from the 4 panel figures in
# report/figures/figure 1/fig1_{A,B,C,D}.png using ggplot2 and patchwork.
#
# Layout:
#   Top row (Initialization):
#     Panel A (left)  : Initialising the landscape
#     Panel B (right) : Initialising the forest
#   Bottom row (Model + simulation outputs):
#     Panel C (left)  : Running the colonisation model
#     Panel D (right) : Generating 50-year simulation results

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
})

#
# 2026-10-23: second layout added. Run with an argument:
#   Rscript scripts/04_results_figures/fig_fig1_landscape.R grid      (default; 2x2, report/figures/)
#   Rscript scripts/04_results_figures/fig_fig1_landscape.R vertical  (A/B/C/D stacked in one column,
#                                                                     "manuscript for paper/figures/")
# Panel titles are 2 pt larger than before in both layouts.

BASE_DIR <- if (dir.exists("report/figures/figure 1")) "." else "/home/s38leste_hpc/seres"
IN_DIR   <- file.path(BASE_DIR, "report/figures/figure 1")
LAYOUT   <- match.arg(commandArgs(trailingOnly = TRUE)[1] %||% "grid", c("grid", "vertical"))
OUT_DIR  <- if (LAYOUT == "vertical") file.path(BASE_DIR, "manuscript for paper/figures") else file.path(BASE_DIR, "report/figures")
OUT_PNG  <- file.path(OUT_DIR, "fig1_landscape.png")
OUT_PDF  <- file.path(OUT_DIR, "fig1_landscape.pdf")
TXT_UP   <- 2 / ggplot2::.pt   # +2 pt in annotate() size units (mm)

# Load PNG, composite over pure white, blank the baked-in top-center title box,
# and crop tightly to the remaining non-white diagram bounding box.
load_and_crop_panel <- function(filename, title_clear_x, title_clear_y) {
  path <- file.path(IN_DIR, filename)
  stopifnot(file.exists(path))
  raw <- readPNG(path)
  h <- dim(raw)[1]
  w <- dim(raw)[2]
  
  if ( length(dim(raw)) == 3 && dim(raw)[3] == 4 ) {
    alpha <- raw[, , 4]
    r <- raw[, , 1] * alpha + (1 - alpha)
    g <- raw[, , 2] * alpha + (1 - alpha)
    b <- raw[, , 3] * alpha + (1 - alpha)
  } else if ( length(dim(raw)) == 3 && dim(raw)[3] == 3 ) {
    r <- raw[, , 1]
    g <- raw[, , 2]
    b <- raw[, , 3]
  } else {
    stop("Unexpected PNG channel dimensions for ", filename)
  }
  
  # Clear baked-in title region so ggplot2 renders unified vector titles
  xs_clr <- seq(max(1L, title_clear_x[1]), min(w, title_clear_x[2]))
  ys_clr <- seq(max(1L, title_clear_y[1]), min(h, title_clear_y[2]))
  r[ys_clr, xs_clr] <- 1
  g[ys_clr, xs_clr] <- 1
  b[ys_clr, xs_clr] <- 1
  
  # Find tight bounding box of remaining diagram content
  non_white <- (r < 0.98) | (g < 0.98) | (b < 0.98)
  rows <- which(rowSums(non_white) > 0)
  cols <- which(colSums(non_white) > 0)
  
  r_crop <- r[min(rows):max(rows), min(cols):max(cols), drop = FALSE]
  g_crop <- g[min(rows):max(rows), min(cols):max(cols), drop = FALSE]
  b_crop <- b[min(rows):max(rows), min(cols):max(cols), drop = FALSE]
  
  out <- array(1, dim = c(nrow(r_crop), ncol(r_crop), 3))
  out[, , 1] <- r_crop
  out[, , 2] <- g_crop
  out[, , 3] <- b_crop
  
  list(
    raster = out,
    w = ncol(r_crop),
    h = nrow(r_crop),
    aspect = ncol(r_crop) / nrow(r_crop)
  )
}

img_A <- load_and_crop_panel("fig1_A.png", title_clear_x = c(580, 980),  title_clear_y = c(1, 50))
img_B <- load_and_crop_panel("fig1_B.png", title_clear_x = c(640, 980),  title_clear_y = c(1, 50))
img_C <- load_and_crop_panel("fig1_C.png", title_clear_x = c(560, 1040), title_clear_y = c(1, 50))
img_D <- load_and_crop_panel("fig1_D.png", title_clear_x = c(720, 1270), title_clear_y = c(1, 55))

cat(sprintf("Cropped diagram aspects (W/H): A=%.3f, B=%.3f, C=%.3f, D=%.3f\n",
            img_A$aspect, img_B$aspect, img_C$aspect, img_D$aspect))

# Build a single ggplot panel with a crisp border, unified header, and
# aspect-preserved raster placement.
make_panel_plot <- function(img, title_text, panel_h, pad_x = 0.025,
                            body_ymin = 0.022, body_ymax = panel_h - 0.052,
                            title_size = 5.0) {
  avail_w <- 1.0 - 2 * pad_x
  avail_h <- body_ymax - body_ymin
  
  # Fit while strictly preserving 1:1 native pixel aspect ratio
  if (avail_w / avail_h > img$aspect) {
    draw_h <- avail_h
    draw_w <- draw_h * img$aspect
  } else {
    draw_w <- avail_w
    draw_h <- draw_w / img$aspect
  }
  
  x_mid <- 0.5
  y_mid <- 0.5 * (body_ymin + body_ymax)
  xmin <- x_mid - 0.5 * draw_w
  xmax <- x_mid + 0.5 * draw_w
  ymin <- y_mid - 0.5 * draw_h
  ymax <- y_mid + 0.5 * draw_h
  
  ggplot() +
    annotation_raster(img$raster, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, interpolate = TRUE) +
    annotate(
      "text",
      x = 0.5, y = panel_h - 0.026,
      label = title_text,
      fontface = "bold", size = title_size, colour = "#000000",
      hjust = 0.5, vjust = 0.5
    ) +
     geom_rect(
      aes(xmin = 0, xmax = 1, ymin = 0, ymax = panel_h),
      fill = NA, colour = "#111111", linewidth = 0.8
    ) +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, panel_h), expand = FALSE, clip = "on") +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

TITLE_SZ <- 4.15 + TXT_UP

# ── Vertical layout (manuscript): A, B, C, D stacked, each at full width ─────
# Each panel's height follows its own drawing (full width minus padding, at
# the drawing's native aspect) plus a title band, so no drawing is shrunk to
# fit a shared row height. Printed at 180 mm wide; at 500 dpi every source
# pixel of the drawings (1583-2059 px wide) is kept.
if (LAYOUT == "vertical") {
  PAD_X <- 0.020; BODY_YMIN <- 0.020; TITLE_BAND <- 0.058
  panel_h_for <- function(img) BODY_YMIN + (1 - 2 * PAD_X) / img$aspect + TITLE_BAND
  specs <- list(
    list(img = img_A, title = "A: Initialising the landscape"),
    list(img = img_B, title = "B: Initialising the forest"),
    list(img = img_C, title = "C: Running the colonisation model"),
    list(img = img_D, title = "D: Generating 50-year simulation results")
  )
  hs <- vapply(specs, function(s) panel_h_for(s$img), numeric(1))
  panels <- Map(function(s, h) make_panel_plot(s$img, title_text = s$title, panel_h = h, pad_x = PAD_X,
                                               body_ymin = BODY_YMIN, body_ymax = h - TITLE_BAND,
                                               title_size = TITLE_SZ), specs, hs)
  p_vert <- wrap_plots(panels, ncol = 1, heights = hs) &
    theme(plot.margin = margin(0, 0, 0, 0))
  p_vert <- p_vert + plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2, unit = "mm"),
                                                   plot.background = element_rect(fill = "white", colour = NA)))
  FIG_W_MM <- 180
  FIG_H_MM <- FIG_W_MM * sum(hs) + 4   # + the 2 mm outer margin top and bottom
  ggsave(OUT_PNG, p_vert, width = FIG_W_MM + 4, height = FIG_H_MM, units = "mm", dpi = 500, bg = "white")
  cat("Saved:", OUT_PNG, sprintf("(%.0f x %.1f mm @ 500 dpi; panel heights %s)\n", FIG_W_MM + 4, FIG_H_MM,
                                 paste(sprintf("%.3f", hs), collapse = ", ")))
  ggsave(OUT_PDF, p_vert, width = FIG_W_MM + 4, height = FIG_H_MM, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  cat("Saved:", OUT_PDF, "\n")
  quit(save = "no")
}

# Row heights in panel-width units (1.0 = one panel width)
H_TOP <- 0.445
H_BOT <- 0.505

pA <- make_panel_plot(
  img_A,
  title_text = "A: Initialising the landscape",
  panel_h    = H_TOP,
  pad_x      = 0.025,
  body_ymin  = 0.020,
  body_ymax  = H_TOP - 0.056,
  title_size = TITLE_SZ
)

pB <- make_panel_plot(
  img_B,
  title_text = "B: Initialising the forest",
  panel_h    = H_TOP,
  pad_x      = 0.022,
  body_ymin  = 0.020,
  body_ymax  = H_TOP - 0.054,
  title_size = TITLE_SZ
)

# In Panel C, input 7 on the left is tall while the center under the title has a pocket,
# so body_ymax can extend slightly higher without touching the centered title.
pC <- make_panel_plot(
  img_C,
  title_text = "C: Running the colonisation model",
  panel_h    = H_BOT,
  pad_x      = 0.020,
  body_ymin  = 0.026,
  body_ymax  = H_BOT - 0.044,
  title_size = TITLE_SZ
)

pD <- make_panel_plot(
  img_D,
  title_text = "D: Generating 50-year simulation results",
  panel_h    = H_BOT,
  pad_x      = 0.020,
  body_ymin  = 0.020,
  body_ymax  = H_BOT - 0.056,
  title_size = TITLE_SZ
)

# Combine into a 2x2 grid: A and B on top, C and D below
p_combined <- ((pA | pB) / (pC | pD)) +
  plot_layout(
    widths  = c(1, 1),
    heights = c(H_TOP, H_BOT)
  ) &
  theme(plot.margin = margin(0, 0, 0, 0))

p_combined <- p_combined +
  plot_annotation(
    theme = theme(
      plot.margin = margin(2, 2, 2, 2, unit = "mm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

# Total aspect ratio: width = 2.0 units, height = H_TOP + H_BOT = 0.920 units
FIG_W_MM <- 360
FIG_H_MM <- FIG_W_MM * ((H_TOP + H_BOT) / 2.0)

ggsave(
  filename = OUT_PNG,
  plot     = p_combined,
  width    = FIG_W_MM,
  height   = FIG_H_MM,
  units    = "mm",
  dpi      = 400,
  bg       = "white"
)
cat("Saved:", OUT_PNG, sprintf("(%d x %.1f mm @ 400 dpi)\n", FIG_W_MM, FIG_H_MM))

ggsave(
  filename = OUT_PDF,
  plot     = p_combined,
  width    = FIG_W_MM,
  height   = FIG_H_MM,
  units    = "mm",
  device   = grDevices::cairo_pdf,
  bg       = "white"
)
cat("Saved:", OUT_PDF, "\n")
