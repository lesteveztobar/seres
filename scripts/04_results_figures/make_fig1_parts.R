suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
})

BASE_DIR <- "/home/s38leste_hpc/seres"
IN_DIR   <- file.path(BASE_DIR, "report/figures/figure 1")
OUT_DIR_REP <- file.path(BASE_DIR, "report/figures")
OUT_DIR_MAN <- file.path(BASE_DIR, "manuscript for paper/figures")
TXT_UP   <- 2 / ggplot2::.pt

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
  
  xs_clr <- seq(max(1L, title_clear_x[1]), min(w, title_clear_x[2]))
  ys_clr <- seq(max(1L, title_clear_y[1]), min(h, title_clear_y[2]))
  r[ys_clr, xs_clr] <- 1
  g[ys_clr, xs_clr] <- 1
  b[ys_clr, xs_clr] <- 1
  
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

make_panel_plot <- function(img, title_text, panel_h, pad_x = 0.025,
                            body_ymin = 0.022, body_ymax = panel_h - 0.052,
                            title_size = 5.0) {
  avail_w <- 1.0 - 2 * pad_x
  avail_h <- body_ymax - body_ymin
  
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

# Part 1: A and B
p_part1 <- wrap_plots(panels[1:2], ncol = 1, heights = hs[1:2]) &
  theme(plot.margin = margin(0, 0, 0, 0))
p_part1 <- p_part1 + plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2, unit = "mm"),
                                                   plot.background = element_rect(fill = "white", colour = NA)))

# Part 2: C and D
p_part2 <- wrap_plots(panels[3:4], ncol = 1, heights = hs[3:4]) &
  theme(plot.margin = margin(0, 0, 0, 0))
p_part2 <- p_part2 + plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2, unit = "mm"),
                                                   plot.background = element_rect(fill = "white", colour = NA)))

FIG_W_MM <- 180
H1_MM <- FIG_W_MM * sum(hs[1:2]) + 4
H2_MM <- FIG_W_MM * sum(hs[3:4]) + 4

for (dir in c(OUT_DIR_REP, OUT_DIR_MAN)) {
  if (dir.exists(dir)) {
    ggsave(file.path(dir, "fig1_part1.png"), p_part1, width = FIG_W_MM + 4, height = H1_MM, units = "mm", dpi = 400, bg = "white")
    ggsave(file.path(dir, "fig1_part1.pdf"), p_part1, width = FIG_W_MM + 4, height = H1_MM, units = "mm", device = grDevices::cairo_pdf, bg = "white")
    ggsave(file.path(dir, "fig1_part2.png"), p_part2, width = FIG_W_MM + 4, height = H2_MM, units = "mm", dpi = 400, bg = "white")
    ggsave(file.path(dir, "fig1_part2.pdf"), p_part2, width = FIG_W_MM + 4, height = H2_MM, units = "mm", device = grDevices::cairo_pdf, bg = "white")
    cat("Saved part1 and part2 in", dir, "\n")
  }
}
cat("Done!\n")

