# fig_common.R -- shared setup for the Results figures in
# output/figures/results/. Colours come only from seres_palette.R; the theme
# keeps every non-data element neutral grey.
suppressPackageStartupMessages(library(ggplot2))
source("scripts/99_shared/seres_palette.R")

FIG_DIR <- "output/figures/results"
EXP_DIR <- "output/tables/exports"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
FULL_W <- 180; HALF_W <- 90   # mm

# 2026-10-23: every text size raised by 2 pt (base 8 -> 10). geom_text /
# annotate sizes are in mm, so +2 pt there is + TXT_UP (2 / ggplot2::.pt).
TXT_UP <- 2 / ggplot2::.pt
theme_seres <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
          axis.line = element_line(colour = "grey40", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
          axis.text = element_text(colour = "grey20"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.caption = element_text(colour = "grey35", hjust = 0, size = base_size - 1, lineheight = 1.1),
          plot.caption.position = "plot",
          plot.background = element_rect(fill = "white", colour = NA),
          legend.key.size = unit(3.5, "mm"))
}

# png + pdf, never over an existing file outside FIG_DIR
save_fig <- function(name, p, width = FULL_W, height) {
  for (dev in c("png", "pdf")) {
    f <- file.path(FIG_DIR, paste0(name, ".", dev))
    ggsave(f, plot = p, width = width, height = height, units = "mm", dpi = 300,
           device = if (dev == "pdf") grDevices::cairo_pdf else dev, bg = "white")
    cat("Saved:", f, "\n")
  }
}

# ggplot stroke = 1 is ~0.53 mm, so the 0.3 mm dark stroke is stroke = 0.57.
STROKE_W <- 0.57
# Point fill/stroke obeying the fill-only rule: FILL_ONLY sites and every
# filled point get a dark stroke; open points of other sites are stroked in
# their own colour.
pt_fill   <- function(site, filled) ifelse(filled, SITE_COL[site], "white")
pt_stroke <- function(site, filled) ifelse(filled | site %in% FILL_ONLY, STROKE_DARK, SITE_COL[site])

# Footprint mean elevation (the manuscript's figures); La Elenita has no
# footprint and falls back to its observation-based elevation. The order is
# the same under either source.
site_elevation <- function() {
  e <- read.csv(file.path(EXP_DIR, "site_elevation.csv"))
  stopifnot(identical(order(e$elevation_m), order(ifelse(is.na(e$footprint_elev_mean_m), e$elevation_m, e$footprint_elev_mean_m))))
  setNames(ifelse(is.na(e$footprint_elev_mean_m), e$elevation_m, e$footprint_elev_mean_m), e$site)
}

# ── Taxon names (2026-10-03) ────────────────────────────────────────────────
# Binomials fully italic; morphospecies ("Maxillaria sp. B", "Ida sp.",
# "Maxillaria sp1") italic genus with an upright "sp." part. Returns plotmath
# strings: use with labeller = label_parsed in facets and
# labels = function(x) parse(text = taxon_expr(x)) in scales.
taxon_expr <- function(x, abbreviate = FALSE) {
  vapply(x, function(s) {
    m <- regmatches(s, regexec("^([A-Z][a-z]+)\\s+(sp\\.?\\s*.*|sp[0-9]+.*)$", s))[[1]]
    if (length(m)) {
      g <- if (abbreviate) paste0(substr(m[2], 1, 1), ".") else m[2]
      sprintf("italic('%s')~'%s'", g, m[3])
    } else {
      if (abbreviate) s <- sub("^([A-Z])[a-z]+\\s+", "\\1. ", s)
      sprintf("italic('%s')", s)
    }
  }, character(1), USE.NAMES = FALSE)
}
taxon_labels <- function(abbreviate = FALSE) function(x) parse(text = taxon_expr(x, abbreviate))

# Yanayacu is the eastern-slope outgroup; the other sites are the western
# (Choco Andino) slope.
OUTGROUP <- "Yanayacu"
slope_group <- function(site) factor(ifelse(site == OUTGROUP, "Eastern-slope outgroup", "Western slope (Chocó Andino)"),
                                     levels = c("Western slope (Chocó Andino)", "Eastern-slope outgroup"))
# Short facet-strip labels for the slope grouping (right-hand strips).
slope_short <- ggplot2::labeller(group = c("Western slope (Chocó Andino)" = "Western\nslope", "Eastern-slope outgroup" = "Outgroup"))
