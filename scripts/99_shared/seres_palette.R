# SERES canonical palette. Verified colourblind-safe: worst pairwise CIE76
# separation is dE 18.3 (SITE_COL mid pink vs magenta under deuteranopia),
# all other pairs >= 21.4, none below the dE 10 collision threshold.
SERES_PAL <- c(
  pale_pink = "#FFD6E5",   # contrast vs white  1.31:1  FILL ONLY
  mid_pink  = "#CC79A7",   #                    3.06:1  marginal
  magenta   = "#AA3377",   #                    6.09:1  safe
  dark_plum = "#5A0F3A",   #                   13.43:1  safe
  cyan      = "#33BBEE",   #                    2.21:1  FILL ONLY
  olive     = "#999933",   #                    3.02:1  marginal
  black     = "#000000"    #                   21.00:1  safe
)

# Sites in elevation order. The three Mindo-cluster sites take the pink family
# so that the one locality reads as one group, darkening with elevation.
SITE_LEVELS <- c("Mashpi", "Saloya", "Maquipucuna", "MindoTarabita",
                 "MindoMirador", "Yanayacu", "LaElenita")
SITE_COL <- c(
  Mashpi        = SERES_PAL[["cyan"]],       #  935 m
  Saloya        = SERES_PAL[["mid_pink"]],   # 1309 m  Mindo cluster
  Maquipucuna   = SERES_PAL[["olive"]],      # 1325 m
  MindoTarabita = SERES_PAL[["magenta"]],    # 1537 m  Mindo cluster
  MindoMirador  = SERES_PAL[["dark_plum"]],  # 1697 m  Mindo cluster
  Yanayacu      = SERES_PAL[["black"]],      # 2202 m
  LaElenita     = SERES_PAL[["pale_pink"]]   # climate analyses only
)

# Colours that fail the WCAG 3:1 non-text contrast threshold against white and
# must never be used for a thin line or a small unstroked marker.
FILL_ONLY <- c("Mashpi", "LaElenita")
STROKE_DARK <- SERES_PAL[["dark_plum"]]

scale_colour_site <- function(...) ggplot2::scale_colour_manual(
  values = SITE_COL, limits = SITE_LEVELS, drop = TRUE, name = "Site", ...)
scale_fill_site <- function(...) ggplot2::scale_fill_manual(
  values = SITE_COL, limits = SITE_LEVELS, drop = TRUE, name = "Site", ...)
