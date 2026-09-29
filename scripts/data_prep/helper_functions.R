# =============================================================================
# Helper functions
# =============================================================================
# Column spec shared by read_csv() calls throughout — keeps types consistent
# when binding rows across sites and prevents pictures from coercing to logical.
COMBINED_COL_TYPES <- cols(
  Source         = col_character(),
  Area_or_Site   = col_character(),
  lat            = col_double(),
  lon            = col_double(),
  Elevation_m    = col_character(),   # may be a range string e.g. "850-890"
  FieldID        = col_character(),
  Abundance      = col_double(),
  Height_m       = col_double(),
  CanopyHeight_m = col_double(),
  pictures       = col_character(),
  note           = col_character(),
  AI_ID          = col_character(),
  FinalID        = col_character(),
  Genus          = col_character(),
  species        = col_character()
)
