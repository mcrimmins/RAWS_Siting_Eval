# ============================================================
# RAWS Exposure from Categorical EVH (LANDFIRE LC23_EVH_240)
# - Keep EVH categorical (factor)
# - Extract per-station buffers, then map classes -> height & metadata
# - Compute "clear area" metrics and X-factor per clearing size
#
# Output:
#   data/RAWS_EVH_summary.csv         (per station × clearing size)
#   data/RAWS_dropped_outside_CONUS.csv (stations outside EVH extent)
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readxl)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(rlang)
})

# ---------------- CONFIG ----------------
# Paths
evh_tif        <- "data/LF2023_EVH_240_CONUS/Tif/LC23_EVH_240.tif"  # categorical EVH with RAT
station_xlsx   <- "data/FEMS_3.0_RAWS_Master_Station_List_and_Metadata.xlsx"
station_sheet  <- NULL  # NULL = first sheet; else sheet name or index

# Column names in your Excel (we’ll normalize safely below)
col_station_id_candidates <- c("station_id", "station_name", "wrcc_id", "nesdis_id")
col_lat_candidates        <- c("lat", "latitude", "LAT", "Latitude")
col_lon_candidates        <- c("lon", "longitude", "LON", "Longitude")

# “Clear” area definition (tune to taste)
clear_height_m          <- 1   # cells <= this are considered open/clear
clear_pct_required      <- 0.95   # require ≥ 95% of cells clear to flag area as clear
treat_nonveg_as_zero    <- TRUE   # non-veg (water, roads, crops, barren, sparse canopy) -> 0 m; FALSE -> NA

# H used in X-factor (choose: "p95", "max", "mean", "median")
H_stat                  <- "p95"

# If TRUE, compute H from a ring outside the clearing (r .. r + ring_width_m)
use_ring_for_H          <- FALSE
ring_width_m            <- 100

# Clearing sizes (acres) to evaluate
clearing_acres <- c(0.1, 2, 10, 40, 100)

# Outputs
out_csv_results <- "data/RAWS_EVH_summary.csv"
out_csv_dropped <- "data/RAWS_dropped_outside_CONUS.csv"

# Terra performance
terraOptions(progress = 1, memfrac = 0.6)
# terraOptions(tempdir = "/fast/ssd/tmp")  # optional: fast scratch

message("terra version: ", as.character(utils::packageVersion("terra")))

# ---------------- LOAD EVH (categorical) ----------------
stopifnot(file.exists(evh_tif))
evh_cat <- terra::rast(evh_tif)
stopifnot(terra::is.factor(evh_cat))

# Ensure meters CRS (EPSG:5070) for correct buffer sizes/areas
if (!grepl("5070", terra::crs(evh_cat))) {
  message("Reprojecting EVH to EPSG:5070 (near-neighbor)...")
  evh_cat <- terra::project(evh_cat, "EPSG:5070", method = "near")
} else {
  message("EVH already in EPSG:5070.")
}

# ---------------- RAT → height + lifeform lookup (row-aligned) ----------------
rat <- terra::levels(evh_cat)[[1]]  # expects columns: Value, CLASSNAMES
stopifnot(all(c("Value","CLASSNAMES") %in% names(rat)))

rat <- rat %>%
  mutate(
    num = suppressWarnings(readr::parse_number(CLASSNAMES)),
    lifeform = case_when(
      str_starts(CLASSNAMES, "Tree Height")  ~ "Tree",
      str_starts(CLASSNAMES, "Shrub Height") ~ "Shrub",
      str_starts(CLASSNAMES, "Herb Height")  ~ "Herb",
      TRUE                                   ~ "NonVeg"
    ),
    height_m = case_when(
      CLASSNAMES == "Fill-NoData" ~ NA_real_,
      lifeform %in% c("Tree","Shrub","Herb") ~ num,
      CLASSNAMES %in% c(
        "Open Water","Snow/Ice","Developed-Roads","Barren",
        "Cultivated Crops","Sparse Vegetation Canopy",
        # common NASS labels within EVH
        "NASS-Row Crop","NASS-Row Crop-Close Grown Crop","NASS-Close Grown Crop","NASS-Wheat",
        "Quarries-Strip Mines-Gravel Pits-Well and Wind Pads",
        "Developed - Low Intensity","Developed - Medium Intensity","Developed - High Intensity"
      ) ~ if (treat_nonveg_as_zero) 0 else NA_real_,
      TRUE ~ NA_real_
    ),
    # Open-ended bins → explicit numeric reps
    height_m = case_when(
      CLASSNAMES == "Tree Height >= 99 meters"   ~ 99,
      CLASSNAMES == "Shrub Height >= 3.0 meters" ~ 3,
      CLASSNAMES == "Herb Height >= 1 meter"     ~ 1,
      TRUE ~ height_m
    )
  )

# Row-aligned vectors (index = category ID 1..nrow(rat))
height_by_id   <- rat$height_m
lifeform_by_id <- rat$lifeform
class_by_id    <- rat$CLASSNAMES

# ---------------- STATIONS ----------------
# Read sheet (respect NULL = first sheet)
st_raw <- if (is.null(station_sheet)) {
  readxl::read_excel(station_xlsx)
} else {
  readxl::read_excel(station_xlsx, sheet = station_sheet)
}

# Normalize column names safely → station_id, lat, lon
nm <- names(st_raw)

# station_id
if (!"station_id" %in% nm) {
  src <- intersect(col_station_id_candidates, nm)
  if (length(src) == 0) stop("No station ID column found. Candidates: ", paste(col_station_id_candidates, collapse=", "))
  st_raw <- st_raw %>% rename(station_id = !!sym(src[1]))
}

# lat
if (!"lat" %in% names(st_raw)) {
  src <- intersect(col_lat_candidates, names(st_raw))
  if (length(src) == 0) stop("No latitude column found. Candidates: ", paste(col_lat_candidates, collapse=", "))
  st_raw <- st_raw %>% rename(lat = !!sym(src[1]))
}

# lon
if (!"lon" %in% names(st_raw)) {
  src <- intersect(col_lon_candidates, names(st_raw))
  if (length(src) == 0) stop("No longitude column found. Candidates: ", paste(col_lon_candidates, collapse=", "))
  st_raw <- st_raw %>% rename(lon = !!sym(src[1]))
}

# Clean & de-duplicate
st_raw <- st_raw %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>%
  filter(is.finite(lat), is.finite(lon)) %>%
  distinct(station_id, .keep_all = TRUE)

# Build vectors and project to EVH CRS
sv_wgs84 <- terra::vect(st_raw, geom = c("lon","lat"), crs = "EPSG:4326")
sv_5070  <- terra::project(sv_wgs84, "EPSG:5070")

# Keep only stations that fall on/inside EVH (drop AK/HI/AS, etc.)
inside_idx <- !is.na(terra::extract(evh_cat, sv_5070)[,2])
dropped    <- st_raw[!inside_idx, , drop = FALSE]
st_raw     <- st_raw[ inside_idx, , drop = FALSE]
sv_5070    <- sv_5070[inside_idx]

message("Stations total: ", nrow(st_raw) + nrow(dropped),
        " | inside EVH: ", nrow(st_raw),
        " | dropped (outside): ", nrow(dropped))

if (nrow(dropped) > 0) {
  readr::write_csv(dropped, out_csv_dropped)
  message("Wrote dropped stations: ", out_csv_dropped)
}

# ---------------- HELPERS ----------------
acres_to_radius_m <- function(ac) sqrt(43560 * ac / pi) * 0.3048  # ft^2 -> m
ft_to_m <- function(ft) ft * 0.3048
m_to_ft <- function(m)  m / 0.3048

pick_stat <- function(x, stat = "p95") {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  switch(stat,
         p95    = quantile(x, 0.95, na.rm = TRUE, names = FALSE),
         max    = max(x, na.rm = TRUE),
         mean   = mean(x, na.rm = TRUE),
         median = median(x, na.rm = TRUE),
         quantile(x, 0.95, na.rm = TRUE, names = FALSE))
}

# Analyze one station at one radius (meters)
analyze_one <- function(pt, r_m, use_ring_for_H = FALSE, ring_width_m = 100) {
  buf  <- terra::buffer(pt, width = r_m)
  ring <- if (use_ring_for_H) terra::buffer(pt, width = r_m + ring_width_m) - buf else NULL

  # Extract categorical values in buffer (df with ID, value)
  df_buf <- terra::extract(evh_cat, buf, df = TRUE)
  if (nrow(df_buf) == 0) {
    return(tibble::tibble(
      radius_m = r_m, clearing_acres = (pi * (m_to_ft(r_m))^2)/43560,
      n_cells = 0, pct_open = NA_real_, is_clear = NA,
      mean_h_m = NA_real_, median_h_m = NA_real_, p95_h_m = NA_real_, max_h_m = NA_real_,
      pct_veg = NA_real_, pct_nonveg = NA_real_, dominant_class = NA_character_,
      H_used_m = NA_real_, H_used_ft = NA_real_, X_factor = NA_real_
    ))
  }

  ids_buf <- as.integer(df_buf[[2]])  # category IDs (row order)
  # Guard invalid IDs (shouldn’t occur often)
  bad <- is.na(ids_buf) | ids_buf < 1 | ids_buf > length(height_by_id)
  ids_buf[bad] <- NA_integer_

  # Map to attributes
  h_buf   <- height_by_id[ids_buf]
  lf_buf  <- lifeform_by_id[ids_buf]
  cls_buf <- class_by_id[ids_buf]

  n_all   <- sum(is.finite(h_buf))
  n_open  <- sum(is.finite(h_buf) & h_buf <= clear_height_m)
  pct_open <- if (n_all > 0) n_open / n_all else NA_real_

  pct_veg    <- if (n_all > 0) sum(lf_buf %in% c("Tree","Shrub","Herb"), na.rm = TRUE) / n_all else NA_real_
  pct_nonveg <- if (n_all > 0) sum(lf_buf == "NonVeg", na.rm = TRUE) / n_all else NA_real_

  mean_h <- if (n_all > 0) mean(h_buf, na.rm = TRUE) else NA_real_
  med_h  <- if (n_all > 0) median(h_buf, na.rm = TRUE) else NA_real_
  p95_h  <- if (n_all > 0) quantile(h_buf, 0.95, na.rm = TRUE, names = FALSE) else NA_real_
  max_h  <- if (n_all > 0) max(h_buf, na.rm = TRUE) else NA_real_

  dom_class <- if (length(cls_buf) && any(!is.na(cls_buf))) {
    tb <- sort(table(cls_buf), decreasing = TRUE)
    names(tb)[1]
  } else NA_character_

  # Choose obstruction H for X-factor
  H_m <- if (!use_ring_for_H) {
    pick_stat(h_buf, H_stat)
  } else {
    df_ring <- terra::extract(evh_cat, ring, df = TRUE)
    ids_ring <- as.integer(df_ring[[2]])
    badr <- is.na(ids_ring) | ids_ring < 1 | ids_ring > length(height_by_id)
    ids_ring[badr] <- NA_integer_
    h_ring <- height_by_id[ids_ring]
    pick_stat(h_ring, H_stat)
  }
  H_ft <- if (is.finite(H_m)) m_to_ft(H_m) else NA_real_

  # Clearing size for this radius (in acres), then X = sqrt(43560*CS)/(2*H_ft)
  CS_ac <- (pi * (m_to_ft(r_m))^2) / 43560
  X     <- if (is.finite(H_ft) && H_ft > 0) sqrt(43560 * CS_ac) / (2 * H_ft) else NA_real_

  tibble::tibble(
    radius_m = r_m,
    clearing_acres = CS_ac,
    n_cells = n_all,
    pct_open = pct_open,
    is_clear = if (is.finite(pct_open)) pct_open >= clear_pct_required else NA,
    mean_h_m = mean_h,
    median_h_m = med_h,
    p95_h_m = p95_h,
    max_h_m = max_h,
    pct_veg = pct_veg,
    pct_nonveg = pct_nonveg,
    dominant_class = dom_class,
    H_used_m = H_m,
    H_used_ft = H_ft,
    X_factor = X
  )
}

# ---------------- RUN ANALYSIS ----------------
radii_m <- acres_to_radius_m(clearing_acres)

results <- purrr::map_dfr(seq_len(nrow(sv_5070)), function(i) {
  pt  <- sv_5070[i]
  row <- st_raw[i, ]

  per_r <- purrr::map_dfr(radii_m, ~ analyze_one(pt, .x, use_ring_for_H, ring_width_m))
  per_r %>%
    mutate(
      station_id   = row$station_id,
      station_name = row$station_name,
      lat = row$lat,
      lon = row$lon,
      x_5070 = terra::geom(pt)[1,1],
      y_5070 = terra::geom(pt)[1,2]
    ) %>%
    select(station_id, station_name, lat, lon, x_5070, y_5070, everything())
  
})

# ---------------- SAVE & SUMMARIES ----------------
readr::write_csv(results, out_csv_results)
message("Wrote: ", out_csv_results)

summary_table <- results %>%
  group_by(clearing_acres) %>%
  summarise(
    stations_evaluated = n_distinct(station_id),
    stations_clear     = sum(is_clear, na.rm = TRUE),
    pct_clear          = stations_clear / stations_evaluated,
    .groups = "drop"
  )
print(summary_table)

best_exposure <- results %>%
  group_by(station_id) %>%
  slice_max(X_factor, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(station_id, clearing_acres, H_used_m, X_factor, dominant_class)

print(best_exposure)
