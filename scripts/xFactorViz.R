# ============================================================
# RAWS × EVH visualization (categorical workflow, no reclass)
# - Select stations by station_name (includes station_id)
# - Maps: EVH height (via factor-index mapping) + open-cells view
# - Histograms: EVH height distribution inside each buffer
# - Writes a tidy CSV of per-buffer metrics
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readxl)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(cowplot)
  library(scales)
  library(rlang)
})

# ---------------- CONFIG ----------------
evh_tif        <- "data/LF2023_EVH_240_CONUS/Tif/LC23_EVH_240.tif"  # categorical EVH with RAT
station_xlsx   <- "data/FEMS_3.0_RAWS_Master_Station_List_and_Metadata.xlsx"
station_sheet  <- NULL

# Column name candidates
col_station_id_candidates   <- c("station_id", "wrcc_id", "nesdis_id", "Station_ID", "StationID")
col_station_name_candidates <- c("station_name", "Station Name", "name", "Name")
col_lat_candidates          <- c("lat", "latitude", "LAT", "Latitude")
col_lon_candidates          <- c("lon", "longitude", "LON", "Longitude")

# Selection by station_name (exact match). If NULL, auto-pick first 2 inside EVH.
#station_names <- NULL  # e.g., c("JETTE", "ANOTHER STATION")
station_names <- c("RINCON","JETTE")

# Clearing sizes (acres) and map frame
clearing_acres <- c(0.1, 2, 10)
margin_m <- 100  # extra padding around largest buffer for map frame (m)

# “Clear” definition & X-factor settings (align w/ analysis script)
clear_height_m       <- 2 # vegetation height threshold for “open” (m)
clear_pct_required   <- 0.95
treat_nonveg_as_zero <- TRUE
H_stat               <- "p95"     # "p95", "max", "mean", "median"
use_ring_for_H       <- FALSE # calculations within buffer (FALSE) or outside ring (TRUE)
ring_width_m         <- 100 # if use_ring_for_H = TRUE

# Outputs
out_dir         <- "figs_raws_evh"
out_csv_metrics <- file.path(out_dir, "EVH_viz_metrics.csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- HELPERS ----------------
acres_to_radius_m <- function(ac) sqrt(43560 * ac / pi) * 0.3048
m_to_ft <- function(m) m / 0.3048

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

format_metrics_label <- function(df_row, H_stat = "p95") {
  sprintf("acres: %s\nopen: %s\nH(%s): %s m\nX: %s",
          number(df_row$clearing_acres, accuracy = 0.1),
          percent(df_row$pct_open, accuracy = 1),
          H_stat,
          number(df_row$H_used_m, accuracy = 0.1),
          number(df_row$X_factor, accuracy = 0.01))
}

# Analyze one point+radius using categorical EVH, mapping via row-aligned vectors
analyze_one <- function(evh_cat, height_by_id, lifeform_by_id, class_by_id,
                        pt, r_m, use_ring_for_H = FALSE, ring_width_m = 100,
                        clear_height_m = 0.5, clear_pct_required = 0.95,
                        H_stat = "p95") {
  buf  <- terra::buffer(pt, width = r_m)
  ring <- if (use_ring_for_H) terra::buffer(pt, width = r_m + ring_width_m) - buf else NULL
  
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
  
  idx <- as.integer(df_buf[[2]])                           # factor level indices
  bad <- is.na(idx) | idx < 1 | idx > length(height_by_id)
  idx[bad] <- NA_integer_
  
  h_buf   <- height_by_id[idx]
  lf_buf  <- lifeform_by_id[idx]
  cls_buf <- class_by_id[idx]
  
  n_all    <- sum(is.finite(h_buf))
  n_open   <- sum(is.finite(h_buf) & h_buf <= clear_height_m)
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
  
  H_m <- if (!use_ring_for_H) {
    pick_stat(h_buf, H_stat)
  } else {
    df_ring <- terra::extract(evh_cat, ring, df = TRUE)
    idx_ring <- as.integer(df_ring[[2]])
    badr <- is.na(idx_ring) | idx_ring < 1 | idx_ring > length(height_by_id)
    idx_ring[badr] <- NA_integer_
    h_ring <- height_by_id[idx_ring]
    pick_stat(h_ring, H_stat)
  }
  H_ft <- if (is.finite(H_m)) m_to_ft(H_m) else NA_real_
  
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

# Build plotting data by mapping factor-level indices → height/class (no joins)
evh_df_for_plot <- function(evh_crop, height_by_id, class_by_id) {
  df <- as.data.frame(evh_crop, xy = TRUE, cells = FALSE, na.rm = FALSE)
  val_col <- names(df)[3]                 # categorical values column (factor)
  idx <- as.integer(df[[val_col]])        # factor level indices: 1..nlevels
  bad <- is.na(idx) | idx < 1 | idx > length(height_by_id)
  
  df$height_m   <- NA_real_
  df$CLASSNAMES <- NA_character_
  
  df$height_m[!bad]   <- height_by_id[idx[!bad]]
  df$CLASSNAMES[!bad] <- class_by_id[idx[!bad]]
  df
}

# ---------------- LOAD EVH ----------------
stopifnot(file.exists(evh_tif))
evh_cat <- terra::rast(evh_tif)
stopifnot(terra::is.factor(evh_cat))

# Ensure meters CRS (EPSG:5070) for correct buffering
if (!grepl("5070", terra::crs(evh_cat))) {
  message("Reprojecting EVH to EPSG:5070 (near-neighbor)...")
  evh_cat <- terra::project(evh_cat, "EPSG:5070", method = "near")
} else {
  message("EVH already in EPSG:5070.")
}

# RAT → row-aligned lookup vectors (index = category ID)
rat <- terra::levels(evh_cat)[[1]]
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
        "NASS-Row Crop","NASS-Row Crop-Close Grown Crop","NASS-Close Grown Crop","NASS-Wheat",
        "Quarries-Strip Mines-Gravel Pits-Well and Wind Pads",
        "Developed - Low Intensity","Developed - Medium Intensity","Developed - High Intensity"
      ) ~ if (treat_nonveg_as_zero) 0 else NA_real_,
      TRUE ~ NA_real_
    ),
    height_m = case_when(
      CLASSNAMES == "Tree Height >= 99 meters"   ~ 99,
      CLASSNAMES == "Shrub Height >= 3.0 meters" ~ 3,
      CLASSNAMES == "Herb Height >= 1 meter"     ~ 1,
      TRUE ~ height_m
    )
  )

height_by_id   <- rat$height_m
lifeform_by_id <- rat$lifeform
class_by_id    <- rat$CLASSNAMES

# ---------------- LOAD STATIONS ----------------
st_raw <- if (is.null(station_sheet)) readxl::read_excel(station_xlsx) else readxl::read_excel(station_xlsx, sheet = station_sheet)
nm <- names(st_raw)

# Normalize fields (keep both name and id)
if (!"station_id" %in% nm) {
  src <- intersect(col_station_id_candidates, nm)
  if (length(src) == 0) stop("No station_id column found. Candidates: ", paste(col_station_id_candidates, collapse=", "))
  st_raw <- st_raw %>% rename(station_id = !!sym(src[1]))
}
if (!"station_name" %in% names(st_raw)) {
  src <- intersect(col_station_name_candidates, names(st_raw))
  if (length(src) == 0) stop("No station_name column found. Candidates: ", paste(col_station_name_candidates, collapse=", "))
  st_raw <- st_raw %>% rename(station_name = !!sym(src[1]))
}
if (!"lat" %in% names(st_raw)) {
  src <- intersect(col_lat_candidates, names(st_raw)); if (!length(src)) stop("No latitude column found.")
  st_raw <- st_raw %>% rename(lat = !!sym(src[1]))
}
if (!"lon" %in% names(st_raw)) {
  src <- intersect(col_lon_candidates, names(st_raw)); if (!length(src)) stop("No longitude column found.")
  st_raw <- st_raw %>% rename(lon = !!sym(src[1]))
}

st_raw <- st_raw %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>%
  filter(is.finite(lat), is.finite(lon)) %>%
  distinct(station_id, .keep_all = TRUE)

sv_wgs84 <- terra::vect(st_raw, geom = c("lon","lat"), crs = "EPSG:4326")
sv_5070  <- terra::project(sv_wgs84, "EPSG:5070")

# Keep only stations that fall on/inside EVH
inside_idx <- !is.na(terra::extract(evh_cat, sv_5070)[,2])
st_keep <- st_raw[inside_idx, , drop = FALSE]
sv_5070 <- sv_5070[inside_idx]
stopifnot(nrow(st_keep) >= 2)

# ---- Selection by station_name ----
if (!is.null(station_names)) {
  sel <- st_keep$station_name %in% station_names
  if (sum(sel) < 2) stop("Fewer than 2 selected stations found inside EVH.")
  st_sel <- st_keep[sel, ]
  sv_sel <- sv_5070[sel]
} else {
  # auto-pick first 2 (change to sample() for random)
  st_sel <- st_keep %>% slice(1:2)
  rows <- match(st_sel$station_id, st_keep$station_id)
  sv_sel <- sv_5070[rows]
}

# ---------------- MAPS + METRICS + HISTOGRAMS (all in EPSG:5070) ----------------
radii_m <- acres_to_radius_m(clearing_acres)
all_metrics <- list()

for (i in seq_len(nrow(st_sel))) {
  st_row <- st_sel[i, ]
  pt     <- sv_sel[i]  # SpatVector POINT in EPSG:5070
  
  # --- projected coords (meters) ---
  coords <- terra::crds(pt, df = TRUE)
  x_5070_val <- coords$x
  y_5070_val <- coords$y
  
  # Buffers (5070)
  bufs <- lapply(radii_m, function(rm) terra::buffer(pt, width = rm))
  bufs_sv <- do.call(rbind, bufs)
  bufs_sf <- sf::st_as_sf(bufs_sv); sf::st_crs(bufs_sf) <- 5070
  bufs_sf$clearing_acres <- clearing_acres
  
  # Frame & crop (5070)
  big_r    <- max(radii_m)
  frame    <- terra::buffer(pt, width = big_r + margin_m)
  evh_crop <- terra::crop(evh_cat, frame)
  
  # Plotting data via factor-index mapping (no joins; 5070)
  evh_df <- evh_df_for_plot(evh_crop, height_by_id, class_by_id) %>%
    dplyr::mutate(is_open_cell = dplyr::case_when(
      is.finite(height_m) & height_m <= clear_height_m ~ "TRUE",
      is.finite(height_m) & height_m >  clear_height_m ~ "FALSE",
      TRUE ~ "NA"
    ))
  
  # Station point (sf, 5070)
  pt_sf <- sf::st_as_sf(coords, coords = c("x","y"), crs = 5070)
  pt_sf$station_id   <- st_row$station_id
  pt_sf$station_name <- st_row$station_name
  
  # Metrics per buffer (5070)
  metrics <- purrr::map_dfr(seq_along(radii_m), function(j) {
    analyze_one(
      evh_cat, height_by_id, lifeform_by_id, class_by_id,
      pt, radii_m[j], use_ring_for_H, ring_width_m,
      clear_height_m, clear_pct_required, H_stat
    )
  }) %>%
    dplyr::mutate(
      station_id   = st_row$station_id,
      station_name = st_row$station_name,
      lat = st_row$lat, lon = st_row$lon,
      x_5070 = x_5070_val, y_5070 = y_5070_val
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(label = format_metrics_label(cur_data(), H_stat)) %>%
    dplyr::ungroup()
  
  all_metrics[[i]] <- metrics
  
  bufs_sf <- dplyr::left_join(
    bufs_sf,
    metrics %>% dplyr::select(clearing_acres, pct_open, H_used_m, X_factor, label),
    by = "clearing_acres"
  )
  
  # --- label placement (auto-distributed, 5070 meters) ---
  n_r <- length(radii_m)
  start_deg  <- 20
  angles_deg <- seq(start_deg, start_deg + 360 - 360/n_r, length.out = n_r)
  angles     <- angles_deg * pi/180
  label_pad  <- 50  # meters beyond ring
  
  lab_xy <- tibble::tibble(
    clearing_acres = clearing_acres,
    x  = x_5070_val + (radii_m + label_pad) * cos(angles),
    y  = y_5070_val + (radii_m + label_pad) * sin(angles),
    xr = x_5070_val +  radii_m            * cos(angles),
    yr = y_5070_val +  radii_m            * sin(angles),
    label = metrics$label
  )
  
  # ======================= CLEAN MAPS + LEGENDS + LAYOUT =======================
  
  # ---- 1) Build categorical EVH df and palette (only classes present) ----
  evh_cat_df <- {
    df <- as.data.frame(evh_crop, xy = TRUE, na.rm = FALSE)
    val_col <- names(df)[3]
    idx <- as.integer(df[[val_col]])
    bad <- is.na(idx) | idx < 1 | idx > length(class_by_id)
    df$CLASSNAMES <- NA_character_
    df$lifeform   <- NA_character_
    df$height_m   <- NA_real_
    df$CLASSNAMES[!bad] <- class_by_id[idx[!bad]]
    df$lifeform[!bad]   <- lifeform_by_id[idx[!bad]]
    df$height_m[!bad]   <- height_by_id[idx[!bad]]
    df
  }
  
  # ===================== EVH categorical palette (drop-in) =====================
  
  # Build a distinct class table for palette construction
  # (one row per EVH class in the crop, with its lifeform + a representative height)
  pal_df <- evh_cat_df %>%
    dplyr::filter(!is.na(CLASSNAMES)) %>%
    dplyr::distinct(CLASSNAMES, lifeform, height_m)
  
  
  # Classes present in this crop (for final filtering)
  present_classes <- sort(unique(stats::na.omit(evh_cat_df$CLASSNAMES)))
  
  # Helper: smooth color ramps
  mk_ramp <- function(n, cols) if (n > 0) grDevices::colorRampPalette(cols)(n) else character(0)
  
  # Fixed colors for clearly non-veg classes
  nonveg_fixed <- c(
    "Open Water" = "#2C7FB8","Snow/Ice" = "#D9F0F3",
    "Developed-Roads" = "#6E6E6E","Developed - Low Intensity"="#969696",
    "Developed - Medium Intensity"="#7A7A7A","Developed - High Intensity"="#4D4D4D",
    "Barren"="#D9C29E","Sparse Vegetation Canopy"="#E6E6E6",
    "Quarries-Strip Mines-Gravel Pits-Well and Wind Pads"="#BDBDBD",
    "Cultivated Crops"="#C2E699","NASS-Row Crop"="#B5E3A1",
    "NASS-Row Crop-Close Grown Crop"="#A1D99B","NASS-Close Grown Crop"="#8BC57F",
    "NASS-Wheat"="#DCCB7A"
  )
  
  # Split the distinct class table by lifeform
  tree_df   <- pal_df %>% dplyr::filter(lifeform == "Tree")  %>% dplyr::arrange(height_m)
  shrub_df  <- pal_df %>% dplyr::filter(lifeform == "Shrub") %>% dplyr::arrange(height_m)
  herb_df   <- pal_df %>% dplyr::filter(lifeform == "Herb")  %>% dplyr::arrange(height_m)
  nonveg_df <- pal_df %>% dplyr::filter(is.na(lifeform) | lifeform == "NonVeg")
  
  # Color ramps by lifeform (taller = darker)
  tree_cols  <- mk_ramp(nrow(tree_df),  c("#D9F0D3","#41AB5D","#005A32"))
  shrub_cols <- mk_ramp(nrow(shrub_df), c("#FEE8C8","#FDD49E","#D7301F"))
  herb_cols  <- mk_ramp(nrow(herb_df),  c("#FFF7BC","#B8E186","#4D9221"))
  
  # Safe lookup for non-veg colors (fallback to light gray if unknown)
  normalize_key <- function(s) {
    s <- as.character(s)
    s <- trimws(s)
    s <- gsub("[\u2013\u2014]", "-", s)   # normalize en/em dashes to hyphen
    s
  }
  
  nonveg_cols <- if (nrow(nonveg_df)) {
    vapply(nonveg_df$CLASSNAMES, function(k) {
      if (is.na(k)) return("#E0E0E0")
      key <- normalize_key(k)
      val <- unname(nonveg_fixed[key])    # name indexing; returns NA if missing
      if (length(val) == 0L || is.na(val)) "#E0E0E0" else val
    }, FUN.VALUE = character(1L))
  } else {
    character(0)
  }
  
  # Optional: log any unmatched non-veg classes (once per station)
  unknown_nonveg <- setdiff(normalize_key(nonveg_df$CLASSNAMES), names(nonveg_fixed))
  unknown_nonveg <- unknown_nonveg[!is.na(unknown_nonveg)]
  if (length(unknown_nonveg)) {
    message("Non-veg classes not in palette (default grey): ",
            paste(unique(unknown_nonveg), collapse = "; "))
  }
  
  # Combine into a single named palette
  pal <- c(
    setNames(tree_cols,   tree_df$CLASSNAMES),
    setNames(shrub_cols,  shrub_df$CLASSNAMES),
    setNames(herb_cols,   herb_df$CLASSNAMES),
    setNames(nonveg_cols, nonveg_df$CLASSNAMES)
  )
  
  # Keep only classes that actually appear in the cropped raster
  pal <- pal[intersect(present_classes, names(pal))]
  # ============================================================================ 
  
  
  # ---- 2) Common theme & helpers ----
  leg_slim <- theme(
    legend.title = element_text(size = 9, hjust = 0),
    legend.text  = element_text(size = 8),
    legend.key.height = unit(0.32, "cm"),
    legend.key.width  = unit(0.32, "cm"),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9),
    plot.margin = margin(4, 4, 4, 4)
  )
  
  get_legend <- function(p) {
    g <- ggplotGrob(p + theme(legend.position = "right"))
    g$grobs[[which(sapply(g$grobs, \(x) x$name) == "guide-box")]]
  }
  
  # ---- 3) Three panels (no legends inside) ----
  p_height_core <- ggplot() +
    geom_raster(data = evh_df, aes(x = x, y = y, fill = height_m)) +
    scale_fill_viridis_c(name = "EVH height (m)", na.value = "grey90") +
    geom_sf(data = bufs_sf, aes(color = factor(round(clearing_acres, 2))),
            fill = NA, linewidth = 0.8, show.legend = FALSE) +
    geom_sf(data = pt_sf, shape = 21, size = 2.3, stroke = 0.5, fill = "white") +
    geom_segment(data = lab_xy, aes(x = xr, y = yr, xend = x, yend = y),
                 linewidth = 0.3, color = "grey30", inherit.aes = FALSE) +
    geom_text(data = lab_xy, aes(x = x, y = y, label = label), size = 3, lineheight = 0.95) +
    coord_sf(crs = 5070, expand = FALSE, clip = "off") +
    labs(
      title = paste0(st_row$station_name, " (", st_row$station_id, ")"),
      subtitle = paste0("EVH height (m); clear ≤ ", clear_height_m,
                        " m; H = ", H_stat, if (use_ring_for_H) " (ring)" else " (inside)")
    ) +
    theme_minimal(base_size = 12) + leg_slim + theme(legend.position = "none")
  
  p_open_core <- ggplot() +
    geom_raster(data = evh_df, aes(x = x, y = y, fill = is_open_cell)) +
    scale_fill_manual(values = c("FALSE"="grey75","TRUE"="white","NA"="grey90"),
                      breaks = c("TRUE","FALSE"),
                      labels = c("≤ clear height", "> clear height"),
                      name = "Cell open?") +
    geom_sf(data = bufs_sf, color = "black", fill = NA, linewidth = 0.6) +
    geom_sf(data = pt_sf, shape = 21, size = 2.3, stroke = 0.5, fill = "white") +
    geom_segment(data = lab_xy, aes(x = xr, y = yr, xend = x, yend = y),
                 linewidth = 0.3, color = "grey30", inherit.aes = FALSE) +
    geom_text(data = lab_xy, aes(x = x, y = y, label = label), size = 3, lineheight = 0.95) +
    coord_sf(crs = 5070, expand = FALSE, clip = "off") +
    labs(title = "Open cells", subtitle = paste0("Threshold ≤ ", clear_height_m, " m")) +
    theme_minimal(base_size = 12) + leg_slim + theme(legend.position = "none")
  
  p_evh_cat_core <- ggplot() +
    geom_raster(data = evh_cat_df, aes(x = x, y = y, fill = CLASSNAMES)) +
    scale_fill_manual(values = pal, breaks = names(pal), drop = TRUE, name = "EVH class") +
    geom_sf(data = bufs_sf, fill = NA, color = "black", linewidth = 0.6) +
    geom_sf(data = pt_sf, shape = 21, size = 2.3, stroke = 0.5, fill = "white") +
    coord_sf(crs = 5070, expand = FALSE, clip = "off") +
    labs(title = "EVH classes", subtitle = "Within map frame") +
    theme_minimal(base_size = 12) + leg_slim + theme(legend.position = "none")
  
  # ---- 4) Build tidy legends ----
  # ring (clearing sizes) + height legend (stacked)
  p_ring_leg <- ggplot() +
    geom_sf(data = bufs_sf, aes(color = factor(round(clearing_acres, 2))),
            fill = NA, linewidth = 0.8) +
    scale_color_discrete(name = "Clearing (ac)") +
    theme_void() + theme(legend.position = "right") + leg_slim
  
  p_height_leg <- ggplot() +
    geom_raster(data = evh_df, aes(x = x, y = y, fill = height_m)) +
    scale_fill_viridis_c(name = "EVH height (m)", na.value = "grey90") +
    theme_void() + theme(legend.position = "right") + leg_slim
  
  p_open_leg <- ggplot() +
    geom_raster(data = evh_df, aes(x = x, y = y, fill = is_open_cell)) +
    scale_fill_manual(values = c("FALSE"="grey75","TRUE"="white","NA"="grey90"),
                      breaks = c("TRUE","FALSE"),
                      labels = c("≤ clear height", "> clear height"),
                      name = "Cell open?") +
    theme_void() + theme(legend.position = "right") + leg_slim
  
  wrap_lab <- function(x, width = 22) stringr::str_wrap(x, width)
  p_cat_leg <- ggplot() +
    geom_raster(data = evh_cat_df, aes(x = x, y = y, fill = CLASSNAMES)) +
    scale_fill_manual(values = pal, breaks = names(pal), labels = wrap_lab(names(pal), 22),
                      name = "EVH class") +
    theme_void() + theme(legend.position = "right") + leg_slim +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
  
  leg_ring   <- get_legend(p_ring_leg)
  leg_height <- get_legend(p_height_leg)
  leg_open   <- get_legend(p_open_leg)
  leg_cat    <- get_legend(p_cat_leg)
  
  # ---- 5) Pair each panel with its legend (keeps legends close to plots) ----
  height_leg_col <- cowplot::plot_grid(leg_ring, leg_height, ncol = 1, rel_heights = c(.55, .45))
  height_pair <- cowplot::plot_grid(p_height_core, height_leg_col, ncol = 2, rel_widths = c(1, 0.33), align = "h")
  open_pair <- cowplot::plot_grid(p_open_core, leg_open, ncol = 2, rel_widths = c(1, 0.30))
  cat_pair  <- cowplot::plot_grid(p_evh_cat_core, leg_cat, ncol = 2, rel_widths = c(1, 0.52))
  
  # ---- 6) Stack vertically (portrait page) and save ----
  g_map <- cowplot::plot_grid(height_pair, open_pair, cat_pair, ncol = 1,
                              rel_heights = c(1, 1, 1.05), align = "v")

  outf_map <- file.path(
    out_dir,
    paste0(
      "EVH_buffers_",
      st_row$station_id, "_",
      gsub("[^A-Za-z0-9]+", "_", st_row$station_name),
      ".png"
    )
  )
  
  ggsave(outf_map, g_map, width = 8.5, height = 11, dpi = 300, bg = "white")
  # ============================================================================
  
  # ---- Histograms (unchanged; 5070 analysis) ----
  hist_df <- purrr::map2_dfr(bufs, clearing_acres, function(b, ac) {
    vals <- terra::extract(evh_cat, b, df = TRUE)[[2]]
    idx  <- as.integer(vals)
    idx[!is.finite(idx)] <- NA_integer_
    h <- height_by_id[idx]
    tibble::tibble(clearing_acres = ac, height_m = h[is.finite(h)])
  })
  
  hist_summ <- hist_df %>%
    dplyr::group_by(clearing_acres) %>%
    dplyr::summarise(
      mean_h = mean(height_m, na.rm = TRUE),
      med_h  = median(height_m, na.rm = TRUE),
      p95_h  = stats::quantile(height_m, 0.95, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
  
  p_hist <- ggplot(hist_df, aes(x = height_m)) +
    geom_histogram(bins = 40) +
    facet_wrap(~ clearing_acres, scales = "free_y", nrow = 2,
               labeller = labeller(clearing_acres = function(x) paste0(x, " ac"))) +
    geom_vline(data = hist_summ, aes(xintercept = mean_h), linetype = "solid") +
    geom_vline(data = hist_summ, aes(xintercept = med_h),  linetype = "dashed") +
    geom_vline(data = hist_summ, aes(xintercept = p95_h),  linetype = "dotted") +
    labs(
      title = paste0("EVH height distributions — ", st_row$station_name, " (", st_row$station_id, ")"),
      subtitle = "Vertical lines: mean (solid), median (dashed), p95 (dotted)",
      x = "EVH height (m)",
      y = "Cell count"
    ) +
    theme_minimal(base_size = 12)
  
  outf_hist <- file.path(
    out_dir,
    paste0(
      "EVH_hist_",
      st_row$station_id, "_",
      gsub("[^A-Za-z0-9]+", "_", st_row$station_name),
      ".png"
    )
  )
  ggsave(outf_hist, p_hist, width = 12, height = 7.0, dpi = 200, bg = "white")
  message("Wrote hist: ", outf_hist)
}


# ---- Tidy CSV for all plotted stations ----
all_metrics_df <- dplyr::bind_rows(all_metrics) %>%
  dplyr::select(
    station_id, station_name, lat, lon, x_5070, y_5070,
    radius_m, clearing_acres, n_cells, pct_open, is_clear,
    mean_h_m, median_h_m, p95_h_m, max_h_m,
    pct_veg, pct_nonveg, dominant_class,
    H_used_m, H_used_ft, X_factor
  )

readr::write_csv(all_metrics_df, out_csv_metrics)
message("Wrote metrics: ", out_csv_metrics)
message("Done.")





