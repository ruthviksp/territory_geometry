rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(tidyverse)

## Get directory where this script lives
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## Lookup table to convert file names to dates
month_lookup <- c(JAN = 1, FEB = 2, MAR = 3, APR = 4,
                  MAY = 5, JUN = 6, JUL = 7, AUG = 8,
                  SEP = 9, OCT = 10, NOV = 11, DEC = 12)

parse_label_to_date <- function(data_label, month_lookup) {
  yy <- as.integer(substr(data_label, 1, 2))
  mm <- unname(month_lookup[substr(data_label, 3, 5)])
  as.Date(sprintf("20%02d-%02d-01", yy, mm))
}

## Root code and data directory
root_dir <- "/Users/vivekhsridhar/Library/Mobile Documents/com~apple~CloudDocs/Documents/Data/SatelliteImagery/GoogleEarth"
setwd(root_dir)

## Lek configuration table
lek_configs <- tibble(
  lek_id   = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
  location = c("Velavadar", "Velavadar", "TalChhapar"),
  suffix   = c("LEK1", "LEK2", "TC"),
  shp_file = c("Velavadar_Lek1_Area.shp",
               "Velavadar_Lek2_Area.shp",
               "TalChhapar_Area.shp")
)

## Output folder
out_dir <- file.path(script_dir, "processed_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Build master table of all files across all leks
files_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  
  cfg <- lek_configs[i, ]
  
  data_dirs <- list.dirs(file.path(root_dir, cfg$location),
                         recursive = FALSE, full.names = TRUE)
  data_dirs <- data_dirs[grepl("_COORDINATES$", basename(data_dirs))]
  
  map_dfr(data_dirs, function(d) {
    
    data_label <- sub("_COORDINATES$", "", basename(d))
    csv_path <- list.files(d,
                           pattern = paste0("_", cfg$suffix, "\\.csv$"),
                           full.names = TRUE)
    
    if (length(csv_path) == 0) return(NULL)
    
    tibble(
      lek_id = cfg$lek_id,
      location = cfg$location,
      suffix = cfg$suffix,
      shp_file = cfg$shp_file,
      data_label = data_label,
      date = parse_label_to_date(data_label, month_lookup),
      csv_path = csv_path[1]
    )
  })
}) %>% arrange(lek_id, date)

## Helper: compute intensity centroid and mode
compute_intensity_features <- function(lek_polygon, lek_points_sf) {
  
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(lek_points_sf)
  X <- ppp(xy[,1], xy[,2], window = W)
  
  sigma <- bw.ppl(X)
  lambda_hat <- density.ppp(X, sigma = sigma, edge = TRUE, at = "pixels")
  
  ## Global mode of intensity surface (robust to NA / flat surfaces)
  v <- lambda_hat$v
  
  max_v <- max(v, na.rm = TRUE)
  idx_all <- which(v == max_v, arr.ind = TRUE)
  
  idx <- idx_all[1, , drop = FALSE]
  mode <- tibble(
    mx = lambda_hat$xcol[idx[2]],
    my = lambda_hat$yrow[idx[1]]
  )
  
  mode
}

## Helper: cross-year nearest-neighbour distances
cross_year_nn <- function(pts_now, pts_prev, W) {
  
  X_now  <- ppp(pts_now[,1],  pts_now[,2],  window = W)
  X_prev <- ppp(pts_prev[,1], pts_prev[,2], window = W)
  
  nn1 <- nncross(X_now,  X_prev)
  nn2 <- nncross(X_prev, X_now)
  
  ## Extract distance component safely
  d1 <- if (is.list(nn1)) nn1$dist else nn1
  d2 <- if (is.list(nn2)) nn2$dist else nn2
  
  d <- c(d1, d2)
  d <- as.numeric(d)
  d <- d[is.finite(d)]
  
  tibble(
    nn_cross_median = median(d),
    nn_cross_mean   = mean(d),
    nn_cross_cv     = sd(d) / mean(d)
  )
}

## Main loop: compute stability metrics for consecutive years
stability_tbl <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  
  sub_tbl <- files_tbl %>% filter(lek_id == lk) %>% arrange(date)
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    
    row_prev <- sub_tbl[i - 1, ]
    row_now  <- sub_tbl[i, ]
    
    lek_polygon <- st_read(file.path(root_dir, row_now$location, row_now$shp_file), quiet = TRUE) |> st_transform(32643) |> st_zm(drop = TRUE)
    
    pts_prev <- read.csv(row_prev$csv_path) |> st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    pts_now <- read.csv(row_now$csv_path) |> st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    
    feat_prev <- compute_intensity_features(lek_polygon, pts_prev)
    feat_now  <- compute_intensity_features(lek_polygon, pts_now)
    
    mode_shift <- sqrt((feat_now$mx - feat_prev$mx)^2 + (feat_now$my - feat_prev$my)^2)
    
    W <- as.owin(st_geometry(lek_polygon))
    
    nn_cross <- cross_year_nn(st_coordinates(pts_now), st_coordinates(pts_prev), W)
    
    tibble(lek_id = lk, date_prev = row_prev$date, date_now  = row_now$date,
           mode_shift = mode_shift) %>% bind_cols(nn_cross)
  })
})

## Save output
out_file <- file.path(out_dir, "stability_ALL.csv")
write.csv(stability_tbl, out_file, row.names = FALSE)

message("Saved stability metrics to:", out_file)
