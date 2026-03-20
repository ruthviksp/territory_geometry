rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(tidyverse)

## SETUP / HOUSEKEEPING
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
  shp_file = c("Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp", "TalChhapar_Area.shp")
)

## Output folder
out_dir <- file.path(script_dir, "processed_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Build master table of all files across all leks
files_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  
  cfg <- lek_configs[i, ]
  
  data_dirs <- list.dirs(file.path(root_dir, cfg$location), recursive = FALSE, full.names = TRUE)
  data_dirs <- data_dirs[grepl("_COORDINATES$", basename(data_dirs))]
  
  map_dfr(data_dirs, function(d) {
    
    data_label <- sub("_COORDINATES$", "", basename(d))
    csv_path <- list.files(d, pattern = paste0("_", cfg$suffix, "\\.csv$"), full.names = TRUE)
    
    if (length(csv_path) == 0) return(NULL)
    
    tibble(lek_id = cfg$lek_id,
           location = cfg$location,
           suffix = cfg$suffix,
           shp_file = cfg$shp_file,
           data_label = data_label,
           date = parse_label_to_date(data_label, month_lookup), 
           csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

## HELPER FUNCTIONS
## Compute KDE intensity surface and extract the global mode
compute_intensity_features <- function(lek_polygon, lek_points_sf) {
  
  # Convert points and polygons to spatstat and create a point pattern object
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(lek_points_sf)
  X <- ppp(xy[,1], xy[,2], window = W)
  
  # Estimate KDE bandwidth
  sigma <- bw.ppl(X)
  lambda_hat <- density.ppp(X, sigma = sigma, edge = TRUE, at = "pixels")
  
  # Extract location of maximum intensity (mode)
  v <- lambda_hat$v
  max_v <- max(v, na.rm = TRUE)
  idx_all <- which(v == max_v, arr.ind = TRUE)
  
  idx <- idx_all[1, , drop = FALSE]
  mode <- tibble(mx = lambda_hat$xcol[idx[2]], my = lambda_hat$yrow[idx[1]])
  
  return(mode)
}

## MAIN 
## Compute the shift in KDE mode between consecutive dates
stability_tbl <- map_dfr(unique(files_tbl$lek_id), function(lek) {
  
  # Subset data specific to lek
  sub_tbl <- files_tbl %>% filter(lek_id == lek) %>% arrange(date)
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    
    # Isolate rows for current and previous time steps
    row_prev <- sub_tbl[i-1,]
    row_now  <- sub_tbl[i,]
    
    # Read lek polygon
    lek_polygon <- st_read(file.path(root_dir, row_now$location, row_now$shp_file), quiet = TRUE) |> st_transform(32643) |> st_zm(drop = TRUE)
    
    # Extract points for current and previous time steps 
    pts_prev <- read.csv(row_prev$csv_path) |> st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    pts_now <- read.csv(row_now$csv_path) |> st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    
    # Calculate KDE at t-1 and t
    feat_prev <- compute_intensity_features(lek_polygon, pts_prev)
    feat_now  <- compute_intensity_features(lek_polygon, pts_now)
    
    # Compute Euclidean distance between the two locations
    mode_shift <- sqrt((feat_now$mx - feat_prev$mx)^2 + (feat_now$my - feat_prev$my)^2)
    
    tibble(lek_id = lek, date_prev = row_prev$date, date_now  = row_now$date, mode_shift = mode_shift)
  })
})

## Save output
out_file <- file.path(out_dir, "stability_ALL.csv")
write.csv(stability_tbl, out_file, row.names = FALSE)

message("Saved stability metrics to:", out_file)
