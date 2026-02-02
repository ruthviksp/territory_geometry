rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(spatstat.geom)
library(dplyr)

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

## Compute NND metrics
compute_nnd <- function(lek_polygon, lek_points) {
  W <- as.owin(st_geometry(lek_polygon))
  pts <- st_coordinates(lek_points)
  X <- ppp(pts[,1], pts[,2], window = W)
  
  nn <- nndist(X)
  
  tibble(nnd_mean = mean(nn), nnd_median = median(nn), nnd_sd = sd(nn), nnd_count = length(nn), nnd_cv = sd(nn)/mean(nn))
}

## Root code and data directory
root_dir <- "/Users/vivekhsridhar/Library/Mobile Documents/com~apple~CloudDocs/Documents/Data/SatelliteImagery/GoogleEarth"
setwd(root_dir)

## Lek configuration table
lek_configs <- tibble(lek_id   = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
                      location = c("Velavadar", "Velavadar", "TalChhapar"), suffix   = c("LEK1", "LEK2", "TC"),
                      shp_file = c("Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp", "TalChhapar_Area.shp"))

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
    
    tibble(lek_id = cfg$lek_id, location = cfg$location, suffix = cfg$suffix, shp_file = cfg$shp_file, 
           data_label = data_label, date = parse_label_to_date(data_label, month_lookup), 
           csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

## Compute NND metrics for all leks and all dates
nnd_results <- map_dfr(seq_len(nrow(files_tbl)), function(i) {
  row <- files_tbl[i, ]
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |> st_transform(32643) |> st_zm(drop = TRUE)
  lek_points <- read.csv(row$csv_path) |> st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
  
  ## Compute NND metrics, append columns for date and data label
  ## map_dfr then stacks rows each time we run this loop
  compute_nnd(lek_polygon, lek_points) |> mutate(lek_id = row$lek_id, data_label = row$data_label, date = row$date)
})

## Save results
out_file <- file.path(out_dir, "nnd_ALL.csv")
write.csv(nnd_results, out_file, row.names = FALSE)

message("Saved NND metrics to: ", out_file)