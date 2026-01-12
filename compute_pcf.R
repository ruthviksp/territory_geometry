rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(tidyverse)
library(zoo)

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

## Compute g_inhom(r)
compute_pcf <- function(lek_polygon, lek_points_sf, r_vals,
                        r_min = 1.2, correction = "translate") {
  
  # Define observation window
  W <- as.owin(st_geometry(lek_polygon))
  
  # Create point pattern object
  pts <- st_coordinates(lek_points_sf)
  X <- ppp(pts[, 1], pts[, 2], window = W)
  stopifnot(all(inside.owin(X$x, X$y, W)))
  
  # Nearest-neighbour distances (for summary)
  nn <- nndist(X)
  nn_median <- median(nn)
  
  # Intensity estimate (inhomogeneous)
  sigma <- bw.ppl(X)
  lambda_hat <- density.ppp(X, sigma = sigma, edge = TRUE, at = "pixels")
  
  # Inhomogeneous pair correlation on a common r-grid
  g <- pcfinhom(X, lambda = lambda_hat, r = r_vals, correction = correction)
  
  # Use translate correction values by default
  g_df <- tibble(r = g$r, g = g$trans) %>% filter(r > r_min)
  
  list(g_df = g_df, nn_median = nn_median, sigma = sigma)
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

## Comparability controls
scale_mode <- "global"
r_max_mult <- 4
n_r <- 240
r_min <- 5.0
correction <- "translate"

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
           data_label = data_label, date = parse_label_to_date(data_label, month_lookup), csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

if (nrow(files_tbl) == 0) {
  stop("No CSV files found across leks.")
}

## First pass: compute a GLOBAL reference median NND for comparability
nn_pass <- map_dfr(seq_len(nrow(files_tbl)), function(i) {
  
  row <- files_tbl[i, ]
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |> st_transform(32643) |> st_zm(drop = TRUE)
  
  df <- read.csv(row$csv_path)
  pts_sf <- st_as_sf(df, coords = c("pos_x", "pos_y"), crs = 32643)
  
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  X <- ppp(xy[, 1], xy[, 2], window = W)
  
  tibble(nn_median = median(nndist(X)))
})

files_tbl <- bind_cols(files_tbl, nn_pass)

ref_median_nn <- median(files_tbl$nn_median, na.rm = TRUE)

r_max  <- r_max_mult * ref_median_nn
r_vals <- seq(0, r_max, length.out = n_r)

message("GLOBAL reference median NND = ", round(ref_median_nn, 2), " m")
message("Using GLOBAL r_max = ", round(r_max, 2), " m")
message("Using ", n_r, " r-values")

## Second pass: compute PCFs for all leks on the SAME r-grid
curve_list <- list()
summary_list <- list()

for (i in seq_len(nrow(files_tbl))) {
  
  row <- files_tbl[i, ]
  message("Processing ", row$lek_id, " : ", row$data_label, " (", i, "/", nrow(files_tbl), ")")
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |> st_transform(32643) |> st_zm(drop = TRUE)
  
  df <- read.csv(row$csv_path)
  lek_points <- st_as_sf(df, coords = c("pos_x", "pos_y"), crs = 32643)
  
  res <- compute_pcf(lek_polygon = lek_polygon, lek_points_sf = lek_points,
                     r_vals = r_vals, r_min = r_min, correction = correction)
  
  curve_list[[i]] <- res$g_df %>% mutate(lek_id = row$lek_id, data_label = row$data_label,
                                         date = row$date, r_max_used = r_max, ref_median_nn = ref_median_nn)
  
  summary_list[[i]] <- tibble(lek_id = row$lek_id, data_label = row$data_label,
                              date = row$date, n_points = nrow(df), nn_median = res$nn_median, 
                              bw_sigma = as.numeric(res$sigma), ref_median_nn = ref_median_nn, r_max_used = r_max)
}

pcf_curves  <- bind_rows(curve_list)
pcf_summary <- bind_rows(summary_list)

## Write outputs
curves_out_file <- file.path(out_dir, "pcf_curve_ALL.csv")
summary_out_file <- file.path(out_dir, "pcf_summary_ALL.csv")
write.csv(pcf_curves, curves_out_file, row.names = FALSE)
write.csv(pcf_summary, summary_out_file, row.names = FALSE)

message("Saved curves to:", curves_out_file)
message("Saved summary to:", summary_out_file)

## Peak detection parameters
lower_nnd_mult <- 0.8
smooth_k <- 5
min_prominence <- 0.02
min_sep_mult <- 0.5

## Peak detection function
detect_peaks <- function(r, g, med_nnd) {
  
  keep <- r >= lower_nnd_mult * med_nnd
  r <- r[keep]
  g <- g[keep]
  
  if (length(r) < 10) return(NULL)
  
  # Smooth only for locating peaks (NOT for measuring shape)
  g_s <- zoo::rollmean(g, k = smooth_k, fill = NA, align = "center")
  
  is_peak <- g_s > dplyr::lag(g_s) & g_s > dplyr::lead(g_s)
  peak_idx <- which(is_peak)
  if (length(peak_idx) == 0) return(NULL)
  
  peaks <- tibble(
    idx = peak_idx,
    r_peak = r[peak_idx],
    g_peak = g[peak_idx]
  ) %>%
    arrange(r_peak) %>%
    filter(c(TRUE, diff(r_peak) >= min_sep_mult * med_nnd))
  
  # Helper: choose a local window (in r units) around each peak
  # Use something proportional to spacing so it scales across patterns
  win_half_width <- 0.75 * med_nnd
  
  dr <- median(diff(r), na.rm = TRUE)
  if (!is.finite(dr) || dr <= 0) dr <- diff(range(r)) / (length(r) - 1)
  
  # Compute prominence, curvature, and width for each peak
  out <- purrr::pmap_dfr(peaks, function(idx, r_peak, g_peak, peak_height_above_mean) {
    
    # Window bounds in index space
    left_limit_r  <- r_peak - win_half_width
    right_limit_r <- r_peak + win_half_width
    
    left_idx  <- which(r >= left_limit_r & r < r_peak)
    right_idx <- which(r > r_peak & r <= right_limit_r)
    
    # If we don't have points on both sides, skip shape metrics
    if (length(left_idx) < 2 || length(right_idx) < 2) {
      return(tibble(
        r_peak = r_peak,
        g_peak = g_peak,
        peak_prominence = NA_real_,
        peak_curvature = NA_real_,
        peak_width = NA_real_
      ))
    }
    
    left_min  <- min(g[left_idx], na.rm = TRUE)
    right_min <- min(g[right_idx], na.rm = TRUE)
    baseline  <- max(left_min, right_min)   # conservative baseline
    
    peak_prominence <- g_peak - baseline
    
    # Curvature: finite-difference 2nd derivative at idx (needs neighbours)
    if (idx <= 1 || idx >= length(g)) {
      peak_curvature <- NA_real_
    } else {
      gpp <- (g_s[idx + 1] - 2 * g_s[idx] + g_s[idx - 1]) / (dr^2)
      peak_curvature <- -gpp
    }
    
    # Width: full width at half prominence (FWHM-like)
    # Target height = baseline + 0.5 * prominence
    if (!is.finite(peak_prominence) || peak_prominence <= 0) {
      peak_width <- NA_real_
    } else {
      target <- baseline + 0.5 * peak_prominence
      
      # Find nearest crossings on left and right
      left_cross_candidates <- left_idx[g[left_idx] <= target]
      right_cross_candidates <- right_idx[g[right_idx] <= target]
      
      if (length(left_cross_candidates) == 0 || length(right_cross_candidates) == 0) {
        peak_width <- NA_real_
      } else {
        left_cross <- max(left_cross_candidates)
        right_cross <- min(right_cross_candidates)
        peak_width <- r[right_cross] - r[left_cross]
      }
    }
    
    tibble(
      r_peak = r_peak,
      g_peak = g_peak,
      peak_prominence = peak_prominence,
      peak_curvature = peak_curvature,
      peak_width = peak_width
    )
  })
  
  # Filter peaks using the NEW strength metric.
  # Replace your old min_prominence threshold with prominence-based threshold.
  out <- out %>%
    filter(!is.na(peak_prominence)) %>%
    filter(peak_prominence >= min_prominence)
  
  if (nrow(out) == 0) return(NULL)
  
  out
}

## Apply peak detection to all lek × date PCFs
peak_table <- pcf_curves %>%
  left_join(pcf_summary %>% select(lek_id, date, nn_median), by = c("lek_id", "date")) %>%
  group_by(lek_id, date) %>%
  group_modify(~{
    
    df <- .x %>% arrange(r)
    med_nnd <- unique(df$nn_median)
    
    if (length(med_nnd) != 1 || is.na(med_nnd)) return(tibble())
    
    peaks <- detect_peaks(df$r, df$g, med_nnd)
    
    if (is.null(peaks)) return(tibble())
    
    peaks
  }) %>% ungroup()

## Save peak table
peaks_out_file <- file.path(out_dir, "pcf_peak_table_ALL.csv")
write.csv(peak_table, peaks_out_file, row.names = FALSE)

message("Saved peak table to:", peaks_out_file)
