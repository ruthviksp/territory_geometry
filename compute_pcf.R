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
  
  # Inhomogeneous pair correlation
  g <- pcfinhom(X, lambda = lambda_hat, r = r_vals, correction = correction)
  
  # Extract translate correction and apply lower r cutoff
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

if (nrow(files_tbl) == 0) stop("No CSV files found across leks.")

## First pass: compute a reference median NND for EACH lek
nn_pass <- map_dfr(seq_len(nrow(files_tbl)), function(i) {
  
  row <- files_tbl[i, ]
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |>
    st_transform(32643) |> st_zm(drop = TRUE)
  
  df <- read.csv(row$csv_path)
  pts_sf <- st_as_sf(df, coords = c("pos_x", "pos_y"), crs = 32643)
  
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  X <- ppp(xy[, 1], xy[, 2], window = W)
  
  tibble(lek_id = row$lek_id, nn_median = median(nndist(X)))
})

lek_ref_nn <- nn_pass %>%
  group_by(lek_id) %>%
  summarise(ref_median_nn = median(nn_median, na.rm = TRUE), .groups = "drop")

files_tbl <- files_tbl %>% left_join(lek_ref_nn, by = "lek_id")

## Second pass: compute PCFs using a lek-specific r-grid
curve_list <- list()
summary_list <- list()

for (i in seq_len(nrow(files_tbl))) {
  
  row <- files_tbl[i, ]
  message("Processing ", row$lek_id, " : ", row$data_label, " (", i, "/", nrow(files_tbl), ")")
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |>
    st_transform(32643) |> st_zm(drop = TRUE)
  
  df <- read.csv(row$csv_path)
  n_pts <- nrow(df)
  lek_points <- st_as_sf(df, coords = c("pos_x", "pos_y"), crs = 32643)
  
  # Define lek-specific r range based on its reference median NND
  r_max <- r_max_mult * row$ref_median_nn
  r_vals <- seq(0, r_max, length.out = n_r)
  
  res <- compute_pcf(
    lek_polygon = lek_polygon,
    lek_points_sf = lek_points,
    r_vals = r_vals,
    r_min = r_min,
    correction = correction
  )
  
  curve_list[[i]] <- res$g_df %>%
    mutate(
      lek_id = row$lek_id,
      data_label = row$data_label,
      date = row$date,
      n_points = n_pts,
      r_max_used = r_max,
      ref_median_nn = row$ref_median_nn
    )
  
  summary_list[[i]] <- tibble(
    lek_id = row$lek_id,
    data_label = row$data_label,
    date = row$date,
    n_points = n_pts,
    nn_median = res$nn_median,
    bw_sigma = as.numeric(res$sigma),
    ref_median_nn = row$ref_median_nn,
    r_max_used = r_max
  )
}

pcf_curves <- bind_rows(curve_list)
pcf_summary <- bind_rows(summary_list)

## Write outputs
curves_out_file <- file.path(out_dir, "pcf_curve_ALL.csv")
summary_out_file <- file.path(out_dir, "pcf_summary_ALL.csv")
write.csv(pcf_curves, curves_out_file, row.names = FALSE)
write.csv(pcf_summary, summary_out_file, row.names = FALSE)

message("Saved curves to: ", curves_out_file)
message("Saved summary to: ", summary_out_file)

## Peak detection parameters
lower_nnd_mult <- 0.8
smooth_k <- 5
min_prominence <- 0.02
min_sep_mult <- 0.5

## Peak detection function
detect_peaks <- function(r, g, med_nnd) {
  
  # Rescale distance by this curve's median nearest-neighbour distance
  s <- r / med_nnd
  
  # Apply lower cutoff in rescaled units
  keep <- s >= lower_nnd_mult
  s <- s[keep]
  r <- r[keep]
  g <- g[keep]
  
  if (length(s) < 10) return(NULL)
  
  # Smooth g for robust peak/curvature estimation
  g_s <- zoo::rollmean(g, k = smooth_k, fill = NA, align = "center")
  
  # Local maxima on the smoothed curve
  is_peak <- g_s > dplyr::lag(g_s) & g_s > dplyr::lead(g_s)
  peak_idx <- which(is_peak)
  if (length(peak_idx) == 0) return(NULL)
  
  # Candidate peaks, filtered by minimum separation in rescaled units
  peaks <- tibble(idx = peak_idx, s_peak = s[peak_idx], r_peak = r[peak_idx],
                  g_peak = g[peak_idx]) %>% arrange(s_peak) %>% filter(c(TRUE, diff(s_peak) >= min_sep_mult))
  
  # Prominence window width in rescaled units
  win_half_width <- 0.75
  
  # Step size on rescaled axis for second-derivative curvature
  ds <- median(diff(s), na.rm = TRUE)
  if (!is.finite(ds) || ds <= 0) ds <- diff(range(s)) / (length(s) - 1)
  
  out <- purrr::pmap_dfr(peaks, function(idx, s_peak, r_peak, g_peak) {
    
    left_limit_s  <- s_peak - win_half_width
    right_limit_s <- s_peak + win_half_width
    
    left_idx  <- which(s >= left_limit_s & s < s_peak)
    right_idx <- which(s > s_peak & s <= right_limit_s)
    
    if (length(left_idx) < 2 || length(right_idx) < 2) {
      return(tibble(
        s_peak = s_peak,
        r_peak = r_peak,
        g_peak = g_peak,
        peak_prominence = NA_real_,
        peak_curvature = NA_real_
      ))
    }
    
    left_min  <- min(g[left_idx], na.rm = TRUE)
    right_min <- min(g[right_idx], na.rm = TRUE)
    baseline  <- max(left_min, right_min)
    
    peak_prominence <- g_peak - baseline
    
    # Curvature: negative second derivative of smoothed g(s)
    if (idx <= 1 || idx >= length(g_s)) {
      peak_curvature <- NA_real_
    } else {
      gpp <- (g_s[idx + 1] - 2 * g_s[idx] + g_s[idx - 1]) / (ds^2)
      peak_curvature <- -gpp
    }
    
    tibble(
      s_peak = s_peak,
      r_peak = r_peak,
      g_peak = g_peak,
      peak_prominence = peak_prominence,
      peak_curvature = peak_curvature
    )
  })
  
  out <- out %>% filter(!is.na(peak_prominence)) %>% filter(peak_prominence >= min_prominence)
  
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
    n_pts   <- unique(df$n_points)
    
    if (length(med_nnd) != 1 || is.na(med_nnd)) {
      return(tibble(
        s_peak = NA_real_,
        r_peak = NA_real_,
        g_peak = NA_real_,
        peak_prominence = NA_real_,
        peak_curvature = NA_real_,
        n_peaks = NA_integer_,
        n_points = n_pts
      ))
    }
    
    peaks <- detect_peaks(df$r, df$g, med_nnd)
    
    if (is.null(peaks) || nrow(peaks) == 0) {
      return(tibble(
        s_peak = NA_real_,
        r_peak = NA_real_,
        g_peak = NA_real_,
        peak_prominence = NA_real_,
        peak_curvature = NA_real_,
        n_peaks = 0L,
        n_points = n_pts
      ))
    }
    
    peaks %>% mutate(n_peaks = n(), n_points = n_pts)
  }) %>% ungroup()

## Save peak table
peaks_out_file <- file.path(out_dir, "pcf_peak_table_ALL.csv")
write.csv(peak_table, peaks_out_file, row.names = FALSE)

message("Saved peak table to: ", peaks_out_file)