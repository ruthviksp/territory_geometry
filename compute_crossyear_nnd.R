rm(list = ls())

## Load libraries
library(sf)
library(rstudioapi)
library(purrr)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(spatstat.geom)
library(spatstat.explore)

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

## Output folder
out_dir <- file.path(script_dir, "processed_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Lek configurations
lek_configs <- tibble(
  lek_id   = c("Tal Chhapar", "Velavadar Lek 1", "Velavadar Lek 2"),
  location = c("TalChhapar", "Velavadar", "Velavadar"),
  suffix   = c("TC", "LEK1", "LEK2"),
  shp_file = c("TalChhapar_Area.shp", "Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp")
)

## Build master table of all files across all leks
files_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  
  cfg <- lek_configs %>% slice(i)
  
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

## ANALYSIS PARAMETERS
## KDE settings
core_prob <- 0.75
kde_dimyx <- 256
min_points_kde <- 10

## Null simulation settings
n_sims <- 1000
set.seed(123)

## HELPER FUNCTIONS
## Consider points within the lek polygon only
clip_points_to_polygon <- function(pts_sf, polygon_sf) {
  inside <- lengths(st_within(pts_sf, polygon_sf)) > 0
  pts_sf[inside, ]
}

## Compute a KDE of the points and make a mask that defines the core region
make_kde_grid <- function(pts_sf, lek_polygon, sigma, core_prob = 0.75, dimyx = 256) {
  
  # Convert the polygon to a spatstat window and extract point coordinates
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  
  # Build the point pattern and estimate the KDE on a grid
  X <- ppp(x = xy[, 1], y = xy[, 2], window = W)
  dens <- density.ppp(X, sigma = sigma, edge = TRUE, dimyx = dimyx)
  
  # Build the point pattern and estimate the KDE on a grid
  kde_df <- as.data.frame(dens)
  names(kde_df) <- c("x", "y", "z")
  
  # Compute KDE mass and normalize so the KDE values so the mass sums up to 1
  cell_area <- dens$xstep * dens$ystep
  total_mass <- sum(kde_df$z, na.rm = TRUE) * cell_area
  kde_df <- kde_df %>% mutate(p = z / total_mass)
  
  # Sort KDE grid cell densities
  vals <- kde_df$p[is.finite(kde_df$p)]
  vals <- sort(vals, decreasing = TRUE)
  
  # Find the density threshold needed to capture the chosen core probability mass
  mass_cum <- cumsum(vals) * cell_area
  idx <- which(mass_cum >= core_prob)[1]
  contour_level <- vals[idx]
  
  # Mark grid cells as inside or outside the KDE core mask
  kde_df <- kde_df %>% mutate(core_prob = core_prob, sigma_used = sigma, contour_level = contour_level, 
                              in_core = is.finite(p) & (p >= contour_level))
  
  list(kde_df = kde_df, contour_level = contour_level, cell_area = cell_area, xstep = dens$xstep, ystep = dens$ystep)
}

## Check which points fall inside the core region (KDE mask)
subset_points_to_kde_core <- function(pts_sf, kde_grid) {
  
  # Extract point coordinates and only keep points that are within the KDE core
  pts_xy <- st_coordinates(pts_sf)
  core_cells <- kde_grid$kde_df %>% filter(in_core) %>% select(x, y)
  if (nrow(core_cells) == 0) return(pts_sf[0, ])
  
  # Match each point to nearest grid-cell centre in x and y separately
  x_vals <- sort(unique(kde_grid$kde_df$x))
  y_vals <- sort(unique(kde_grid$kde_df$y))
  
  # Get grid cell centres along each axis
  nearest_x <- vapply(pts_xy[, 1], function(xx) x_vals[which.min(abs(x_vals - xx))], numeric(1))
  nearest_y <- vapply(pts_xy[, 2], function(yy) y_vals[which.min(abs(y_vals - yy))], numeric(1))
  
  pts_key <- tibble(x = nearest_x, y = nearest_y)
  core_key <- core_cells %>% mutate(in_core = TRUE)
  
  # Mark and keep points whose matched grid cell lies inside the KDE core
  inside <- pts_key %>% left_join(core_key, by = c("x", "y")) %>%
    mutate(in_core = ifelse(is.na(in_core), FALSE, in_core)) %>% pull(in_core)
  
  pts_sf[inside, ]
}

## Compute distance of points from the KDE mode
distance_to_kde_mode <- function(pts_sf, kde_grid, crs_use) {
  mode_cell <- kde_grid$kde_df %>% filter(is.finite(p)) %>% slice_max(order_by = p, n = 1, with_ties = FALSE)
  mode_pt <- st_as_sf(mode_cell, coords = c("x", "y"), crs = crs_use)
  as.numeric(st_distance(pts_sf, mode_pt))
}

## Estimate KDE bandwidth for a given point pattern
get_kde_sigma <- function(pts_sf, lek_polygon) {
  
  # Convert polygon and points to spatstat objects and create a point pattern object
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  
  X <- ppp(x = xy[, 1], y = xy[, 2], window = W)
  
  # Estimate the smoothing bandwidth for the KDE
  sigma <- bw.diggle(X, hmax = diameter(W))
  sigma <- as.numeric(sigma)[1]
  
  return(sigma)
}

## Compute nearest-neighbour distance between points
nnd <- function(from_pts, to_pts) {
  dmat <- st_distance(from_pts, to_pts)
  dmin <- apply(dmat, 1, min)
  
  as.numeric(dmin)
}

## Translate points by a fixed distance in a random direction
translate_points_random <- function(pts_sf, shift_dist) {
  
  # Draw a random direction and shift all points by the same offset
  theta <- runif(1, 0, 2 * pi)
  dx <- shift_dist * cos(theta)
  dy <- shift_dist * sin(theta)
  
  pts_xy <- st_coordinates(pts_sf)
  pts_shift <- tibble(x = pts_xy[, 1] + dx, y = pts_xy[, 2] + dy)
  
  st_as_sf(pts_shift, coords = c("x", "y"), crs = st_crs(pts_sf))
}

## Rotate points by a random angle around their centroid
rotate_points_random <- function(pts_sf, angle_max = pi / 6) {
  
  # Draw a random rotation angle and rotate all points around the centroid
  theta <- runif(1, -angle_max, angle_max)
  pts_xy <- st_coordinates(pts_sf)
  ctr <- colMeans(pts_xy)
  pts_ctr <- sweep(pts_xy, 2, ctr, "-")
  
  rot_mat <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE)
  pts_rot <- pts_ctr %*% t(rot_mat)
  pts_rot <- sweep(pts_rot, 2, ctr, "+")
  pts_rot <- tibble(x = pts_rot[, 1], y = pts_rot[, 2])
  
  st_as_sf(pts_rot, coords = c("x", "y"), crs = st_crs(pts_sf))
}

## Translate and rotate points
transform_points_random <- function(pts_sf, shift_dist, angle_max = pi / 6) {
  pts_trans <- translate_points_random(pts_sf, shift_dist = shift_dist)
  rotate_points_random(pts_trans, angle_max = angle_max)
}

## Simulate transformed previous-year points and compute nearest-neighbour distances
simulate_transform_crossyear_nnd <- function(prev_pts, curr_pts, n_sims = 999) {

  # Compute the shift distance as half of the median nearest-neighbour distance within the previous image
  prev_dmat <- st_distance(prev_pts, prev_pts)
  diag(prev_dmat) <- Inf
  prev_nnd <- apply(prev_dmat, 1, min)
  shift_dist <- median(as.numeric(prev_nnd)) / 2

  # Repeat the translate-rotate randomisation and store point-level NNDs from each simulation
  sim_pointwise <- map_dfr(seq_len(n_sims), function(s) {
    prev_transform <- transform_points_random(prev_pts, shift_dist = shift_dist, angle_max = pi / 6)
    transform_nnd <- nnd(curr_pts, prev_transform)

    tibble(sim = s, point_id = seq_along(transform_nnd),
           nnd_to_prev = transform_nnd)
  })

  # Summarise point-level NNDs from each simulation
  sim_summary <- sim_pointwise %>%
    group_by(sim) %>%
    summarise(mean_nnd_transform_rand = mean(nnd_to_prev),
              median_nnd_transform_rand = median(nnd_to_prev),
              .groups = "drop")

  list(summary = sim_summary, pointwise = sim_pointwise)
}

## MAIN
## Read and store information about the polygons for each lek
lek_polygons <- map(seq_len(nrow(lek_configs)), function(i) {
  cfg <- lek_configs %>% slice(i)
  poly_sf <- st_read(file.path(root_dir, cfg$location, cfg$shp_file),
                     quiet = TRUE) %>% st_transform(32643) %>% st_zm(drop = TRUE) %>% st_make_valid()
})

names(lek_polygons) <- lek_configs$lek_id

## Compute a single representative KDE bandwidth per lek
sigma_fixed_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  # Extract config for this lek and get all files ordered by date
  cfg <- lek_configs %>% slice(i)
  files_lek <- files_tbl %>% filter(lek_id == cfg$lek_id) %>% arrange(date)
  
  lek_polygon <- lek_polygons[[cfg$lek_id]]
  
  # Loop over dates and estimate the bandwidth value for each time point
  sigma_tbl <- map_dfr(seq_len(nrow(files_lek)), function(j) {
    
    row_j <- files_lek[j, ]
    
    # Read the point coordinates and only keep ones inside the lek polygon
    pts_sf <- read.csv(row_j$csv_path) %>% st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    pts_sf <- clip_points_to_polygon(pts_sf, lek_polygon)
    
    # Store the number of points and the KDE bandwidth 
    tibble(date = row_j$date, n_points = nrow(pts_sf),
           sigma_year = if (nrow(pts_sf) >= min_points_kde) get_kde_sigma(pts_sf, lek_polygon) else NA_real_)
  })
  
  # Take the median bandwidth across years as the fixed value for each lek
  tibble(lek_id = cfg$lek_id, sigma_fixed = median(sigma_tbl$sigma_year, na.rm = TRUE))
})

## Read and store all points (lek id x year)
pts_by_lek_year <- map(seq_len(nrow(files_tbl)), function(i) {
  
  row_i <- files_tbl[i, ]
  lek_polygon <- lek_polygons[[row_i$lek_id]]
  
  pts_sf <- read.csv(row_i$csv_path) %>% st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
  pts_sf <- clip_points_to_polygon(pts_sf, lek_polygon)
  pts_sf %>% mutate(lek_id = row_i$lek_id, date = row_i$date)
})

names(pts_by_lek_year) <- paste(files_tbl$lek_id, files_tbl$date, sep = "__")

## Initialise output containers
summary_list <- list()
pointwise_list <- list()
sim_list <- list()
pointwise_randomisation_list <- list()

## Split the file table by lek
files_by_lek <- split(files_tbl, files_tbl$lek_id)

## For each lek
for (lek in names(files_by_lek)) {
  
  # Get all files ordered by date
  files_lek <- files_by_lek[[lek]] %>% arrange(date)
  
  # Retreive the lek polygon and the KDE bandwidth for this lek
  lek_polygon <- lek_polygons[[lek]]
  sigma_fixed <- as.numeric(sigma_fixed_tbl[sigma_fixed_tbl$lek_id == lek,'sigma_fixed'])
  
  # Loop over consecutive dates
  for (i in 2:nrow(files_lek)) {
    
    # Extract points and metadata for previous and current year
    row_prev <- files_lek[i-1,]
    row_curr <- files_lek[i,]
    
    date_prev <- row_prev$date
    date_curr <- row_curr$date
    
    pts_prev <- pts_by_lek_year[[paste(row_prev$lek_id, row_prev$date, sep = "__")]]
    pts_curr <- pts_by_lek_year[[paste(row_curr$lek_id, row_curr$date, sep = "__")]]
    
    # Compute KDE grid from t-1
    kde_prev <- make_kde_grid(pts_sf = pts_prev, lek_polygon = lek_polygon, sigma = sigma_fixed, core_prob = core_prob, dimyx = kde_dimyx)
    
    # Keep only current-year points inside previous-year KDE core
    pts_curr_in_core <- subset_points_to_kde_core(pts_curr, kde_prev)
    
    # Compute nearest-neighbour distance from current to previous points
    obs_nnd <- nnd(pts_curr_in_core, pts_prev)
    dist_from_centre <- distance_to_kde_mode(pts_curr_in_core, kde_prev, crs_use = st_crs(pts_prev))
    
    # Simulate transformed previous-year points and compute NNDs
    sim_out <- simulate_transform_crossyear_nnd(prev_pts = pts_prev, curr_pts = pts_curr_in_core, n_sims = n_sims)
    sim_tbl <- sim_out$summary
    sim_pointwise_tbl <- sim_out$pointwise
    
    # Observed summary statistics
    obs_mean <- mean(obs_nnd)
    obs_median <- median(obs_nnd)
    
    # Qunatiles of simulated transformed mean NNDs
    rand_mean_transform_q025 <- quantile(sim_tbl$mean_nnd_transform_rand, 0.025, na.rm = TRUE)
    rand_mean_transform_q25  <- quantile(sim_tbl$mean_nnd_transform_rand, 0.25, na.rm = TRUE)
    rand_mean_transform_q50  <- quantile(sim_tbl$mean_nnd_transform_rand, 0.50, na.rm = TRUE)
    rand_mean_transform_q75  <- quantile(sim_tbl$mean_nnd_transform_rand, 0.75, na.rm = TRUE)
    rand_mean_transform_q975 <- quantile(sim_tbl$mean_nnd_transform_rand, 0.975, na.rm = TRUE)

    # Quantiles of simulated transformed median NNDs
    rand_median_transform_q025 <- quantile(sim_tbl$median_nnd_transform_rand, 0.025, na.rm = TRUE)
    rand_median_transform_q25  <- quantile(sim_tbl$median_nnd_transform_rand, 0.25, na.rm = TRUE)
    rand_median_transform_q50  <- quantile(sim_tbl$median_nnd_transform_rand, 0.50, na.rm = TRUE)
    rand_median_transform_q75  <- quantile(sim_tbl$median_nnd_transform_rand, 0.75, na.rm = TRUE)
    rand_median_transform_q975 <- quantile(sim_tbl$median_nnd_transform_rand, 0.975, na.rm = TRUE)

    # Store summary statistics for this time step
    summary_list[[length(summary_list) + 1]] <- tibble(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
                                                       n_prev = nrow(pts_prev), n_curr = nrow(pts_curr), 
                                                       n_curr_in_prev_core = nrow(pts_curr_in_core),
                                                       obs_mean_nnd = obs_mean, obs_median_nnd = obs_median,
                                                       rand_mean_transform_nnd_q025 = as.numeric(rand_mean_transform_q025),
                                                       rand_mean_transform_nnd_q25 = as.numeric(rand_mean_transform_q25),
                                                       rand_mean_transform_nnd_q50 = as.numeric(rand_mean_transform_q50),
                                                       rand_mean_transform_nnd_q75 = as.numeric(rand_mean_transform_q75),
                                                       rand_mean_transform_nnd_q975 = as.numeric(rand_mean_transform_q975),
                                                       rand_median_transform_nnd_q025 = as.numeric(rand_median_transform_q025),
                                                       rand_median_transform_nnd_q25 = as.numeric(rand_median_transform_q25),
                                                       rand_median_transform_nnd_q50 = as.numeric(rand_median_transform_q50),
                                                       rand_median_transform_nnd_q75 = as.numeric(rand_median_transform_q75),
                                                       rand_median_transform_nnd_q975 = as.numeric(rand_median_transform_q975))
    
    # Store point-level NNDs for current points within the KDE core from previous year's points
    pointwise_list[[length(pointwise_list) + 1]] <- pts_curr_in_core %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr, point_id = seq_along(obs_nnd),
             nnd_to_prev = obs_nnd, dist_from_centre = dist_from_centre) %>% st_drop_geometry()
    
    # Store simulation output for the same transition
    sim_list[[length(sim_list) + 1]] <- sim_tbl %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr)

    # Store point-level simulation output for the same transition
    pointwise_randomisation_list[[length(pointwise_randomisation_list) + 1]] <- sim_pointwise_tbl %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
             dist_from_centre = dist_from_centre[point_id])
  }
}

summary_tbl <- bind_rows(summary_list) %>% arrange(lek_id, date_curr)
pointwise_tbl <- bind_rows(pointwise_list) %>% arrange(lek_id, date_curr)
sim_tbl_all <- bind_rows(sim_list) %>% arrange(lek_id, date_curr, sim)
pointwise_randomisation_tbl <- bind_rows(pointwise_randomisation_list) %>% arrange(lek_id, date_curr, sim, point_id)

## Export dataframes
summary_tbl <- summary_tbl %>%
  transmute(lek_id, date_prev, date_now = date_curr, n_prev = n_prev, n_curr = n_curr, n_curr_in_prev_core = n_curr_in_prev_core,
            mean_crossyear_nnd = obs_mean_nnd, median_crossyear_nnd = obs_median_nnd)

simulation_tbl <- sim_tbl_all %>%
  transmute(lek_id, date_prev, date_now = date_curr, sim, 
            mean_crossyear_nnd_transform_rand = mean_nnd_transform_rand,
            median_crossyear_nnd_transform_rand = median_nnd_transform_rand)

pointwise_randomisation_tbl <- pointwise_randomisation_tbl %>%
  transmute(lek_id, date_prev, date_now = date_curr, sim, point_id, nnd_to_prev, dist_from_centre)

crossyear_nnd_tbl <- simulation_tbl %>%
  left_join(summary_tbl, by = c("lek_id", "date_prev", "date_now")) %>%
  relocate(lek_id, date_prev, date_now, n_curr_in_prev_core, sim, mean_crossyear_nnd, median_crossyear_nnd,
           mean_crossyear_nnd_transform_rand, median_crossyear_nnd_transform_rand)

write_csv(summary_tbl, file.path(out_dir, "crossyear_nnd_summary.csv"))
write_csv(pointwise_tbl, file.path(out_dir, "crossyear_nnd_pointwise.csv"))
write_csv(pointwise_randomisation_tbl, file.path(out_dir, "crossyear_nnd_pointwise_randomisations.csv"))
write_csv(crossyear_nnd_tbl, file.path(out_dir, "crossyear_nnd_with_randomisations.csv"))
