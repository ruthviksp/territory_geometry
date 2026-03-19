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
out_dir <- file.path(script_dir, "processed_data", "crossyear_nnd_gridmask_all_leks")
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
min_points_kde <- 5

## Null simulation settings
n_sims <- 9
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

## Sample random points uniformly within the core region (KDE mask)
sample_random_points_in_kde_core <- function(n_pts, kde_grid, crs_use) {
  
  # Only keep cells within the KDE core mask
  core_cells <- kde_grid$kde_df %>% filter(in_core)
  
  # Sample core cells with replacement
  samp_idx <- sample(seq_len(nrow(core_cells)), size = n_pts, replace = TRUE)
  samp_cells <- core_cells[samp_idx, ]
  
  # Draw a random point within the sampled cell
  x_rand <- runif(n_pts, min = samp_cells$x - kde_grid$xstep / 2, max = samp_cells$x + kde_grid$xstep / 2)
  y_rand <- runif(n_pts, min = samp_cells$y - kde_grid$ystep / 2, max = samp_cells$y + kde_grid$ystep / 2)
  
  st_as_sf(tibble(x = x_rand, y = y_rand), coords = c("x", "y"), crs = crs_use)
}

## Compute nearest-neighbour distance between points
nnd <- function(from_pts, to_pts) {
  dmat <- st_distance(from_pts, to_pts)
  dmin <- apply(dmat, 1, min)
  
  as.numeric(dmin)
}

## Simulate random points within the core region from the previous year and compute nearest-neighbour distances
simulate_random_crossyear_nnd <- function(n_pts, kde_grid, prev_pts, n_sims = 999, crs_use) {
  
  # Repeat the randomisation and summarise the NNDs from each simulation
  map_dfr(seq_len(n_sims), function(s) {
    rand_pts <- sample_random_points_in_kde_core(n_pts = n_pts, kde_grid = kde_grid, crs_use = crs_use)
    rand_nnd <- nnd(rand_pts, prev_pts)
    
    tibble(sim = s, mean_nnd_rand = mean(rand_nnd), median_nnd_rand = median(rand_nnd))
  })
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
    
    # Simulate random points within the same KDE core and compute NNDs
    sim_tbl <- simulate_random_crossyear_nnd(n_pts = nrow(pts_curr_in_core), kde_grid = kde_prev,
                                             prev_pts = pts_prev, n_sims = n_sims, crs_use = st_crs(pts_prev))
    
    # Observed summary statistics
    obs_mean <- mean(obs_nnd)
    obs_median <- median(obs_nnd)
    
    # Qunatiles of simulated mean NNDs
    rand_mean_q025 <- quantile(sim_tbl$mean_nnd_rand, 0.025, na.rm = TRUE)
    rand_mean_q25  <- quantile(sim_tbl$mean_nnd_rand, 0.25, na.rm = TRUE)
    rand_mean_q50  <- quantile(sim_tbl$mean_nnd_rand, 0.50, na.rm = TRUE)
    rand_mean_q75  <- quantile(sim_tbl$mean_nnd_rand, 0.75, na.rm = TRUE)
    rand_mean_q975 <- quantile(sim_tbl$mean_nnd_rand, 0.975, na.rm = TRUE)
    
    # Quantiles of simulated median NNDs
    rand_median_q025 <- quantile(sim_tbl$median_nnd_rand, 0.025, na.rm = TRUE)
    rand_median_q25  <- quantile(sim_tbl$median_nnd_rand, 0.25, na.rm = TRUE)
    rand_median_q50  <- quantile(sim_tbl$median_nnd_rand, 0.50, na.rm = TRUE)
    rand_median_q75  <- quantile(sim_tbl$median_nnd_rand, 0.75, na.rm = TRUE)
    rand_median_q975 <- quantile(sim_tbl$median_nnd_rand, 0.975, na.rm = TRUE)
    
    # Store summary statistics for this time step
    summary_list[[length(summary_list) + 1]] <- tibble(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
                                                       n_prev = nrow(pts_prev), n_curr = nrow(pts_curr), 
                                                       n_curr_in_prev_core = nrow(pts_curr_in_core),
                                                       obs_mean_nnd = obs_mean, obs_median_nnd = obs_median,
                                                       rand_mean_nnd_q025 = as.numeric(rand_mean_q025),
                                                       rand_mean_nnd_q25 = as.numeric(rand_mean_q25),
                                                       rand_mean_nnd_q50 = as.numeric(rand_mean_q50),
                                                       rand_mean_nnd_q75 = as.numeric(rand_mean_q75),
                                                       rand_mean_nnd_q975 = as.numeric(rand_mean_q975),
                                                       rand_median_nnd_q025 = as.numeric(rand_median_q025),
                                                       rand_median_nnd_q25 = as.numeric(rand_median_q25),
                                                       rand_median_nnd_q50 = as.numeric(rand_median_q50),
                                                       rand_median_nnd_q75 = as.numeric(rand_median_q75),
                                                       rand_median_nnd_q975 = as.numeric(rand_median_q975))
    
    # Store point-level NNDs for current points within the KDE core from previous year's points
    pointwise_list[[length(pointwise_list) + 1]] <- pts_curr_in_core %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr, nnd_to_prev = obs_nnd) %>% st_drop_geometry()
    
    # Store simulation output for the same transition
    sim_list[[length(sim_list) + 1]] <- sim_tbl %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr)
  }
}

summary_tbl <- bind_rows(summary_list) %>% arrange(lek_id, date_curr)
pointwise_tbl <- bind_rows(pointwise_list) %>% arrange(lek_id, date_curr)
sim_tbl_all <- bind_rows(sim_list) %>% arrange(lek_id, date_curr, sim)

## Export dataframes
summary_tbl <- summary_tbl %>%
  transmute(lek_id, date_prev, date_now = date_curr, 
            mean_crossyear_nnd = obs_mean_nnd, median_crossyear_nnd = obs_median_nnd)

simulation_tbl <- sim_tbl_all %>%
  transmute(lek_id, date_prev, date_now = date_curr, sim, 
            mean_crossyear_nnd_rand = mean_nnd_rand, median_crossyear_nnd_rand = median_nnd_rand)

crossyear_nnd_tbl <- simulation_tbl %>%
  left_join(summary_tbl, by = c("lek_id", "date_prev", "date_now")) %>%
  relocate(lek_id, date_prev, date_now, sim, mean_crossyear_nnd, median_crossyear_nnd,
           mean_crossyear_nnd_rand, median_crossyear_nnd_rand)

write_csv(summary_tbl, file.path(out_dir, "crossyear_nnd_summary.csv"))
write_csv(pointwise_tbl, file.path(out_dir, "crossyear_nnd_pointwise.csv"))
write_csv(crossyear_nnd_tbl, file.path(out_dir, "crossyear_nnd_export_for_plotting.csv"))

## ---------- Plot mean NND over time ----------

p_mean <- ggplot(summary_tbl, aes(x = date_curr)) +
  geom_ribbon(
    aes(ymin = rand_mean_nnd_q025, ymax = rand_mean_nnd_q975),
    fill = "grey85"
  ) +
  geom_ribbon(
    aes(ymin = rand_mean_nnd_q25, ymax = rand_mean_nnd_q75),
    fill = "grey70"
  ) +
  geom_line(
    aes(y = rand_mean_nnd_q50),
    color = "grey30",
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_line(
    aes(y = obs_mean_nnd),
    color = "#0072B2",
    linewidth = 1.2
  ) +
  geom_point(
    aes(y = obs_mean_nnd),
    color = "#0072B2",
    size = 2.8
  ) +
  labs(
    title = "Cross-year mean NND across leks",
    subtitle = paste0(
      "Observed year t to t-1 nearest-neighbour distance after restricting year t points to the previous-year KDE core grid mask (core_prob = ",
      core_prob, "). Grey ribbons show the null from random placement within that same mask."
    ),
    x = "Year transition",
    y = "Mean cross-year NND (m)"
  ) +
  scale_x_date(
    breaks = sort(unique(summary_tbl$date_curr)),
    labels = format(sort(unique(summary_tbl$date_curr)), "%Y")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(~ lek_id, scales = "free_x")

p_mean


## ---------- Plot median NND over time ----------

p_median <- ggplot(summary_tbl, aes(x = date_curr)) +
  geom_ribbon(
    aes(ymin = rand_median_nnd_q025, ymax = rand_median_nnd_q975),
    fill = "grey85"
  ) +
  geom_ribbon(
    aes(ymin = rand_median_nnd_q25, ymax = rand_median_nnd_q75),
    fill = "grey70"
  ) +
  geom_line(
    aes(y = rand_median_nnd_q50),
    color = "grey30",
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_line(
    aes(y = obs_median_nnd),
    color = "#D55E00",
    linewidth = 1.2
  ) +
  geom_point(
    aes(y = obs_median_nnd),
    color = "#D55E00",
    size = 2.8
  ) +
  labs(
    title = "Cross-year median NND across leks",
    subtitle = paste0(
      "Observed year t to t-1 nearest-neighbour distance after restricting year t points to the previous-year KDE core grid mask (core_prob = ",
      core_prob, "). Grey ribbons show the null from random placement within that same mask."
    ),
    x = "Year transition",
    y = "Median cross-year NND (m)"
  ) +
  scale_x_date(
    breaks = sort(unique(summary_tbl$date_curr)),
    labels = format(sort(unique(summary_tbl$date_curr)), "%Y")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(~ lek_id, scales = "free_x")

p_median

message("Saved grid-mask cross-year NND outputs and plots to: ", out_dir)
