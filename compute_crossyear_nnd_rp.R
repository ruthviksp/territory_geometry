rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(tidyverse)
library(rstudioapi)
library(patchwork)
library(cowplot)

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
setwd(dirname(getActiveDocumentContext()$path))
root_dir <- paste0(script_dir,"/rawdata")
# setwd(root_dir)

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
out_dir <- file.path(script_dir, "processed_data_rp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Analysis parameters
core_prob      <- 0.75   # probability mass captured by the KDE core region
kde_dimyx      <- 256    # KDE grid resolution
min_points_kde <- 5      # minimum points needed to estimate KDE
set.seed(123)

# Master table
files_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  
  cfg <- lek_configs[i, ]
  
  data_dirs <- list.dirs(file.path(root_dir, cfg$location),
                         recursive = FALSE, full.names = TRUE)
  data_dirs <- data_dirs[grepl("_COORDINATES$", basename(data_dirs))]
  
  map_dfr(data_dirs, function(d) {
    data_label <- sub("_COORDINATES$", "", basename(d))
    csv_path   <- list.files(d, pattern = paste0("_", cfg$suffix, "\\.csv$"),
                             full.names = TRUE)
    if (length(csv_path) == 0) return(NULL)
    tibble(lek_id = cfg$lek_id, location = cfg$location, suffix = cfg$suffix,
           shp_file = cfg$shp_file, data_label = data_label,
           date = parse_label_to_date(data_label, month_lookup), csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

## KDE HELPER FUNCTIONS

# Keep only points that fall strictly inside the lek polygon
clip_points_to_polygon <- function(pts_sf, polygon_sf) {
  pts_sf[lengths(st_within(pts_sf, polygon_sf)) > 0, ]
}

# Build a probability-normalised KDE grid and flag cells inside the core region
make_kde_grid <- function(pts_sf, lek_polygon, sigma, core_prob = 0.75, dimyx = 256) {
  W    <- as.owin(st_geometry(lek_polygon))
  xy   <- st_coordinates(pts_sf)
  X    <- ppp(x = xy[, 1], y = xy[, 2], window = W)
  dens <- density.ppp(X, sigma = sigma, edge = TRUE, dimyx = dimyx)
  
  kde_df     <- as.data.frame(dens) %>% setNames(c("x", "y", "z"))
  cell_area  <- dens$xstep * dens$ystep
  total_mass <- sum(kde_df$z, na.rm = TRUE) * cell_area
  kde_df     <- kde_df %>% mutate(p = z / total_mass)
  
  vals      <- sort(kde_df$p[is.finite(kde_df$p)], decreasing = TRUE)
  idx       <- which(cumsum(vals) * cell_area >= core_prob)[1]
  thresh    <- vals[idx]
  kde_df    <- kde_df %>% mutate(in_core = is.finite(p) & p >= thresh)
  
  list(kde_df = kde_df, contour_level = thresh,
       cell_area = cell_area, xstep = dens$xstep, ystep = dens$ystep)
}

# Return the subset of pts_sf whose nearest KDE grid cell lies inside the core mask
subset_points_to_kde_core <- function(pts_sf, kde_grid) {
  pts_xy     <- st_coordinates(pts_sf)
  core_cells <- kde_grid$kde_df %>% filter(in_core) %>% select(x, y)
  if (nrow(core_cells) == 0) return(pts_sf[0, ])
  
  x_vals    <- sort(unique(kde_grid$kde_df$x))
  y_vals    <- sort(unique(kde_grid$kde_df$y))
  nearest_x <- vapply(pts_xy[, 1], function(xx) x_vals[which.min(abs(x_vals - xx))], numeric(1))
  nearest_y <- vapply(pts_xy[, 2], function(yy) y_vals[which.min(abs(y_vals - yy))], numeric(1))
  
  inside <- tibble(x = nearest_x, y = nearest_y) %>%
    left_join(core_cells %>% mutate(in_core = TRUE), by = c("x", "y")) %>%
    mutate(in_core = !is.na(in_core)) %>%
    pull(in_core)
  pts_sf[inside, ]
}

# Logical flag: which pts_sf fall inside the KDE core mask from previous year
flag_points_in_kde_core <- function(pts_sf, kde_grid) {
  pts_xy     <- st_coordinates(pts_sf)
  core_cells <- kde_grid$kde_df %>% filter(in_core) %>% select(x, y)
  if (nrow(core_cells) == 0) return(rep(FALSE, nrow(pts_sf)))
  
  x_vals    <- sort(unique(kde_grid$kde_df$x))
  y_vals    <- sort(unique(kde_grid$kde_df$y))
  nearest_x <- vapply(pts_xy[, 1], function(xx) x_vals[which.min(abs(x_vals - xx))], numeric(1))
  nearest_y <- vapply(pts_xy[, 2], function(yy) y_vals[which.min(abs(y_vals - yy))], numeric(1))
  
  tibble(x = nearest_x, y = nearest_y) %>%
    left_join(core_cells %>% mutate(hit = TRUE), by = c("x", "y")) %>%
    mutate(hit = !is.na(hit)) %>%
    pull(hit)
}

subset_points_to_kde_core <- function(pts_sf, kde_grid) {
  pts_sf[flag_points_in_kde_core(pts_sf, kde_grid), ]
}

# Estimate a per-year KDE bandwidth with bw.diggle
get_kde_sigma <- function(pts_sf, lek_polygon) {
  W  <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  X  <- ppp(x = xy[, 1], y = xy[, 2], window = W)
  as.numeric(bw.diggle(X, hmax = diameter(W)))[1]
}

# Draw n_pts uniformly at random from within the KDE core cells - simulating csr in core polygon
sample_random_points_in_kde_core <- function(n_pts, kde_grid, crs_use) {
  core_cells <- kde_grid$kde_df %>% filter(in_core)
  samp       <- core_cells[sample(nrow(core_cells), n_pts, replace = TRUE), ]
  x_rand     <- runif(n_pts, samp$x - kde_grid$xstep / 2, samp$x + kde_grid$xstep / 2)
  y_rand     <- runif(n_pts, samp$y - kde_grid$ystep / 2, samp$y + kde_grid$ystep / 2)
  st_as_sf(tibble(x = x_rand, y = y_rand), coords = c("x", "y"), crs = crs_use)
}

# Convex hull owin of a KDE core (used for the grid null model)
core_to_owin <- function(kde_grid) {
  tryCatch(
    kde_grid$kde_df %>%
      filter(in_core) %>%
      st_as_sf(coords = c("x", "y"), crs = 32643) %>%
      st_union() %>% st_convex_hull() %>% as.owin(),
    error = function(e) NULL
  )
}


## PRE-COMPUTE LEK POLYGONS AND FIXED KDE BANDWIDTH PER LEK

# Creating lek-polygons of the outer boundary of each of leks
lek_polygons <- setNames(
  map(seq_len(nrow(lek_configs)), function(i) {
    cfg <- lek_configs[i, ]
    st_read(file.path(root_dir, cfg$location, cfg$shp_file), quiet = TRUE) %>%
      st_transform(32643) %>% st_zm(drop = TRUE) %>% st_make_valid()
  }),
  lek_configs$lek_id
)

# Median bandwidth across all years for each lek
sigma_fixed_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  cfg         <- lek_configs[i, ]
  lek_polygon <- lek_polygons[[cfg$lek_id]]
  files_lek   <- files_tbl %>% filter(lek_id == cfg$lek_id) %>% arrange(date)
  
  sigma_tbl <- map_dfr(seq_len(nrow(files_lek)), function(j) {
    row_j  <- files_lek[j, ]
    pts_sf <- read.csv(row_j$csv_path) %>%
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) %>%
      clip_points_to_polygon(lek_polygon)
    tibble(date = row_j$date, n_points = nrow(pts_sf),
           sigma_year = if (nrow(pts_sf) >= min_points_kde) get_kde_sigma(pts_sf, lek_polygon) else NA_real_)
  })
  
  tibble(lek_id = cfg$lek_id, sigma_fixed = median(sigma_tbl$sigma_year, na.rm = TRUE))
})


#### OBSERVED NNDs ####
## current points compared to previous year points
point_level_nnd <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  sub_tbl     <- files_tbl %>% filter(lek_id == lk) %>% arrange(date)
  lek_polygon <- lek_polygons[[lk]]
  sigma       <- sigma_fixed_tbl$sigma_fixed[sigma_fixed_tbl$lek_id == lk]
  W           <- as.owin(st_geometry(lek_polygon))
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    row_prev <- sub_tbl[i-1, ]
    row_now  <- sub_tbl[i, ]
    
    pts_prev <- read.csv(row_prev$csv_path) %>%
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) %>%
      clip_points_to_polygon(lek_polygon)
    pts_now  <- read.csv(row_now$csv_path) %>%
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) %>%
      clip_points_to_polygon(lek_polygon)
    
    if (nrow(pts_prev) < min_points_kde || nrow(pts_now) == 0) return(NULL)
    
    kde_prev     <- make_kde_grid(pts_prev, lek_polygon, sigma, core_prob, kde_dimyx)
    in_prev_core <- flag_points_in_kde_core(pts_now, kde_prev)
    peak         <- kde_prev$kde_df |> filter(is.finite(p)) |> slice_max(p, n = 1, with_ties = FALSE)
    cx <- peak$x; cy <- peak$y

    xy_prev <- st_coordinates(pts_prev)
    xy_now  <- st_coordinates(pts_now)
    X_prev  <- ppp(xy_prev[, 1], xy_prev[, 2], window = W)
    X_now   <- ppp(xy_now[, 1],  xy_now[, 2],  window = W)
    nn      <- nncross(X_now, X_prev)

    tibble(nnd_to_prev    = as.numeric(nn$dist),
           x = xy_now[, 1], y = xy_now[, 2],
           nn_x = xy_prev[nn$which, 1], nn_y = xy_prev[nn$which, 2],
           in_prev_core   = in_prev_core,
           dist_to_centre = sqrt((xy_now[, 1] - cx)^2 + (xy_now[, 2] - cy)^2),
           lek_id = lk,
           date_prev = row_prev$date, date_now = row_now$date,
           period = paste(format(row_prev$date, "%b%y"), format(row_now$date, "%b%y"), sep = "→"))
  })
})

nrow(point_level_nnd)
table(point_level_nnd$lek_id)


## CSR SIMULATION
## Random points drawn uniformly within the previous-year KDE core → nearest neighbour in pts_prev
point_level_nnd_sim <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  sub_tbl     <- files_tbl %>% filter(lek_id == lk) %>% arrange(date)
  lek_polygon <- lek_polygons[[lk]]
  sigma       <- sigma_fixed_tbl$sigma_fixed[sigma_fixed_tbl$lek_id == lk]
  W           <- as.owin(st_geometry(lek_polygon))
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    row_prev <- sub_tbl[i-1, ]
    row_now  <- sub_tbl[i, ]
    
    pts_prev <- read.csv(row_prev$csv_path) %>%
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) %>%
      clip_points_to_polygon(lek_polygon)
    pts_now  <- read.csv(row_now$csv_path) %>%
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) %>%
      clip_points_to_polygon(lek_polygon)
    
    if (nrow(pts_prev) < min_points_kde || nrow(pts_now) == 0) return(NULL)
    
    kde_prev     <- make_kde_grid(pts_prev, lek_polygon, sigma, core_prob, kde_dimyx)
    in_prev_core <- flag_points_in_kde_core(pts_now, kde_prev)
    peak         <- kde_prev$kde_df |> filter(is.finite(p)) |> slice_max(p, n = 1, with_ties = FALSE)
    cx <- peak$x; cy <- peak$y
    n_sim        <- nrow(pts_now)

    # Simulate n_sim random points inside the previous-year KDE core (matched to all current points)
    sim_pts <- sample_random_points_in_kde_core(n_sim, kde_prev, crs_use = st_crs(pts_prev))
    xy_prev <- st_coordinates(pts_prev)
    xy_sim  <- st_coordinates(sim_pts)
    X_prev  <- ppp(xy_prev[, 1], xy_prev[, 2], window = W)
    X_sim   <- ppp(xy_sim[, 1],  xy_sim[, 2],  window = W)
    nn      <- nncross(X_sim, X_prev)

    tibble(nnd_to_prev    = as.numeric(nn$dist),
           x = xy_sim[, 1], y = xy_sim[, 2],
           nn_x = xy_prev[nn$which, 1], nn_y = xy_prev[nn$which, 2],
           in_prev_core   = in_prev_core,
           dist_to_centre = sqrt((xy_sim[, 1] - cx)^2 + (xy_sim[, 2] - cy)^2),
           lek_id = lk,
           date_prev = row_prev$date, date_now = row_now$date,
           period = paste(format(row_prev$date, "%b%y"), format(row_now$date, "%b%y"), sep = "→"))
  })
})

write_csv(point_level_nnd, file.path(out_dir, "crossyear_nnd_obs.csv"))
write_csv(point_level_nnd_sim, file.path(out_dir, "crossyear_nnd_sim.csv"))



#### KDE core visualisation across years ####

walk(lek_configs$lek_id, function(lk) {
  cat("Building KDE diagnostic for", lk, "...\n")
  
  files_lek   <- files_tbl %>% filter(lek_id == lk) %>% arrange(date)
  lek_polygon <- lek_polygons[[lk]]
  sigma       <- sigma_fixed_tbl$sigma_fixed[sigma_fixed_tbl$lek_id == lk]
  
  year_data <- keep(
    map(seq_len(nrow(files_lek)), function(i) {
      row <- files_lek[i, ]
      
      pts_clipped <- read.csv(row$csv_path) %>%
        st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) %>%
        clip_points_to_polygon(lek_polygon)
      
      if (nrow(pts_clipped) < min_points_kde) return(NULL)
      
      kde_grid <- make_kde_grid(pts_clipped, lek_polygon, sigma, core_prob, kde_dimyx)
      
      # Classify each clipped point as inside / outside the KDE core
      pts_xy <- st_coordinates(pts_clipped)
      x_vals <- sort(unique(kde_grid$kde_df$x))
      y_vals <- sort(unique(kde_grid$kde_df$y))
      nx     <- vapply(pts_xy[, 1], function(xx) x_vals[which.min(abs(x_vals - xx))], numeric(1))
      ny     <- vapply(pts_xy[, 2], function(yy) y_vals[which.min(abs(y_vals - yy))], numeric(1))
      in_core_pts <- tibble(x = nx, y = ny) %>%
        left_join(kde_grid$kde_df %>% filter(in_core) %>%
                    select(x, y) %>% mutate(in_core = TRUE),
                  by = c("x", "y")) %>%
        mutate(in_core = !is.na(in_core)) %>% pull(in_core)
      
      yr_label <- paste0(format(row$date, "%b %Y"), "\nn=", nrow(pts_clipped),
                         " (", sum(in_core_pts), " core)")
      
      list(
        kde = kde_grid$kde_df %>% mutate(yr_label = yr_label),
        pts = tibble(x = pts_xy[, 1], y = pts_xy[, 2],
                     in_core = in_core_pts, yr_label = yr_label)
      )
    }),
    Negate(is.null)
  )
  
  all_kde <- bind_rows(map(year_data, "kde"))
  all_pts <- bind_rows(map(year_data, "pts"))
  
  # Preserve chronological order of facets
  yr_levels <- unique(all_kde$yr_label)
  all_kde$yr_label <- factor(all_kde$yr_label, levels = yr_levels)
  all_pts$yr_label <- factor(all_pts$yr_label, levels = yr_levels)
  
  p <- ggplot() +
    geom_tile(data = all_kde,
              aes(x = x, y = y, fill = in_core), alpha = 0.65) +
    geom_sf(data = lek_polygon, fill = NA, colour = "black", linewidth = 0.8) +
    geom_point(data = all_pts %>% filter(!in_core),
               aes(x = x, y = y), colour = "#c0392b", shape = 4,
               size = 1.5, stroke = 0.8) +
    geom_point(data = all_pts %>% filter(in_core),
               aes(x = x, y = y), colour = "#27ae60", shape = 16, size = 1.8) +
    facet_wrap(~ yr_label, ncol = 5) +
    scale_fill_manual(values = c("TRUE" = "#2980b9", "FALSE" = "grey88"),
                      labels = c("TRUE" = "KDE core", "FALSE" = "Outside core"),
                      name   = "KDE region") +
    coord_sf(datum = sf::st_crs(32643)) +
    labs(title    = paste("KDE core diagnostic —", lk),
         subtitle = paste0("Fixed sigma = ", round(sigma, 1),
                           " m  |  Core = ", round(core_prob * 100),
                           "% KDE mass  |  green = in core, red × = outside"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          axis.text       = element_text(size = 5),
          strip.text      = element_text(size = 6.5, face = "bold"),
          plot.title      = element_text(size = 11, face = "bold"),
          plot.subtitle   = element_text(size = 8),
          panel.grid      = element_blank())
  
  print(p)
})

