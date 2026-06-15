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
library(lubridate)
library(deldir)
library(mgcv)

## Get directory where this script lives
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(getActiveDocumentContext()$path))
root_dir <- paste0(script_dir,"/rawdata")

## Lookup table to convert file names to dates
month_lookup <- c(JAN = 1, FEB = 2, MAR = 3, APR = 4,
                  MAY = 5, JUN = 6, JUL = 7, AUG = 8,
                  SEP = 9, OCT = 10, NOV = 11, DEC = 12)

parse_label_to_date <- function(data_label, month_lookup) {
  yy <- as.integer(substr(data_label, 1, 2))
  mm <- unname(month_lookup[substr(data_label, 3, 5)])
  as.Date(sprintf("20%02d-%02d-01", yy, mm))
}

## Lek configuration table
lek_configs <- tibble(
  lek_id   = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
  location = c("Velavadar", "Velavadar", "TalChhapar"),
  suffix   = c("LEK1", "LEK2", "TC"),
  shp_file = c("Velavadar_Lek1_Area.shp",
               "Velavadar_Lek2_Area.shp",
               "TalChhapar_Area.shp")
)

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


## Output folder
out_dir <- file.path(script_dir, "processed_data_plots_rp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Reading data
point_level_nnd <- read.csv("./processed_data_rp/crossyear_nnd_obs.csv")
point_level_nnd_sim <- read.csv("./processed_data_rp/crossyear_nnd_sim.csv")

nnd_combined <- bind_rows(
  point_level_nnd     %>% mutate(source = "observed"),
  point_level_nnd_sim %>% mutate(source = "CSR")
) 

#selecting only core territories
nnd_combined <- nnd_combined[nnd_combined$in_prev_core,]

cols_to_change <- c("date_prev", "date_now")
nnd_combined[cols_to_change] <- lapply(
  nnd_combined[cols_to_change],
  as.Date
)

period_levels <- unique(point_level_nnd$period[order(point_level_nnd$date_prev)])

## Colour palette
n_periods <- length(period_levels)
colour_vals <- c(
  colorRampPalette(c("#3d1f0f", "#c06d4c", "#f5d5c5"))(n_periods),
  colorRampPalette(c("#0d2a33", "#41849c", "#c2e0ea"))(n_periods),
  colorRampPalette(c("#0f3533", "#5cb4b0", "#c5ecea"))(n_periods)
) |> setNames(c(
  paste("Velavadar_LEK1", period_levels, sep = "."),
  paste("Velavadar_LEK2", period_levels, sep = "."),
  paste("TalChhapar_TC",  period_levels, sep = ".")
))


#### DENSITY PLOT: TWO SOURCES ####
nnd_combined |>
  mutate(source = factor(source, levels = c("observed", "CSR"))) |>
  ggplot(aes(x = nnd_to_prev, colour = source, fill = source)) +
  geom_density(adjust = 2, linewidth = 0.8, alpha = 0.25) +
  geom_vline(xintercept = 4.7, linewidth = 0.5, linetype = "dashed", colour = "grey30") +
  geom_text(
    data = data.frame(source = factor("observed", levels = c("observed", "CSR"))),
    aes(x = 4.7, y = Inf, label = "4.7 m"),
    hjust = -0.15, vjust = 1.5, size = 3, colour = "grey30", fontface = "italic",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ source, ncol = 1, strip.position = "right") +
  coord_cartesian(xlim = c(1, 50), ylim = c(0,0.08)) +
  scale_colour_manual(values = c("observed" = "#2c3e6b", "CSR" = "#c06d4c"), guide = "none") +
  scale_fill_manual(values  = c("observed" = "#2c3e6b", "CSR" = "#c06d4c"), guide = "none") +
  scale_x_continuous(breaks = c(1, 10, 25, 50, 75, 100), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Nearest neighbour distance to previous survey (m)",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text        = element_text(size = 10, face = "bold"),
    strip.background  = element_blank(),
    axis.line         = element_line(colour = "black", linewidth = 0.4),
    axis.ticks        = element_line(colour = "black", linewidth = 0.4),
    axis.text         = element_text(colour = "black", size = 10),
    axis.title        = element_text(colour = "black", size = 11),
    panel.spacing     = unit(0.8, "lines"),
    plot.margin       = margin(10, 15, 10, 10)
  )




### Three metrics ####
#cross year nnd
#nearest number of neighbouts 
#distance to lek centre

#### PAIRWISE PLOTS ####

add_voronoi_neighbours <- function(df) {
  df |>
    group_by(lek_id, period, source) |>
    group_modify(~{
      d <- .x
      if (nrow(d) < 3) { d$n_neighbours <- NA_integer_; return(d) }
      dd    <- deldir(d$x, d$y)
      cnt   <- table(c(dd$delsgs$ind1, dd$delsgs$ind2))
      n_nbr <- as.integer(cnt[as.character(seq_len(nrow(d)))])
      n_nbr[is.na(n_nbr)] <- 0L
      d$n_neighbours <- n_nbr
      d
    }) |>
    ungroup()
}

## Source
source <- "CSR"
#source <- "observed"

tal_nbr <- nnd_combined |>
  filter(lek_id == "TalChhapar_TC", source == source) |>
  add_voronoi_neighbours()

tal_theme <- theme_classic(base_size = 11) +
  theme(axis.line   = element_line(colour = "black", linewidth = 0.35),
        axis.ticks  = element_line(colour = "black", linewidth = 0.35),
        axis.text   = element_text(colour = "black", size = 9),
        axis.title  = element_text(colour = "black", size = 10),
        plot.tag    = element_text(face = "bold", size = 11),
        plot.margin = margin(6, 10, 6, 10))

## cross-year NND ~ distance to lek centre
p_tal1 <- tal_nbr |>
  ggplot(aes(x = dist_to_centre, y = nnd_to_prev)) +
  geom_point(alpha = 0.25, size = 0.8, colour = "#2c3e6b") +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,
              linewidth = 0.85, colour = "grey15", fill = "grey75", alpha = 0.35) +
  labs(tag = "A", x = "Distance to lek centre (m)", y = "Cross-year NND (m)") +
  tal_theme

## cross-year NND ~ number of Voronoi neighbours
p_tal2 <- tal_nbr |>
  ggplot(aes(x = n_neighbours, y = nnd_to_prev)) +
  geom_jitter(alpha = 0.25, size = 0.8, width = 0.2, colour = "#2c3e6b") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE,
              linewidth = 0.85, colour = "grey15", fill = "grey75", alpha = 0.35) +
  scale_x_continuous(breaks = 1:12) +
  labs(tag = "B", x = "Number of Voronoi neighbours", y = "Cross-year NND (m)") +
  tal_theme

## distance to lek centre ~ number of Voronoi neighbours
p_tal3 <- tal_nbr |>
  ggplot(aes(x = n_neighbours, y = dist_to_centre)) +
  geom_jitter(alpha = 0.25, size = 0.8, width = 0.2, colour = "#2c3e6b") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE,
              linewidth = 0.85, colour = "grey15", fill = "grey75", alpha = 0.35) +
  scale_x_continuous(breaks = 1:12) +
  labs(tag = "C", x = "Number of Voronoi neighbours", y = "Distance to lek centre (m)") +
  tal_theme

fig_tal_pairwise <- (p_tal1 | p_tal2 | p_tal3) +
  plot_annotation(
    title = paste0("Tal Chhapar (", source, ")"),
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  )

print(fig_tal_pairwise)
ggsave(file.path(out_dir, paste0("talchhapar_pairwise_metrics_", source, ".pdf")),
       fig_tal_pairwise, width = 13, height = 5, dpi = 300)



#### PLOTTING PREV-YEAR vs CURRENT-YEAR KDE PEAKS ####

## Build peak locations for every consecutive survey pair
kde_peaks_diag <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  sub_tbl     <- files_tbl |> filter(lek_id == lk) |> arrange(date)
  lek_polygon <- lek_polygons[[lk]]
  sigma       <- sigma_fixed_tbl$sigma_fixed[sigma_fixed_tbl$lek_id == lk]

  map_dfr(2:nrow(sub_tbl), function(i) {
    row_prev <- sub_tbl[i - 1, ]
    row_now  <- sub_tbl[i, ]

    pts_prev <- read.csv(row_prev$csv_path) |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) |>
      clip_points_to_polygon(lek_polygon)
    pts_now  <- read.csv(row_now$csv_path) |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) |>
      clip_points_to_polygon(lek_polygon)

    if (nrow(pts_prev) < min_points_kde || nrow(pts_now) < min_points_kde) return(NULL)

    kde_prev  <- make_kde_grid(pts_prev, lek_polygon, sigma, core_prob, kde_dimyx)
    kde_now   <- make_kde_grid(pts_now,  lek_polygon, sigma, core_prob, kde_dimyx)
    pk_prev   <- kde_prev$kde_df |> filter(is.finite(p)) |> slice_max(p, n = 1, with_ties = FALSE)
    pk_now    <- kde_now$kde_df  |> filter(is.finite(p)) |> slice_max(p, n = 1, with_ties = FALSE)
    prd       <- paste(format(row_prev$date, "%b%y"), format(row_now$date, "%b%y"), sep = "→")

    tibble(
      cx          = c(pk_prev$x,        pk_now$x),
      cy          = c(pk_prev$y,        pk_now$y),
      centre_type = c("Prev-year peak", "Current-year peak"),
      lek_id      = lk,
      period      = prd
    )
  })
})

## One page per lek inside a single PDF
pdf(file.path(out_dir, "kde_peaks_diagnostic.pdf"), width = 14, height = 10)

for (lk in unique(files_tbl$lek_id)) {

  period_order <- unique(
    nnd_combined$period[nnd_combined$lek_id == lk][
      order(nnd_combined$date_prev[nnd_combined$lek_id == lk])
    ]
  )

  pts_lk <- nnd_combined |>
    filter(lek_id == lk, source == "observed") |>
    mutate(
      in_core_lab = if_else(in_prev_core, "In prev-year core", "Outside core"),
      period      = factor(period, levels = period_order)
    )

  peaks_lk <- kde_peaks_diag |>
    filter(lek_id == lk) |>
    mutate(
      centre_type = factor(centre_type, levels = c("Prev-year peak", "Current-year peak")),
      period      = factor(period, levels = period_order)
    )

  lek_geom <- st_geometry(lek_polygons[[lk]])

  p <- ggplot() +
    geom_sf(data = lek_geom, fill = "grey96", colour = "black", linewidth = 0.5) +
    geom_point(
      data  = pts_lk,
      aes(x = x, y = y, colour = in_core_lab),
      size  = 1, alpha = 0.55, shape = 16
    ) +
    geom_point(
      data   = peaks_lk,
      aes(x  = cx, y = cy, fill = centre_type, shape = centre_type),
      size   = 4.5, stroke = 0.8, colour = "black"
    ) +
    facet_wrap(~ period, ncol = 5) +
    scale_colour_manual(
      values = c("In prev-year core" = "#2c3e6b", "Outside core" = "grey70"),
      name   = "Territory point"
    ) +
    scale_fill_manual(
      values = c("Prev-year peak" = "#c06d4c", "Current-year peak" = "#2ecc71"),
      name   = "KDE peak"
    ) +
    scale_shape_manual(
      values = c("Prev-year peak" = 23, "Current-year peak" = 24),
      name   = "KDE peak"
    ) +
    coord_sf(datum = st_crs(32643)) +
    labs(
      title    = paste("KDE peaks diagnostic —", lk),
      subtitle = paste0(
        "Orange diamond = prev-year KDE peak  |  ",
        "Green triangle = current-year KDE peak  |  ",
        "Blue dots = current-year territories in prev-year core"
      ),
      x = NULL, y = NULL
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 3)),
      fill   = guide_legend(override.aes = list(size = 4)),
      shape  = guide_legend(override.aes = list(size = 4))
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position  = "bottom",
      legend.key.size  = unit(0.5, "cm"),
      legend.text      = element_text(size = 8),
      axis.text        = element_text(size = 5),
      strip.text       = element_text(size = 7, face = "bold"),
      plot.title       = element_text(size = 12, face = "bold"),
      plot.subtitle    = element_text(size = 8, colour = "grey40"),
      panel.grid       = element_blank()
    )

  print(p)
}

dev.off()

  

#### VORONOI NEIGHBOURS ~ DISTANCE TO CENTRE ####

nnd_full <- bind_rows(
  point_level_nnd     |> mutate(source = "observed"),
  point_level_nnd_sim |> mutate(source = "CSR")
) |>
  mutate(across(c(date_prev, date_now), as.Date))

## Voronoi neighbour count on ALL points, per (lek, period, source)
nnd_nbr_full <- add_voronoi_neighbours(nnd_full)

## Keep core points only for the model + plot
nbr_core <- nnd_nbr_full |>
  filter(in_prev_core, !is.na(n_neighbours)) |>
  mutate(source = factor(source, levels = c("observed", "CSR")))

## Linear models: n_neighbours ~ dist_to_centre, per lek and source
mods <- nbr_core |>
  group_by(lek_id, source) |>
  group_modify(~ broom::tidy(lm(n_neighbours ~ dist_to_centre, data = .x))) |>
  ungroup()
print(mods, n = Inf)

## Plot
p_nbr_dist <- nbr_core |>
  ggplot(aes(x = dist_to_centre, y = n_neighbours, colour = source, fill = source)) +
  geom_hline(yintercept = 6, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
  geom_jitter(height = 0.2, alpha = 0.18, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.25) +
  facet_wrap(~ lek_id, scales = "free_x") +
  scale_colour_manual(values = c("observed" = "#2c3e6b", "CSR" = "#c06d4c")) +
  scale_fill_manual(values  = c("observed" = "#2c3e6b", "CSR" = "#c06d4c")) +
  labs(x = "Distance to lek centre (m)",
       y = "Number of Voronoi neighbours",
       colour = "Source", fill = "Source") +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black", linewidth = 0.4),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(colour = "black", size = 11),
    legend.position  = "bottom",
    panel.spacing    = unit(0.8, "lines")
  )

print(p_nbr_dist)
ggsave(file.path(out_dir, "voronoi_neighbours_vs_centre.pdf"),
       p_nbr_dist, width = 12, height = 4.5, dpi = 300)


#### LOGNORMAL MODEL: CROSS-YEAR NND ~ DISTANCE TO LEK CENTRE ####

nnd_mod <- nnd_combined |>
  filter(source == "observed", nnd_to_prev > 0, is.finite(dist_to_centre))

m_lnorm <- lm(log(nnd_to_prev) ~ dist_to_centre, data = nnd_mod)
summary(m_lnorm)

## Back-transformed slope: multiplicative change in NND per +1 m from centre
exp(coef(m_lnorm)["dist_to_centre"])

## Plot on log y-axis (data + lognormal fit)
p_lnorm <- nnd_mod |>
  ggplot(aes(x = dist_to_centre, y = nnd_to_prev)) +
  geom_point(alpha = 0.18, size = 0.7, colour = "#2c3e6b") +
  geom_smooth(method = "lm", se = TRUE,
              linewidth = 0.9, colour = "grey15", fill = "grey75", alpha = 0.35) +
  scale_y_log10() +
  labs(x = "Distance to lek centre (m)",
       y = "Cross-year NND (m, log scale)") +
  theme_classic(base_size = 12) +
  theme(
    axis.line  = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    axis.text  = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 11)
  )

print(p_lnorm)
ggsave(file.path(out_dir, "nnd_lognormal_vs_centre.pdf"),
       p_lnorm, width = 6, height = 4.5, dpi = 300)


#### VORONOI NEIGHBOUR-COUNT DISTRIBUTION: OBSERVED vs CSR, BY LEK ####

if (!exists("nbr_core")) {
  nbr_core <- bind_rows(
    point_level_nnd     |> mutate(source = "observed"),
    point_level_nnd_sim |> mutate(source = "CSR")
  ) |>
    mutate(across(c(date_prev, date_now), as.Date)) |>
    add_voronoi_neighbours() |>
    filter(in_prev_core, !is.na(n_neighbours)) |>
    mutate(source = factor(source, levels = c("observed", "CSR")))
}

## Proportion of core points at each neighbour count, per lek and source
nbr_dist <- nbr_core |>
  count(lek_id, source, n_neighbours) |>
  group_by(lek_id, source) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

## Mean neighbour count per lek and source
nbr_mean <- nbr_core |>
  group_by(lek_id, source) |>
  summarise(mean_nbr = mean(n_neighbours), .groups = "drop")

p_nbr_dist <- nbr_dist |>
  ggplot(aes(x = n_neighbours, y = prop, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, alpha = 0.85) +
  geom_vline(data = nbr_mean,
             aes(xintercept = mean_nbr, colour = source),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = 6, linetype = "dotted", linewidth = 0.4, colour = "grey50") +
  facet_wrap(~ lek_id, ncol = 1) +
  scale_x_continuous(breaks = 1:14) +
  scale_fill_manual(values   = c("observed" = "#2c3e6b", "CSR" = "#c06d4c")) +
  scale_colour_manual(values = c("observed" = "#2c3e6b", "CSR" = "#c06d4c")) +
  labs(x = "Number of Voronoi neighbours",
       y = "Proportion of core territories",
       fill = "Source") +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black", linewidth = 0.4),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(colour = "black", size = 11),
    legend.position  = "bottom",
    panel.spacing    = unit(0.8, "lines")
  )

print(p_nbr_dist)
ggsave(file.path(out_dir, "voronoi_neighbour_distribution_by_lek.pdf"),
       p_nbr_dist, width = 7, height = 8, dpi = 300)