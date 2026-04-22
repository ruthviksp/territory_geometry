rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(tidyverse)
library(rstudioapi)

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
wd <- getwd()
root_dir <- "C:/Users/ruthv/OneDrive/Desktop/Uni/post_MEME/blackbuck_leks/blackbuck_data_wrangling/territory_geometry/rawdata"
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
out_dir <- file.path(script_dir, "processed_data_rp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Master table
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

compute_point_to_prev_nnd <- function(pts_now_sf, pts_prev_sf, lek_polygon) {
  W <- as.owin(st_geometry(lek_polygon))
  
  xy_now <- st_coordinates(pts_now_sf)
  xy_prev <- st_coordinates(pts_prev_sf)
  
  # Only inside lek
  keep_now <- inside.owin(xy_now[,1], xy_now[,2], W)
  keep_prev <- inside.owin(xy_prev[,1], xy_prev[,2], W)
  
  xy_now <- xy_now[keep_now, , drop = FALSE]
  xy_prev <- xy_prev[keep_prev, , drop = FALSE]
  
  if(nrow(xy_prev) == 0 || nrow(xy_now) == 0) {
    return(tibble(nnd_to_prev = numeric(0)))
  }
  
  X_now <- ppp(xy_now[,1], xy_now[,2], window = W)
  X_prev <- ppp(xy_prev[,1], xy_prev[,2], window = W)
  
  dists <- nncross(X_now, X_prev)$dist
  dists <- as.numeric(dists)
  
  tibble(nnd_to_prev = dists)
}


point_level_nnd <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  sub_tbl <- files_tbl %>% filter(lek_id == lk) %>% arrange(date)
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    row_prev <- sub_tbl[i-1,]
    row_now <- sub_tbl[i,]
    
    lek_polygon <- st_read(file.path(root_dir, row_now$location, row_now$shp_file), 
                           quiet = TRUE) |> 
      st_transform(32643) |> st_zm(drop = TRUE)
    
    pts_prev <- read.csv(row_prev$csv_path) |> 
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    pts_now <- read.csv(row_now$csv_path) |> 
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    
    nnd_table <- compute_point_to_prev_nnd(pts_now, pts_prev, lek_polygon)
    
    nnd_table %>%
      mutate(
        lek_id = lk,
        date_prev = row_prev$date,
        date_now = row_now$date,
        period = paste(year(row_prev$date), year(row_now$date), sep = "-")
      )
  })
})

# Verify
nrow(point_level_nnd) 
table(point_level_nnd$lek_id) 

# Simulate CSR for each timestep pair, matched n and area

point_level_nnd_sim <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  sub_tbl <- files_tbl |> filter(lek_id == lk) |> arrange(date)
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    row_prev <- sub_tbl[i-1,]
    row_now  <- sub_tbl[i,]
    
    lek_polygon <- st_read(file.path(root_dir, row_now$location, row_now$shp_file),
                           quiet = TRUE) |>
      st_transform(32643) |> st_zm(drop = TRUE)
    W_broad <- as.owin(st_geometry(lek_polygon))
    
    # Read and clip to shapefile first
    xy_prev <- read.csv(row_prev$csv_path) |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) |>
      st_coordinates()
    xy_now <- read.csv(row_now$csv_path) |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) |>
      st_coordinates()
    
    xy_prev <- xy_prev[inside.owin(xy_prev[,1], xy_prev[,2], W_broad), , drop = FALSE]
    xy_now  <- xy_now[inside.owin(xy_now[,1],  xy_now[,2],  W_broad), , drop = FALSE]
    
    if (nrow(xy_prev) < 3 || nrow(xy_now) < 3) return(NULL)
    
    # MCP of each timestep's clipped points as owin
    W_prev <- xy_prev |>
      as.data.frame() |>
      st_as_sf(coords = c("X", "Y"), crs = 32643) |>
      st_union() |>
      st_convex_hull() |>
      as.owin()
    
    W_now <- xy_now |>
      as.data.frame() |>
      st_as_sf(coords = c("X", "Y"), crs = 32643) |>
      st_union() |>
      st_convex_hull() |>
      as.owin()
    
    # Simulate matched n inside each timestep's own MCP window
    sim_prev <- runifpoint(nrow(xy_prev), win = W_prev)
    sim_now  <- runifpoint(nrow(xy_now),  win = W_now)
    
    # Unrestricted nncross, like-for-like with observed
    dists <- nncross(sim_now, sim_prev)$dist
    
    tibble(nnd_to_prev = as.numeric(dists),
           lek_id    = lk,
           date_prev = row_prev$date,
           date_now  = row_now$date,
           period    = paste(year(row_prev$date), year(row_now$date), sep = "-"))
  })
})


point_level_nnd_sim <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  sub_tbl <- files_tbl |> filter(lek_id == lk) |> arrange(date)
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    row_prev <- sub_tbl[i-1,]
    row_now  <- sub_tbl[i,]
    
    lek_polygon <- st_read(file.path(root_dir, row_now$location, row_now$shp_file),
                           quiet = TRUE) |>
      st_transform(32643) |> st_zm(drop = TRUE)
    W <- as.owin(st_geometry(lek_polygon))
    
    # Get observed n for each timestep
    obs_prev <- read.csv(row_prev$csv_path) |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    obs_now  <- read.csv(row_now$csv_path)  |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    
    xy_prev <- st_coordinates(obs_prev)
    xy_now  <- st_coordinates(obs_now)
    n_prev  <- sum(inside.owin(xy_prev[,1], xy_prev[,2], W))
    n_now   <- sum(inside.owin(xy_now[,1],  xy_now[,2],  W))
    
    if (n_prev == 0 || n_now == 0) return(NULL)
    
    # Simulate random points inside W
    sim_prev <- runifpoint(n_prev, win = W)
    sim_now  <- runifpoint(n_now,  win = W)
    
    dists <- nncross(sim_now, sim_prev)$dist
    
    tibble(nnd_to_prev = as.numeric(dists),
           lek_id    = lk,
           date_prev = row_prev$date,
           date_now  = row_now$date,
           period    = paste(year(row_prev$date), year(row_now$date), sep = "-"))
  })
})

point_level_nnd_grid <- map_dfr(unique(files_tbl$lek_id), function(lk) {
  sub_tbl <- files_tbl |> filter(lek_id == lk) |> arrange(date)
  
  map_dfr(2:nrow(sub_tbl), function(i) {
    row_prev <- sub_tbl[i-1,]
    row_now  <- sub_tbl[i,]
    
    lek_polygon <- st_read(file.path(root_dir, row_now$location, row_now$shp_file),
                           quiet = TRUE) |>
      st_transform(32643) |> st_zm(drop = TRUE)
    W_broad <- as.owin(st_geometry(lek_polygon))
    
    xy_prev <- read.csv(row_prev$csv_path) |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) |>
      st_coordinates()
    xy_now <- read.csv(row_now$csv_path) |>
      st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643) |>
      st_coordinates()
    
    xy_prev <- xy_prev[inside.owin(xy_prev[,1], xy_prev[,2], W_broad), , drop = FALSE]
    xy_now  <- xy_now[inside.owin(xy_now[,1],  xy_now[,2],  W_broad), , drop = FALSE]
    
    if (nrow(xy_prev) < 3 || nrow(xy_now) < 3) return(NULL)
    
    W_prev <- xy_prev |>
      as.data.frame() |>
      st_as_sf(coords = c("X", "Y"), crs = 32643) |>
      st_union() |>
      st_convex_hull() |>
      as.owin()
    
    W_now <- xy_now |>
      as.data.frame() |>
      st_as_sf(coords = c("X", "Y"), crs = 32643) |>
      st_union() |>
      st_convex_hull() |>
      as.owin()
    
    # aspect <- diff(W_now$yrange) / diff(W_now$xrange)
    # rsyst(nx = round(sqrt(n / aspect)), ny = round(sqrt(n * aspect)), win = W_now)
    
    grid_prev <- rsyst(nx = round(sqrt(nrow(xy_prev))), win = W_prev)
    grid_now  <- rsyst(nx = round(sqrt(nrow(xy_now))),  win = W_now)
    
    dists <- nncross(grid_now, grid_prev)$dist
    
    tibble(nnd_to_prev = as.numeric(dists),
           lek_id    = lk,
           date_prev = row_prev$date,
           date_now  = row_now$date,
           period    = paste(year(row_prev$date), year(row_now$date), sep = "-"))
  })
})

# Combine all three
nnd_combined <- bind_rows(
  point_level_nnd      |> mutate(source = "observed"),
  point_level_nnd_sim  |> mutate(source = "CSR"),
  point_level_nnd_grid |> mutate(source = "grid")
)

# Combine observed and simulated with a source label
nnd_combined <- bind_rows(
  point_level_nnd |> mutate(source = "observed"),
  point_level_nnd_sim |> mutate(source = "CSR")
)

period_levels <- unique(point_level_nnd$period[order(point_level_nnd$date_prev)])

#colour palatte
colour_vals <- c(
  colorRampPalette(c("#3d1f0f", "#c06d4c", "#f5d5c5"))(17),
  colorRampPalette(c("#0d2a33", "#41849c", "#c2e0ea"))(17),
  colorRampPalette(c("#0f3533", "#5cb4b0", "#c5ecea"))(17)
) |> setNames(c(
  paste("Velavadar_LEK1", period_levels, sep = "."),
  paste("Velavadar_LEK2", period_levels, sep = "."),
  paste("TalChhapar_TC",  period_levels, sep = ".")
))

median_points <- nnd_combined |>
  mutate(period     = factor(period, levels = period_levels),
         colour_key = interaction(lek_id, period)) |>
  group_by(lek_id, period, source, colour_key) |>
  summarise(median_nnd = median(nnd_to_prev), .groups = "drop")

# Cross year nnd distribution and median plot for csr vs grid vs observed
nnd_combined |>
  mutate(period     = factor(period, levels = period_levels),
         colour_key = interaction(lek_id, period)) |>
  ggplot(aes(x = nnd_to_prev, colour = colour_key, linetype = source, linewidth = source)) +
  facet_grid(source ~ lek_id) +
  geom_density(adjust = 2) +
  geom_point(data = median_points,
             aes(x = median_nnd, y = as.numeric(factor(period)) * 0.005,
                 fill = colour_key, shape = source),
             colour = "black", stroke = 0.5, size = 2, inherit.aes = FALSE) +
  coord_cartesian(xlim = c(1, 100)) +
  scale_colour_manual(values = colour_vals) +
  scale_linetype_manual(values = c("observed" = "solid", "CSR" = "dashed", "grid" = "dotted")) +
  scale_linewidth_manual(values = c("observed" = 1, "CSR" = 0.5, "grid" = 0.5)) +
  scale_shape_manual(values = c("observed" = 21, "CSR" = 24, "grid" = 22)) +
  scale_fill_manual(values = colour_vals) +
  guides(colour    = guide_legend(ncol = 1),
         fill      = "none",
         linewidth = "none") +
  labs(x = "NND to previous time step", y = "Density",
       linetype = "Source", shape = "Source") +
  theme_bw(base_size = 11) +
  theme(legend.position   = "right",
        legend.key.size   = unit(0.4, "cm"),
        legend.text       = element_text(size = 8),
        legend.title      = element_blank(),
        strip.text.y      = element_blank())


# #calculate median old
# median_points <- nnd_combined |>
#   mutate(period     = factor(period, levels = period_levels),
#          colour_key = interaction(lek_id, as.integer(period))) |>
#   group_by(lek_id, period, source, colour_key) |>
#   summarise(median_nnd = median(nnd_to_prev), .groups = "drop")
# 
# # combined distribution plot onlt for csr vs observed
# nnd_combined |>
#   mutate(period     = factor(period, levels = period_levels),
#          colour_key = interaction(lek_id, as.integer(period))) |>
#   ggplot(aes(x = nnd_to_prev, colour = colour_key, linetype = source, linewidth = source)) +
#   facet_wrap(source ~ lek_id) +
#   geom_density(adjust = 2) +
#   geom_point(data = median_points,
#              aes(x = median_nnd, y = as.numeric(factor(period)) * 0.005,
#                  fill = colour_key, shape = source),
#              colour = "black", stroke = 0.5, size = 2, inherit.aes = FALSE) +
#   coord_cartesian(xlim = c(1, 100)) +
#   scale_colour_manual(values = colour_vals) +
#   scale_linetype_manual(values = c("observed" = "solid", "CSR" = "dashed", "grid" = "solid")) +
#   scale_linewidth_manual(values = c("observed" = 1, "CSR" = 0.5)) +
#   scale_shape_manual(values = c("observed" = 21, "CSR" = 24)) +
#   scale_fill_manual(values = colour_vals) +
#   guides(colour = "none", linewidth = "none") +
#   labs(x = "NND to previous time step", y = "Density",
#        linetype = "Source", shape = "Source") +
#   theme_bw(base_size = 11) +
#   theme(legend.position   = "right",
#         legend.key.size   = unit(0.4, "cm"),
#         legend.text       = element_text(size = 8),
#         legend.title = element_text(""))


# Median shift plot for cross year nnd
# Median observed - median csr for the same time period 
nnd_combined |>
  mutate(period = factor(period, levels = period_levels)) |>
  group_by(lek_id, period, source) |>
  summarise(median_nnd = median(nnd_to_prev), .groups = "drop") |>
  pivot_wider(names_from = source, values_from = median_nnd) |>
  mutate(shift = observed - CSR) |>
  ggplot(aes(x = period, y = shift, colour = lek_id, group = lek_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 0.8) +
  scale_colour_manual(values = c("#c06d4c", "#41849c", "#5cb4b0")) +
  labs(x = "Period", y = "Cross year NND median observed − median CSR (m)",
       colour = "Lek") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")


# ECDF plot - Emperical Cumulative Distrivution FUnction
# With this, you basically see that for a given value of x (nnd) what proportion of values fall below that value, eventually reaching to 1 within increasing threshold of nnd. 
nnd_combined |>
  mutate(period = factor(period, levels = period_levels)) |>
  ggplot(aes(x = nnd_to_prev, colour = source, linetype = source)) +
  facet_grid(lek_id ~ period) +
  stat_ecdf(linewidth = 0.7) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_colour_manual(values = c("observed" = "#2c3e6b", "CSR" = "#c06d4c")) +
  scale_linetype_manual(values = c("observed" = "solid", "CSR" = "dashed")) +
  labs(x = "NND to previous time step (m)", y = "Cumulative proportion",
       colour = "Source", linetype = "Source") +
  theme_bw(base_size = 9) +
  theme(strip.text.x    = element_text(size = 7),
        strip.text.y    = element_text(size = 8),
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")


#### PLOT TRIALS ####

# Prepare for plotting
point_level_nnd_plot <- point_level_nnd %>%
  mutate(
    period_label = paste(year(date_prev), "\n→", year(date_now)),
    lek_id = factor(lek_id, levels = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"))
  ) %>%
  filter(!is.na(nnd_to_prev) & nnd_to_prev > 0)  # Clean data

# 1. TIME SERIES BOXPLOT (your main request)
p1 <- ggplot(point_level_nnd_plot, aes(x = reorder(period_label, date_prev), 
                                       y = nnd_to_prev, fill = lek_id)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.5, outlier.alpha = 0.6,
               position = position_dodge2(width = 0.9, padding = 0.1)) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Lek Site") +
  labs(
    title = "Territory Stability: Distance to Previous Year",
    #subtitle = "Box = IQR, Whiskers = 1.5×IQR, Points = Outliers",
    x = "Time Period (Prev → Now)", 
    y = "Nearest Neighbor Distance (meters)"
  ) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# 2. Complementary violin plot
p2 <- ggplot(point_level_nnd_plot, aes(x = reorder(period_label, date_prev), 
                                       y = nnd_to_prev, fill = lek_id)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.2, outlier.size = 1, fill = "white", 
               position = position_dodge(width = 0.9)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(title = "Distribution Shape", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# 3. Summary stats table
summary_stats <- point_level_nnd_plot %>%
  group_by(lek_id, period_label) %>%
  summarise(
    n_territories = n(),
    median_nnd = round(median(nnd_to_prev), 1),
    prop_stable = paste0(round(mean(nnd_to_prev <= 5)*100, 0), "%"),
    .groups = "drop"
  )

# DISPLAY & SAVE
print(p1)  # Main plot
print(p2)  # Violin complement

# High-res publication PNGs
ggsave("territory_stability_timeseries.png", p1, width = 14, height = 8, dpi = 300)
ggsave("territory_violin_timeseries.png", p2, width = 14, height = 8, dpi = 300)

# Summary CSV
write_csv(summary_stats, file.path(out_dir, "stability_summary_table.csv"))

cat("✅ SAVED:\n",
    "• territory_stability_timeseries.png (BOXPLOT TIME SERIES)\n",
    "• territory_violin_timeseries.png\n",
    "• stability_summary_table.csv\n")

print("Summary stats preview:")
print(summary_stats)



median_trends <- point_level_nnd %>%
  filter(!is.na(nnd_to_prev) & nnd_to_prev > 0) %>%
  group_by(lek_id, date_prev, period = paste(year(date_prev), "→", year(date_now))) %>%
  summarise(
    median_nnd = median(nnd_to_prev),
    n_territories = n(),
    .groups = "drop"
  ) %>%
  mutate(lek_id = factor(lek_id, levels = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC")))

# MEDIAN TIME SERIES (super clean)
p_median <- ggplot(median_trends, aes(date_prev, median_nnd, color = lek_id, group = lek_id)) +
  geom_line(linewidth = 2, alpha = 0.9) +
  geom_point(size = 4) +
  scale_color_brewer(type = "qual", palette = "Set2", name = "Lek") +
  labs(
    title = "Median Territory Stability Over Time",
    subtitle = "Line = median NN distance to previous year",
    x = "Period Start Date", 
    y = "Median Distance (meters)"
  ) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_median)

print("Medians preview:")
print(median_trends)

write.csv(point_level_nnd, "../processed_data_plots_rp/stability_ALL_new.csv")


# distribution of cross year nnds
point_level_nnd |>
  group_by(period, lek_id) |>
  filter(n() > 0) |> 
  ungroup() |>
  mutate(period = droplevels(factor(period)),
         lek_id = droplevels(factor(lek_id))) |>
  ggplot(aes(x = nnd_to_prev, fill = period, colour = period)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  coord_cartesian(xlim = c(1, 10)) +
  facet_wrap(~ period + lek_id,
             labeller = labeller(.multi_line = FALSE)) +
  facet_grid(period ~ lek_id) +
  scale_fill_viridis_d(alpha = 0.7) +
  scale_colour_viridis_d() +
  labs(x = "NND to previous survey", y = "Density") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9),
        legend.position = "none")


# All distribution in one plot
point_level_nnd |>
  mutate(lek_period = paste(lek_id, period, sep = " | "),
         period_order = as.numeric(factor(date_prev))) |>
  arrange(lek_id, period_order) |>
  mutate(lek_period = factor(lek_period, levels = unique(lek_period))) |>
  ggplot(aes(x = nnd_to_prev, fill = lek_period, colour = lek_period)) +
  geom_density(alpha = 0.05, linewidth = 1) +
  coord_cartesian(xlim = c(1, 20)) +
  scale_fill_manual(values = {
    leks <- unique(point_level_nnd$lek_id)
    periods_per_lek <- point_level_nnd |>
      group_by(lek_id) |>
      summarise(n = n_distinct(date_prev)) |>
      deframe()
    lek_colours <- setNames(RColorBrewer::brewer.pal(length(leks), "Set1"), leks)
    unlist(lapply(leks, function(l) {
      n <- periods_per_lek[[l]]
      colorRampPalette(c(lek_colours[[l]], "white"))(n + 2)[1:n]
    }))
  }) +
  scale_colour_manual(values = {
    leks <- unique(point_level_nnd$lek_id)
    periods_per_lek <- point_level_nnd |>
      group_by(lek_id) |>
      summarise(n = n_distinct(date_prev)) |>
      deframe()
    lek_colours <- setNames(RColorBrewer::brewer.pal(length(leks), "Set1"), leks)
    unlist(lapply(leks, function(l) {
      n <- periods_per_lek[[l]]
      colorRampPalette(c(lek_colours[[l]], "white"))(n + 2)[1:n]
    }))
  }) +
  labs(x = "NND to previous time step", y = "Density",
       fill = "Lek | Period", colour = "Lek | Period") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8))


# separate distribution plots by lek_id

point_level_nnd |>
  mutate(period = factor(period, levels = unique(period[order(date_prev)]))) |>
  ggplot(aes(x = nnd_to_prev, linetype = "solid",
             colour = interaction(lek_id, as.integer(period)))) +
  facet_wrap(~ lek_id) +
  geom_density(linewidth = 1, adjust = 2) +
  coord_cartesian(xlim = c(1, 100)) +
  scale_colour_manual(values = c(
    colorRampPalette(c("#3d1f0f", "#c06d4c", "#f5d5c5"))(17),
    colorRampPalette(c("#0d2a33", "#41849c", "#c2e0ea"))(17),
    colorRampPalette(c("#0f3533", "#5cb4b0", "#c5ecea"))(17)
  ) |> setNames(c(
    paste("Velavadar_LEK1", 1:17, sep = "."),
    paste("Velavadar_LEK2", 1:17, sep = "."),
    paste("TalChhapar_TC",  1:17, sep = ".")
  ))) +
  guides(colour = "none") +
  labs(x = "NND to previous time step", y = "Density", linetype = "Period") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8))

# for loop with each layer being added.
period_levels <- unique(point_level_nnd$period[order(point_level_nnd$date_prev)])

colour_vals <- c(
  colorRampPalette(c("#3d1f0f", "#c06d4c", "#f5d5c5"))(17),
  colorRampPalette(c("#0d2a33", "#41849c", "#c2e0ea"))(17),
  colorRampPalette(c("#0f3533", "#5cb4b0", "#c5ecea"))(17)
) |> setNames(c(
  paste("Velavadar_LEK1", 1:17, sep = "."),
  paste("Velavadar_LEK2", 1:17, sep = "."),
  paste("TalChhapar_TC",  1:17, sep = ".")
))

plot_data <- point_level_nnd |>
  mutate(period = factor(period, levels = period_levels),
         colour_key = interaction(lek_id, as.integer(period)))

for (i in seq_along(period_levels)) {
  p <- plot_data |>
    filter(as.integer(period) <= i) |>
    ggplot(aes(x = nnd_to_prev, colour = colour_key, group = colour_key)) +
    facet_wrap(~ lek_id) +
    geom_density(linewidth = 1, adjust = 2) +
    coord_cartesian(xlim = c(1, 100)) +
    scale_colour_manual(
      values = colour_vals,
      breaks = paste(c("Velavadar_LEK1","Velavadar_LEK2","TalChhapar_TC"),
                     rep(1:i, each = 3), sep = "."),
      labels = rep(period_levels[1:i], each = 3),
      name   = "Period"
    ) +
    labs(x = "NND to previous time step", y = "Density",
         title = paste("Periods shown: 1 to", i, "—", period_levels[i])) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          legend.key.size = unit(0.4, "cm"),
          legend.text     = element_text(size = 8))
  
  print(p)
  if (i < length(period_levels)) readline(prompt = "Press Enter to add next period...")
}


# gif 
library(magick)

period_levels <- unique(point_level_nnd$period[order(point_level_nnd$date_prev)])

colour_vals <- c(
  colorRampPalette(c("#3d1f0f", "#c06d4c", "#f5d5c5"))(17),
  colorRampPalette(c("#0d2a33", "#41849c", "#c2e0ea"))(17),
  colorRampPalette(c("#0f3533", "#5cb4b0", "#c5ecea"))(17)
) |> setNames(c(
  paste("Velavadar_LEK1", 1:17, sep = "."),
  paste("Velavadar_LEK2", 1:17, sep = "."),
  paste("TalChhapar_TC",  1:17, sep = ".")
))

plot_data <- point_level_nnd |>
  mutate(period = factor(period, levels = period_levels),
         colour_key = interaction(lek_id, as.integer(period)))

frames <- image_graph(width = 1400, height = 600, res = 150)

for (i in seq_along(period_levels)) {
  p <- plot_data |>
    filter(as.integer(period) <= i) |>
    ggplot(aes(x = nnd_to_prev, colour = colour_key, group = colour_key)) +
    facet_wrap(~ lek_id) +
    geom_density(linewidth = 1, adjust = 2) +
    coord_cartesian(xlim = c(1, 100)) +
    scale_colour_manual(
      values = colour_vals,
      breaks = paste(c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
                     rep(1:i, each = 3), sep = "."),
      labels = rep(period_levels[1:i], each = 3),
      name   = "Period"
    ) +
    labs(x = "NND to previous time step", y = "Density",
         title = period_levels[i]) +
    theme_bw(base_size = 11) +
    theme(legend.position   = "right",
          legend.key.size   = unit(0.4, "cm"),
          legend.text       = element_text(size = 8))
  print(p)
}

dev.off()

#making a gif

image_animate(frames, fps = 2, optimize = TRUE) |>
  image_write("../processed_data_plots_rp/cross_year_nnd_dist.gif")



