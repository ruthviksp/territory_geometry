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

library(spatstat.geom)  # Already loaded

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

# YOUR EXACT LOOP (works perfectly)
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


library(ggplot2)
library(dplyr)

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
