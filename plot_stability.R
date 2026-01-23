rm(list = ls())

## Load libraries
library(tidyverse)

## Save figure parameters
save_pub_fig <- function(plot, filename_base, width = 7, height = 5) {
  ggsave(paste0(filename_base, ".png"),
         plot = plot,
         width = width,
         height = height,
         units = "in",
         dpi = 600)
}

## Load stability metrics
stab <- read.csv('/Users/vivekhsridhar/Library/Mobile Documents/com~apple~CloudDocs/Documents/Code/territory_geometry/processed_data/stability_ALL.csv')

## Factor ordering + labels
stab$lek_id <- factor(stab$lek_id,
                      levels = c("Velavadar_LEK1",
                                 "Velavadar_LEK2",
                                 "TalChhapar_TC"),
                      labels = c("Velavadar Lek 1",
                                 "Velavadar Lek 2",
                                 "Tal Chhapar"))

## Colour palettes
fill_cols <- c(
  "Velavadar Lek 1" = "#4DAF4A",
  "Velavadar Lek 2" = "#377EB8",
  "Tal Chhapar"     = "#D6604D"
)

point_cols <- c(
  "Velavadar Lek 1" = "#1B7837",
  "Velavadar Lek 2" = "#2166AC",
  "Tal Chhapar"     = "#8B1A1A"
)

## ---- Time series: cross-year NN distance ----

p_nn_ts <- ggplot(stab, aes(x = date_now, y = nn_cross_median, colour = lek_id, group = lek_id)) +
  geom_point(size = 2.6) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Year", y = "Cross-year nearest-neighbour distance (m)", colour = "Lek")

save_pub_fig(p_nn_ts, "fig_stability_crossyear_nn_timeseries", width = 8, height = 4.8)

## ---- Time series: mode displacement ----

p_modes_ts <- ggplot(stab, aes(x = date_now, y = mode_shift, colour = lek_id, group = lek_id)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Year", y = "Intensity centroid displacement (m)", colour = "Lek")

save_pub_fig(p_modes_ts, "fig_stability_mode_timeseries", width = 8, height = 4.8)
