rm(list = ls())

## Load libraries
library(tidyverse)

## Save figure parameters
save_pub_fig <- function(plot, filename_base, width = 7, height = 5) {
  ggsave(paste0(filename_base, ".png"), plot = plot,
         width = width, height = height, units = "in", dpi = 600)
}

## Load PCF curves
pcf_all <- read.csv("derived_metrics/pcf_curve_ALL.csv")

## Load peak table
peaks <- read.csv("derived_metrics/pcf_peak_table_ALL.csv")

## Factor ordering + labels (identical to NND plots)
pcf_all$lek_id <- factor(pcf_all$lek_id,
                         levels = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
                         labels = c("Velavadar Lek 1", "Velavadar Lek 2", "Tal Chhapar"))

peaks$lek_id <- factor(peaks$lek_id,
                       levels = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
                       labels = c("Velavadar Lek 1", "Velavadar Lek 2", "Tal Chhapar"))

## Colour palettes
fill_cols <- c(
  "Velavadar Lek 1" = "#4DAF4A",
  "Velavadar Lek 2" = "#377EB8",
  "Tal Chhapar"     = "#D6604D"
)

line_cols <- c(
  "Velavadar Lek 1" = "#1B7837",
  "Velavadar Lek 2" = "#2166AC",
  "Tal Chhapar"     = "#8B1A1A"
)

point_cols <- line_cols

## Mean PCF curve across dates within each lek
pcf_mean <- pcf_all %>%
  group_by(lek_id, r) %>%
  summarise(g_mean = mean(g, na.rm = TRUE), .groups = "drop")

## Plot mean PCF curve with individual curves
p_pcf_mean <- ggplot() +
  geom_line(data = pcf_all, aes(x = r, y = g, group = interaction(lek_id, date), colour = lek_id), linewidth = 0.4, alpha = 0.3) +
  geom_line(data = pcf_mean, aes(x = r, y = g_mean, colour = lek_id), linewidth = 1) +
  scale_colour_manual(values = line_cols) + theme_classic(base_size = 13) +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Distance r (m)", y = expression(g[inhom](r)))

save_pub_fig(p_pcf_mean, "fig_pcf_mean_with_ribbon")

## Plot individual PCF curves with detected peaks overlaid
p_pcf_peaks <- ggplot() +
  geom_line(data = pcf_all, aes(x = r, y = g, group = interaction(lek_id, date), colour = lek_id), linewidth = 0.6, alpha = 0.35) +
  geom_point(data = peaks, aes(x = r_peak, y = g_peak, colour = lek_id), size = 2.2, alpha = 0.9) +
  facet_wrap(~ lek_id, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = line_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Distance r (m)", y = expression(g[inhom](r)))

save_pub_fig(p_pcf_peaks, "fig_pcf_curves_with_detected_peaks", width = 7, height = 9)

## Number of peaks per date
peak_counts <- peaks %>%
  count(lek_id, date, name = "n_peaks")

p_npeaks <- ggplot(peak_counts, aes(x = lek_id, y = n_peaks)) +
  geom_boxplot(aes(fill = lek_id), width = 0.55, outlier.shape = NA, linewidth = 0.8, colour = "black") +
  geom_jitter(aes(colour = lek_id), width = 0.10, height = 0, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = fill_cols) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Number of meso-scale PCF peaks")

save_pub_fig(p_npeaks, "fig_pcf_n_peaks_by_lek")

## Distance of first peak
first_peaks <- peaks %>%
  group_by(lek_id, date) %>%
  summarise(r_first_peak = min(r_peak), .groups = "drop")

p_rpeak <- ggplot(first_peaks, aes(x = lek_id, y = r_first_peak)) +
  geom_boxplot(aes(fill = lek_id), width = 0.55, outlier.shape = NA, linewidth = 0.8, colour = "black") +
  geom_jitter(aes(colour = lek_id), width = 0.10, height = 0, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = fill_cols) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Distance to first meso-scale PCF peak (m)")

save_pub_fig(p_rpeak, "fig_pcf_first_peak_distance_by_lek")

## Peak curvature
p_curvature <- ggplot(peaks, aes(x = lek_id, y = peak_curvature)) +
  geom_boxplot(aes(fill = lek_id), width = 0.55, outlier.shape = NA, linewidth = 0.8, colour = "black") +
  geom_jitter(aes(colour = lek_id), width = 0.10, height = 0, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = point_cols) + theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)), 
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "PCF peak curvature")

save_pub_fig(p_curvature, "fig_pcf_peak_curvature_by_lek")

## Peak width
p_width <- ggplot(peaks, aes(x = lek_id, y = peak_width)) +
  geom_boxplot(aes(fill = lek_id), width = 0.55, outlier.shape = NA, linewidth = 0.8, colour = "black") +
  geom_jitter(aes(colour = lek_id), width = 0.10, height = 0, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = point_cols) + theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)), 
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "PCF peak width")

save_pub_fig(p_width, "fig_pcf_peak_width_by_lek")