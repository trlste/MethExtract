library(data.table)
library(ggplot2)

p_cutoff <- 0.0001
diff_cutoff <- 30

# result one
tile_dt_1 <- fread("diff_methyl_0518/dm_results_cd4-cd8_0518_2_tile_level.csv")

random_dt <- copy(tile_dt_1)

random_summary_dt <- random_dt[
  ,
  .(
    n_total = .N,
    n_beyond_cutoff = sum(pvalue < p_cutoff & abs(methylation_difference) > diff_cutoff, na.rm = TRUE),
    fraction_beyond_cutoff = mean(pvalue < p_cutoff & abs(methylation_difference) > diff_cutoff, na.rm = TRUE)
  ),
  by = segment_count
]

random_avg_width <- random_dt[
  ,
  .(avg_width = mean(width, na.rm = TRUE)),
  by = segment_count
]

width_fraction_random_dt <- merge(
  random_avg_width,
  random_summary_dt,
  by = "segment_count"
)


# result two
tile_dt_2_1 <- fread("coseg/dm_fraction_sig_results_cd4-cd8_0528_tile_level.csv")
tile_dt_2_2 <- fread("coseg/dm_fraction_sig_results_cd4-cd8_0607_tile_level.csv")
tile_dt_2 <- rbindlist(list(tile_dt_2_1, tile_dt_2_2))

coseg_dt <- copy(tile_dt_2)
coseg_summary_dt <- coseg_dt[
  ,
  .(
    n_total = .N,
    n_beyond_cutoff = sum(pvalue < p_cutoff & abs(methylation_difference) > diff_cutoff, na.rm = TRUE),
    fraction_beyond_cutoff = mean(pvalue < p_cutoff & abs(methylation_difference) > diff_cutoff, na.rm = TRUE)
  ),
  by = coseg_file
]

coseg_avg_width <- coseg_dt[
  ,
  .(avg_width = mean(width, na.rm = TRUE)),
  by = coseg_file
]

width_fraction_coseg_dt <- merge(
  coseg_avg_width,
  coseg_summary_dt,
  by = "coseg_file"
)

# result 3
tile_dt_3 <- fread("coseg/dm_results_coseg_cd4-cd8_0629_tile_level.csv")

atlas_dt <- copy(tile_dt_3)

atlas_summary_dt <- atlas_dt[
  ,
  .(
    n_total = .N,
    n_beyond_cutoff = sum(pvalue < p_cutoff & abs(methylation_difference) > diff_cutoff, na.rm = TRUE),
    fraction_beyond_cutoff = mean(pvalue < p_cutoff & abs(methylation_difference) > diff_cutoff, na.rm = TRUE)
  )
]

atlas_avg_width <- atlas_dt[
  ,
  .(avg_width = mean(width, na.rm = TRUE))
]

width_fraction_atlas_dt <- cbind(
  atlas_avg_width,
  atlas_summary_dt
)


# combine results
random_plot_dt <- copy(width_fraction_random_dt)

random_plot_dt[, `:=`(
  method = "Random tiles",
  label = as.character(segment_count)
)]

random_plot_dt <- random_plot_dt[
  ,
  .(
    method,
    label,
    avg_width,
    n_total,
    n_beyond_cutoff,
    fraction_beyond_cutoff
  )
]

coseg_plot_dt <- copy(width_fraction_coseg_dt)

label_map <- c(
  "output" = "10%",
  "output_9_per" = "9%",
  "output_8_per" = "8%",
  "output_7_per" = "7%",
  "output_6_per" = "6%",
  "output_5_per" = "5%",
  "output_2.5_per" = "2.5%",
  "output_1.25_per" = "1.25%",
  "output_0.625_per" = "0.625%"
)

coseg_plot_dt[, `:=`(
  method = "Co-segmentation",
  label = label_map[coseg_file]
)]

coseg_plot_dt <- coseg_plot_dt[
  ,
  .(
    method,
    label,
    avg_width,
    n_total,
    n_beyond_cutoff,
    fraction_beyond_cutoff
  )
]

atlas_plot_dt <- copy(width_fraction_atlas_dt)

atlas_plot_dt[, `:=`(
  method = "Atlas blocks",
  label = "Atlas"
)]

atlas_plot_dt <- atlas_plot_dt[
  ,
  .(
    method,
    label,
    avg_width,
    n_total,
    n_beyond_cutoff,
    fraction_beyond_cutoff
  )
]

plot_dt <- rbindlist(
  list(atlas_plot_dt, random_plot_dt, coseg_plot_dt),
  use.names = TRUE
)

p <- ggplot(
  plot_dt,
  aes(
    x = avg_width,
    y = fraction_beyond_cutoff,
    color = method,
    shape = method,
    group = method
  )
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_text(
    aes(label = label),
    vjust = -0.6,
    size = 4,
    show.legend = FALSE
  ) +
  labs(
    x = "Average width",
    y = "Fraction",
    color = "Method",
    shape = "Method",
    title = "Fraction of p < 1e-4 & |meth diff| > 0.3"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

p <- p +
  scale_x_log10()

p

ggsave(
  "avg_width_vs_fraction_random_vs_coseg_0607_atlas.pdf",
  p,
  width = 8,
  height = 6
)

# Remove invalid/missing widths
atlas_width_dt <- atlas_dt[!is.na(width) & width > 0]

p1 <- ggplot(atlas_width_dt, aes(x = width)) +
  geom_histogram(bins = 100, color = "black", fill = "steelblue") +
  scale_x_log10() +
  labs(
    title = "Distribution of Atlas Blocks Widths",
    x = "Segment width (bp, log10 scale)",
    y = "Number of segments"
  ) +
  theme_classic()

p1

ggsave(
  "distribution_width_atlas.pdf",
  p1,
  width = 8,
  height = 6
)