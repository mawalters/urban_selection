library(tidyverse)
library(scales)

fst_results <- read_csv("fst_results.csv", show_col_types = FALSE)

# label Z scaffolds
z_scaffolds <- c(
  "CM042602.1", "QZWM02001482.1", "QZWM02002897.1", "QZWM02003046.1", "QZWM02003073.1",
  "QZWM02004069.1", "QZWM02004073.1", "QZWM02004074.1", "QZWM02004077.1")

fst_results <- fst_results %>%
  mutate(chrom_type = ifelse(CHROM %in% z_scaffolds, "Z", "Autosome"))

# put Z last in order
chrom_order <- fst_results %>%
  select(CHROM, chrom_type) %>%
  distinct() %>%
  arrange(chrom_type == "Z", CHROM) %>% 
  pull(CHROM)

chrom_key <- tibble(CHROM = chrom_order) %>%
  mutate(chrom_idx = row_number()) %>%
  mutate(is_odd = chrom_idx %% 2 == 1)

# Create a helper to identify odd/even chromosomes for alternating colors
chrom_key <- tibble(CHROM = chrom_order) %>%
  mutate(chrom_idx = row_number()) %>%
  mutate(is_odd = chrom_idx %% 2 == 1)

# Calculate cumulative positions
data_cum <- fst_results %>%
  mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
  arrange(Comparison, CHROM, step_bin) %>%
  group_by(Comparison, CHROM) %>%
  summarise(max_bp = max(step_bin), .groups = "drop") %>%
  group_by(Comparison) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(Comparison, CHROM, bp_add) %>%
  left_join(fst_results, by = c("Comparison", "CHROM")) %>%
  mutate(bp_cum = step_bin + bp_add) %>%
  left_join(chrom_key, by = "CHROM")

# Calculate 99th percentile threshold per Comparison
thresholds <- data_cum %>%
  group_by(Comparison) %>%
  summarise(threshold_99 = quantile(window_fst, 0.99, na.rm = TRUE), .groups = "drop")

# Assign Color Categories
# Logic: Significant > Z > Autosome (Odd/Even)
data_plot <- data_cum %>%
  left_join(thresholds, by = "Comparison") %>%
  mutate(
    color_group = case_when(
      window_fst >= threshold_99 ~ "Significant",
      chrom_type == "Z" ~ "Z",
      is_odd ~ "Auto_Odd",
      TRUE ~ "Auto_Even"
    )
  ) %>%
  # Sort so 'Significant' points are plotted LAST (on top)
  arrange(color_group != "Significant") 

# Calculate axis labels
scaffold_sizes <- fst_results %>%
  group_by(CHROM) %>%
  summarise(len_mb = max(step_bin) / 1e6) %>%
  arrange(desc(len_mb)) %>% 
  left_join(chrom_key)

top_scaffolds <- scaffold_sizes %>% slice_head(n = 15) %>% pull(chrom_idx)

axis_set <- data_plot %>%
  group_by(Comparison, chrom_idx) %>%
  summarize(center = mean(bp_cum), .groups = "drop") %>%
  mutate(label = ifelse(chrom_idx %in% top_scaffolds, as.character(chrom_idx), ""))

# 3. Generate the Plot
ggplot(data_plot, aes(x = bp_cum, y = window_fst)) +
  geom_point(aes(color = color_group), alpha = 0.9, size = 0.5) +
  geom_hline(data = thresholds, aes(yintercept = threshold_99), 
             linetype = "dashed", color = "firebrick", alpha = 0.5) +
  facet_grid(Comparison ~ .) +
  scale_x_continuous(label = axis_set$label, breaks = axis_set$center) +
  scale_color_manual(
    values = c(
      "Significant" = "firebrick",
      "Z"           = "gold",
      "Auto_Odd"    = "gray30",
      "Auto_Even"   = "gray70"),
    breaks = c("Significant", "Z", "Auto_Odd"), 
    labels = c(
      "Significant" = "Top 1% Outlier",
      "Z"           = "Z Chromosome",
      "Auto_Odd"    = "Autosome")) +
  labs(title = "FST Differentiation",
       x = "Chromosome",
       y = "Mean FST (200kb Window)",
       color = "Legend") +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  )
