library(tidyverse)
library(zoo)
library(scales)

z_scaffolds <- c(
  "CM042602.1", "QZWM02001482.1", "QZWM02002897.1", "QZWM02003046.1", "QZWM02003073.1",
  "QZWM02004069.1", "QZWM02004073.1", "QZWM02004074.1", "QZWM02004077.1")

analyze_fst <- function(file_path, comparison_name) {
  
  message(paste("Processing:", comparison_name))

  WINDOW_SIZE <- 200000
  STEP_SIZE   <- 50000
  STEPS_PER_WINDOW <- WINDOW_SIZE / STEP_SIZE 

  raw_data <- read_tsv(file_path, show_col_types = FALSE) %>%
    rename(FST = WEIR_AND_COCKERHAM_FST) %>%
    mutate(FST = replace_na(FST, 0)) %>%
    mutate(FST = pmax(FST, 0))
  
  windowed_data <- raw_data %>%
    mutate(step_bin = floor(POS / STEP_SIZE) * STEP_SIZE) %>%
    group_by(CHROM, step_bin) %>%
    summarise(step_fst = mean(FST), .groups = "drop") %>%
    group_by(CHROM) %>%
    arrange(CHROM, step_bin) %>%
    mutate(window_fst = rollmean(step_fst, k = STEPS_PER_WINDOW, fill = NA, align = "center")) %>%
    filter(!is.na(window_fst))

  final_stats <- windowed_data %>%
    ungroup() %>%
    # Label each window as "Z" or "Autosome"
    mutate(chrom_type = ifelse(CHROM %in% z_scaffolds, "Z", "Autosome")) %>%
    group_by(chrom_type) %>%
    mutate(Z_FST = (window_fst - mean(window_fst)) / sd(window_fst)) %>%
    ungroup() %>%
    
    mutate(Comparison = comparison_name)
  
  return(final_stats)
}

pop1_pop2_fst <- analyze_fst("pop1_pop2_fst.weir.fst.gz", "pop1_pop2")
pop1_pop3_fst <- analyze_fst("pop1_pop3_fst.weir.fst.gz", "pop1_pop3")
pop2_pop3_fst <- analyze_fst("pop2_pop3_fst.weir.fst.gz", "pop2_pop3")

fst_results <- bind_rows(pop1_pop2_fst, pop1_pop3_fst, pop2_pop3_fst)

write_csv(fst_results, "fst_results.csv")

# Find significant hits. Z-score threshold = 4.

candidates <- fst_results %>%
  filter(Z_FST > 4) %>%
  rename(Start = step_bin) %>%
  mutate(End = Start + 200000) %>%
  select(CHROM, Start, End, Comparison, Z_FST) %>%
  arrange(desc(Z_FST))

# Merge significant adjacent/overlapping windows

merged_candidates <- candidates %>%
  arrange(CHROM, Start) %>%
  group_by(CHROM) %>%
  mutate(
    prev_end = lag(End, default = 0),
    is_new_region = Start > prev_end, 
    region_id = cumsum(is_new_region)) %>%
  group_by(CHROM, region_id) %>%
  summarise(
    Start = min(Start),
    End = max(End),
    Max_Z_FST = max(Z_FST),
    Top_Comparison = Comparison[which.max(Z_FST)], # Which pair drove this signal?
    n_windows = n(),
    .groups = "drop") %>%
  select(CHROM, Start, End, Max_Z_FST, Top_Comparison, n_windows) %>%
  arrange(desc(Max_Z_FST))

print(head(merged_candidates))

write_csv(merged_candidates, "candidate_regions_z4.csv")

# BED files must be tab-separated with NO header.
# Use on the cluster with bedtools to see which genes these are.
merged_candidates %>%
  select(CHROM, Start, End) %>%
  write_tsv("candidates.bed", col_names = FALSE)