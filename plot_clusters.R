library(tidyverse)

pca <- read_tsv("QC/VCF.eigenvec") %>% 
  select(INDV = `#IID`,
         PC1,
         PC2
         )

inbreeding <- read_tsv("QC/VCF.het") %>% 
  select(INDV,
         F_stat = `F`,
         )

combo_data <- left_join(pca, inbreeding)

# K-means clustering

set.seed(121)

clusters <- kmeans(combo_data %>% select(PC1, PC2), centers = 3)

cluster_data <- combo_data %>%
  mutate(cluster = as.factor(clusters$cluster))

ggplot(cluster_data, aes(x = PC1, y = PC2, fill = cluster)) +
  geom_point(size = 3, pch = 21, alpha = 0.7) +
  theme_bw() +
  labs(title = "PCA with K-Means Clustering")

# Add sites

checksum <- read_delim("MD5.txt", delim = "/", col_names = FALSE)

sites <- checksum %>%
  select(X2) %>% 
  separate(col = X2, into = c("INDV", "site"), sep = "_") %>% 
  distinct(INDV, .keep_all = TRUE) %>% 
  na.omit()

cluster_data <- left_join(cluster_data, sites)

write_csv(cluster_data, "clusters.csv")

# Plot sites

ggplot(cluster_data, aes(x = PC1, y = PC2, fill = site)) +
  geom_point(size = 3, pch = 21, alpha = 0.7) +
  theme_bw() +
  labs(title = "PCA by Site")

# Pull site sample lists

wild <- cluster_data %>% 
  filter(site == "SMM" | site == "ANF") %>% 
  pull(INDV)

write_lines(wild, "wild.txt")


ucla <- cluster_data %>% 
  filter(site == "UCLA") %>% 
  pull(INDV)

write_lines(ucla, "ucla.txt")


occ <- cluster_data %>% 
  filter(site == "OCC") %>% 
  pull(INDV)

write_lines(occ, "occ.txt")


ucsb <- cluster_data %>% 
  filter(site == "UCSB") %>% 
  pull(INDV)

write_lines(ucsb, "ucsb.txt")


ucsd <- cluster_data %>% 
  filter(site == "UCSD") %>% 
  pull(INDV)

write_lines(ucsd, "ucsd.txt")


sfsu <- cluster_data %>% 
  filter(site == "SFSU") %>% 
  pull(INDV)

write_lines(sfsu, "sfsu.txt")

# separate WL for validation

smm <- cluster_data %>% 
  filter(site == "SMM") %>% 
  pull(INDV)

write_lines(smm, "smm.txt")

anf <- cluster_data %>% 
  filter(site == "ANF") %>% 
  pull(INDV)

write_lines(anf, "anf.txt")

# For fastPHASE, pop_map.txt file with samples and sites.

pop_map <- cluster_data %>% 
  select(INDV, site) %>% 
  mutate(site = recode(site, 
                       "ANF" = "WL", 
                       "SMM" = "WL"))

write_tsv(pop_map, "pop_map.txt", col_names = FALSE)

# Pull cluster sample lists

# cluster1 <- cluster_data %>% 
#   filter(cluster == 1) %>% 
#   pull(INDV)
# 
# write_lines(cluster1, "pop1.txt")
# 
# cluster2 <- cluster_data %>% 
#   filter(cluster == 2) %>% 
#   pull(INDV)
# 
# write_lines(cluster2, "pop2.txt")
# 
# cluster3 <- cluster_data %>% 
#   filter(cluster == 3) %>% 
#   pull(INDV)
# 
# write_lines(cluster3, "pop3.txt")

