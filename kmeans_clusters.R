library(tidyverse)
library(knitr)

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

set.seed(121)

clusters <- kmeans(combo_data %>% select(PC1, PC2), centers = 3)

cluster_data <- combo_data %>%
  mutate(cluster = as.factor(clusters$cluster))

ggplot(cluster_data, aes(x = PC1, y = PC2, fill = cluster)) +
  geom_point(size = 3, pch = 21, alpha = 0.7) +
  theme_bw() +
  labs(title = "PCA with K-Means Clustering")

write_csv(cluster_data, "clusters.csv")

cluster1 <- cluster_data %>% 
  filter(cluster == 1) %>% 
  pull(INDV)

write_lines(cluster1, "pop1.txt")

cluster2 <- cluster_data %>% 
  filter(cluster == 2) %>% 
  pull(INDV)

write_lines(cluster2, "pop2.txt")

cluster3 <- cluster_data %>% 
  filter(cluster == 3) %>% 
  pull(INDV)

write_lines(cluster3, "pop3.txt")