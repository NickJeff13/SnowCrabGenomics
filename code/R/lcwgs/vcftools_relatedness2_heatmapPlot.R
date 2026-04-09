library(data.table)
library(dplyr)
library(ggplot2)

crabkin <- read.table("/mnt/sdb/SnowCrab_LCWGS/vcfs/allcrabs.maffiltered.relatedness2", sep="\t", header = T) %>% glimpse()

#remove some text to make easier to read
crabkin$INDV1 <- gsub(".realigned.bam","",crabkin$INDV1)
crabkin$INDV2 <- gsub(".realigned.bam","",crabkin$INDV2)

phi <- crabkin$RELATEDNESS_PHI
crabkin2 <- crabkin %>% 
  dplyr::select(INDV1, INDV2, RELATEDNESS_PHI)

crabkin2$INDV1 <- gsub("NS.1780.003.UDP0193_i7---UDP0193_i5.","", crabkin2$INDV1)
crabkin2$INDV2 <- gsub("_i5.","", crabkin2$INDV2)

df <- crabkin2

# Ensure consistent ordering
df$Var1 <- factor(df$INDV1, levels = unique(df$INDV1))
df$Var2 <- factor(df$INDV2, levels = unique(df$INDV1))

ggplot(crabkin2, aes(x = INDV1, y = INDV2, fill = RELATEDNESS_PHI)) +
  geom_tile() +
  scale_fill_gradientn(
  colours=topo.colors(5)
  ) +
  theme_minimal() +
  theme(axis.text.y=element_blank())+
  labs(
    x = "Individual",
    y = "Individual",
    fill = "Relatedness"
  )

ggsave("figures/crab_maffiltered_relatedness2Heatmap.png", plot=last_plot(),
       device = "png", width=8, height=8, dpi = 300)
