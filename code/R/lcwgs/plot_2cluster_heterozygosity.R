#Plot heterozygosity between UMAP2 axis clusters 

# Read heterozygosity output
het1 <- read.table("MAF_Filtered_Plink/het_cluster1.het", header = TRUE)
het2 <- read.table("MAF_Filtered_Plink/het_cluster2.het", header = TRUE)

# Plink's .het columns: FID, IID, O(HOM), E(HOM), N(NM), F
# Observed heterozygosity = (N(NM) - O(HOM)) / N(NM)
het1$ObsHet <- (het1$N.NM. - het1$O.HOM.) / het1$N.NM.
het2$ObsHet <- (het2$N.NM. - het2$O.HOM.) / het2$N.NM.

het1$Cluster <- "cluster1"
het2$Cluster <- "cluster2"

het_all <- rbind(het1, het2)

# Summary statistics
tapply(het_all$ObsHet, het_all$Cluster, summary)

# Plot
ggplot(het_all, aes(x = Cluster, y = ObsHet, fill = Cluster)) +
  geom_boxplot(alpha = 0.7, colour="black") +
  scale_fill_manual(values=c("grey80", "grey20"))+
  labs(title = "Observed Heterozygosity by UMAP Cluster",
       x = "Cluster", y = "Observed Heterozygosity")+
  theme_bw()+
  theme(legend.position = "none")

#replot with less labels for inset plot
het.plot <- ggplot(het_all, aes(x = Cluster, y = ObsHet, fill = Cluster)) +
  geom_boxplot(alpha = 0.7, colour="black") +
  scale_fill_manual(values=c("grey80", "grey20"))+
  labs(
       x = "Cluster", 
       y = "Observed Heterozygosity")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x=element_blank());het.plot

ggsave(filename = "~/Documents/GitHub/SnowCrabGenomics/figures/UMAP2cluster_heterozygosity_Comparison.png", plot = last_plot(),
       device="png", dpi=320, width=10, height=8)
  