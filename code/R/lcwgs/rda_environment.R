##########################################################################################
############Redundancy Analyses of Genotypes and Environment############################
## Written by Nick Jeffery in winter 2025

# Load libraries ----------------------------------------------------------
library(adegenet)
library(hierfstat)
library(vegan)
library(psych)
library(dplyr)
library(lfmm)
library(qvalue)
library(ggplot2)
library(FactoMineR)
library(factoextra)

#We need 1) genotypes per individual, 2) environmental data per site and individual, and 3) GPS coordinates per individual
# We will be using the exome dataset as the full snp dataset will take a long time

########RDA Pop Plot Function########
#ggplot
RDA_plot <- function(envdata,sitedata, xaxislab, yaxislab, nudgeX,nudgeY,r2x, r2y, r2,vjust,hjust){ #xaxis and yaxis are the RDAs to be plotted from envdata
  ggplot(sitedata, aes(x=RDA1, y=RDA2))+
    geom_hline(linewidth=0.3, colour = "grey20", yintercept = 0) + #horizontal line through the origin
    geom_vline(linewidth=0.3, colour = "grey20", xintercept = 0) + #vertical line through the origin
    # geom_point(size=3, alpha = 0.7,shape=20,colour = "grey70", data = snpdata) + #SNPs
    geom_segment(data = envdata, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 arrow = ggplot2::arrow(length = unit(0.25, "cm")),
                 colour = "navy", linewidth = 0.4) +#env arrows
    geom_text_repel(data = envdata, aes(x=RDA1, y=RDA2,label = Labels), 
                    size = 3,colour = "navy",parse = T,
                    segment.color=NA,box.padding=0,
                    nudge_x =nudgeX,#PC1,PC2,Linf
                    nudge_y = nudgeY) + #labels for env arrows
    geom_point(data = sitedata, size=1, alpha = 0.7, aes(shape=Pop,fill=Pop,  colour=Pop)) + #sample  points
    scale_shape_manual(values = shapes) +
    scale_fill_manual(values=palette)+
    scale_colour_manual(values=palette)+
    annotate("text",x=r2x,y=r2y,label=r2, parse=T,vjust=vjust,hjust=hjust,size=3) +
    xlab(xaxislab) +
    ylab(yaxislab) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme_bw() +
    theme(plot.title=element_text(size = 10),
          legend.position = "right", 
          legend.title=element_blank(),
          legend.text = element_text(size=12),
          legend.key = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(r=20),
          axis.text = element_text(size = 8), 
          axis.title.x = element_text(size=10, vjust =3),
          axis.title.y = element_text(size=10, vjust=-3),
          panel.border = element_rect(linewidth =0.5),
          plot.margin = unit(c(1,1,1,0), 'lines'), #t,r,b,l
          panel.grid=element_blank()) + #t,r,b,l
    NULL
}


# Read in plink genotypes -------------------------------------------------
# these large files are on my local Linux computer so the file paths only work for me, but the data will be available on NCBI and other sources later
# first converted imputed exon snps file to RAW with the following stats: 1072 individuals, 134,087 variants passed QC, genotyping rate is 0.9162

exon.snps <- read.PLINK("/mnt/sdb/SnowCrab_LCWGS/vcfs/rawexon.raw", parallel = TRUE)

# find clusters while we have this loaded to see if there are any
#Note, this took 3 days without finishing so we'll skip it
#grp <- find.clusters(exon.snps, n.pca = 10, method = "kmeans", scale = FALSE)

dim(exon.snps)
sum(is.na(exon.snps)) # zero NAs so we're good to use this for RDA 

exon.snps@ind.names <- gsub(".realigned.bam", "", exon.snps@ind.names)

exon.snps@loc.names <- c(paste0("Locus","_",1:length(exon.snps$loc.names)))

poplist <- as.factor(c(rep("Mar C",45), rep("Mar E", 45), rep("Mar B",33), rep("Quebec", 33), rep("Lilly Canyon", 22), rep("CMA 3B", 34), rep("CMA 6B", 34), 
                   rep("CMA N5440", 21), rep("CMA 10B",6), "Chaleurs", rep("CMA 10B", 7), "Chaleurs", rep("CMA 10B", 7), "Chaleurs", rep("CMA 10B", 7),
                   "Chaleurs",rep("CMA 10B", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs",
                   rep("CMA 4", 7), "Chaleurs",rep("CMA 4", 6), rep("Chaleurs",12), rep("CMA 5A", 33), rep("CMA 8A", 32), rep("Lilly Canyon",9), 
                   rep("Chaleurs", 4), rep("CMA 3N200", 10), rep("NAFO 3L", 33), rep("CMA 3D", 35), rep("CMA 10A", 34), rep("St Marys Bay", 27),
                   rep("Fortune Bay", 35), rep("Trinity Bay", 13), rep("West Cape Breton", 8), rep("Bradelle Bank", 44), rep("CMA 3N200", 15),
                   rep("CMA N5440", 3), "St Marys Bay", rep("CMA N5440", 8), rep("West Cape Breton", 4), rep("Trinity Bay", 20), rep("NENSin", 25), 
                   rep("NAFO 4X", 32), rep("NENSout", 35), rep("CMA 12G", 33), rep("CMA 12C", 35), rep("Laurentian Chan", 35), rep("West Cape Breton", 54), 
                   rep("Mar D", 8), rep("West Cape Breton", 65), rep("Mar D", 26)) )
#the first few reps of West Cape Breton are female, and the last 65 are the males

exon.snps@pop <- poplist
sum(is.na(exon.snps))

exon2 <- as.data.frame.genlight(exon.snps)

#impute NAs
exon.imp <-  as.data.frame(apply(exon2, 2, 
                                 function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))))

exon.imp$ind <- rownames(exon.imp)
exon.imp$pop <- poplist
imp.pca <- dudi.pca(exon.imp, scannf = F, nf = 5)
scatter(imp.pca)

#make the genetic data in order by pop
exon.ordered <- exon.imp[order(exon.imp$pop), ]
sumstats <- adegenet::summary(exon.imp)


# Read in environmental data ----------------------------------------------
# here we have seasonal temperature and salinity and depth data from the Digital Twins of the Ocean GLORYS model, and coordinates 

env <- read.csv("data/DTO extractions/Pop_Coords2024_with_GLORYS_info.csv", sep="\t", header = T) %>% glimpse()

pairs.panels(env[,12:length(colnames(env))], scale=T)
#this plot shows high correlations between all salinity variables, unsurprisingly, so we'll just select one. Also drop spring temp min and summer temp mean and fall temp max

env.filt <- env[,c(3:5,10,13,16:17, 19, 21,22)]

env.ordered <- env.filt[order(env.filt$Rpop), ]

#make sure the pops are in the same order between the genetics and the env dataframe
identical(env.ordered$Rpop, unique(exon.ordered$pop)) #TRUE
rownames(env.ordered) <- env.ordered$Rpop

# Run a PCA on the environment data
crab_enviro_PCA <- dudi.pca(env.ordered %>% dplyr::select(c("GLORYS.depth..m.", 
                                                            "summer_sal", 
                                                            "winter_temp_min", 
                                                            "winter_temp_mean", 
                                                            "spring_temp_mean",
                                                            "summer_temp_max",
                                                            "fall_temp_min")),
                                               center=T,scale=T,
                                               scannf=F, nf=4) #view screeplot to see how many axes to keep
#use factoextra to get scores etc and make a nice plot
eigs <- get_eigenvalue(crab_enviro_PCA)

#scree plot and biplot of results
fviz_eig(crab_enviro_PCA, addlabels = T)

fviz_pca_var(crab_enviro_PCA)

ind <- get_pca_ind(crab_enviro_PCA)

fviz_pca_ind(crab_enviro_PCA, repel = T,
             pointsize=3,
             pointshape =21, fill="blue")

fviz_pca_biplot(crab_enviro_PCA, repel=T,
                col.var="#2E9FDF",
                col.ind="#696969",
                pointsize=3) #looks like there's no real relationship between general location and environment



# Run the RDAs ------------------------------------------------------------

# we will run a full environmental PCA and a partial RDA controlling for lat/long

# need to merge genotype and enviro data, then split again in the rda formula 

env.snps.merged <- left_join(x = exon.ordered, y=env.ordered, by = c("pop"="Rpop"))
env.dat.inds <- env.snps.merged %>% 
  select(-contains("Locus"))

#1. Start with RDA with environment only

crab.env.rda <- rda(exon.ordered[,1:134087] ~ GLORYS.depth..m. + summer_sal + winter_temp_min + winter_temp_mean + spring_temp_mean + summer_temp_max + fall_temp_min, data=env.dat.inds, scale =T)

summary(crab.env.rda) #constrained proportion =0.007953 
RsquareAdj(crab.env.rda) #r square=0.00785, rsq adj=0.0014

signif.env.rda <- anova.cca(crab.env.rda, parallel = 40, permutations = 999)
anova.cca(crab.env.rda, by="terms", parallel=40)
anova.cca(crab.env.rda, by="axis", parallel=40)

#2. Run the same RDA but controlling for lat and long


#3. Run an RDA on just lat and long

