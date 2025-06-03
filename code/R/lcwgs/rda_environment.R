##########################################################################################
############Redundancy Analyses of Genotypes and Environment############################
## Written by Nick Jeffery in winter 2025

# Load libraries ----------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")

library(adegenet)
library(hierfstat)
library(vegan)
library(psych)
library(dplyr)
library(lfmm)
library(qvalue)
library(ggplot2)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(Polychrome)

#We need 1) genotypes per individual, 2) environmental data per site and individual, and 3) GPS coordinates per individual
# We will be using the exome dataset as the full snp dataset will take a long time

########RDA Pop Plot Function########
palette <- c("Bradelle Bank" ="#FFFFFF", "Chaleurs"= "#0000FF", "CMA 10A"= "#FF0000", "CMA 10B"= "#00FF00", 
             "CMA 12C"= "#000033", "CMA 12G"= "#FF00B6", "CMA 3B"="#005300", "CMA 3D"= "#FFD300", 
             "CMA 3N200"= "#009FFF", "CMA 4" ="#9A4D42", "CMA 5A"= "#00FFBE", "CMA 6B"= "#783FC1", 
             "CMA 8A"= "#1F9698", "CMA N5440"="#FFACFD", "Fortune Bay"= "#B1CC71",  "Laurentian Chan" ="#F1085C",
             "Lilly Canyon"= "#FE8F42","Mar B"= "#DD00FF", "Mar C"= "#201A01", "Mar D"= "#720055", "Mar E"= "#766C95", 
             "NAFO 3L" ="#02AD24", "NAFO 4X" = "#C8FF00", "NENSin"= "#886C00", "NENSout"="#FFB79F", "Quebec"= "#858567", 
             "St Marys Bay"= "#A10300", "Trinity Bay"= "#14F9FF", "West Cape Breton"= "#00479E")
#ggplot
RDA_plot <- function(envdata,sitedata, xaxislab, yaxislab, nudgeX,nudgeY,r2x, r2y, r2,vjust,hjust){ #xaxis and yaxis are the RDAs to be plotted from envdata
  ggplot(sitedata, aes(x=RDA1, y=RDA2))+
    geom_hline(linewidth=0.3, colour = "grey20", yintercept = 0) + #horizontal line through the origin
    geom_vline(linewidth=0.3, colour = "grey20", xintercept = 0) + #vertical line through the origin
    # geom_point(size=3, alpha = 0.7,shape=20,colour = "grey70", data = snpdata) + #SNPs
    geom_segment(data = envdata, aes(x=0, xend=RDA1*3, y=0, yend=RDA2*3), 
                 arrow = ggplot2::arrow(length = unit(0.25, "cm")),
                 colour = "navy", linewidth = 1) +#env arrows
    geom_text_repel(data = envdata, aes(x=RDA1*3, y=RDA2*3,label = Labels), 
                    size = 5,colour = "navy",parse = T,
                    segment.color=NA,box.padding=0,
                    nudge_x =nudgeX,#however many enviro variables we have
                    nudge_y = nudgeY) + #labels for env arrows
    geom_point(data = sitedata,  aes(fill=Pop),shape=21, color="black", size=3, alpha = 0.7) + #sample  points
    #scale_shape_manual(values = shapes) +
    scale_fill_manual(values=palette)+
    #scale_colour_manual(values=palette)+
    annotate("text",x=r2x,y=r2y,label=r2, parse=T,vjust=vjust,hjust=hjust,size=5) +
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

exon.snps <- read.PLINK("/mnt/sdb/SnowCrab_LCWGS/rawexon.raw", parallel = TRUE)

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
imp.pca <- dudi.pca(exon.imp[,1:134087], scannf = F, nf = 5)
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
identical(as.factor(env.ordered$Rpop), unique(exon.ordered$pop)) #TRUE
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
                col.ind="black",
                pointsize=3) #looks like there's no real relationship between general location and environment


ggsave(filename = "PCA_EnvOnly_Black.png", plot = last_plot(), device = "png", path = "figures/", width = 10, height=8, dpi=320)

# Run the RDAs ------------------------------------------------------------
gc()
# we will run a full environmental PCA and a partial RDA controlling for lat/long

# need to merge genotype and enviro data, then split again in the rda formula 

env.snps.merged <- left_join(x = exon.ordered, y=env.ordered, by = c("pop"="Rpop"))

env.dat.inds <- env.snps.merged %>% 
  select(-contains("Locus"))

#1. Start with RDA with environment only

crab.env.rda <- rda(exon.ordered[,1:134087] ~ GLORYS.depth..m. + summer_sal + winter_temp_min + winter_temp_mean + spring_temp_mean + summer_temp_max + fall_temp_min, data=env.dat.inds, scale =T)

summary(crab.env.rda) #constrained proportion =0.007953 
RsquareAdj(crab.env.rda) #r square=0.00785, rsq adj=0.0014

signif.env.rda <- anova.cca(crab.env.rda, parallel = 40, permutations = 999) #p=0.001
anova.cca(crab.env.rda, by="terms", parallel=40)
anova.cca(crab.env.rda, by="axis", parallel=40)

#2. Run the same RDA but controlling for lat and long
crab.env.rda.latloncond <- rda(exon.ordered[,1:134087] ~ GLORYS.depth..m. + summer_sal + winter_temp_min + 
                                 winter_temp_mean + spring_temp_mean + summer_temp_max + fall_temp_min + 
                                 Condition(Lat + Long), 
                               data=env.dat.inds, scale =T)

summary(crab.env.rda.latloncond) #proportion conditioned =0.0023, proportion constrained =0.00766
RsquareAdj(crab.env.rda.latloncond) #R2 =0.00766, r2 adj =0.0011

#3. Run an RDA on just lat and long
crab.latlong.rda <- rda(exon.ordered[,1:134087] ~ Lat + Long, 
                               data=env.dat.inds, scale =T)
summary(crab.latlong.rda) #constrained=0.00229

# Plot the RDAs -----------------------------------------------------------

#extract data for ggplot
SNPPoints_GenEnvPCA_RDA <- as.data.frame(scores(crab.env.rda, scaling = 3, display = "species"))
IndivPoints_GenEnvPCA_RDA <- as.data.frame(scores(crab.env.rda, scaling = 3, display = "sites")) %>% 
  cbind(env.snps.merged %>% dplyr::select(c(ind,pop))) %>% 
  dplyr::mutate(Pop = factor(pop, levels = c("Bradelle Bank", "Chaleurs", "CMA 10A", "CMA 10B", "CMA 12C",  "CMA 12G",  "CMA 3B", "CMA 3D",          
                                             "CMA 3N200",  "CMA 4", "CMA 5A", "CMA 6B", "CMA 8A", "CMA N5440", "Fortune Bay",  "Laurentian Chan", 
                                             "Lilly Canyon","Mar B", "Mar C","Mar D", "Mar E", "NAFO 3L", "NAFO 4X", "NENSin", "NENSout",
                                             "Quebec", "St Marys Bay", "Trinity Bay", "West Cape Breton")))

                EnvArrows_GenEnvPCA_RDA <- as.data.frame(scores(crab.env.rda, scaling = 2, display = "bp")) %>%
                                                           dplyr::mutate(Labels = c("Depth","Summer_salinity","Winter_temp_min","Winter_temp_mean",
                                                                                    "Spring_temp_mean","Summer_temp_max", "Fall_temp_min"))

# use the RDA_plot function made at the start
crab.env.rda.plot1 <- RDA_plot(envdata = EnvArrows_GenEnvPCA_RDA,
                               #snpdata= SNPPoints_GenEnvPCA_RDA,
                               sitedata=IndivPoints_GenEnvPCA_RDA,
                               xaxislab="RDA 1",
                               yaxislab="RDA 2",
                               nudgeX = 0,
                               nudgeY = 0,
                               r2x=3.2,
                               r2y=2.5,
                               r2=paste("Adjusted~R^2","== 0.001"),
                               hjust=1, 
                               vjust=0)
crab.env.rda.plot1

ggsave("figures/Crab_RDA_EnvOnly.png", plot = crab.env.rda.plot1, height=8, width=12, units = "in", dpi=300)



# Identification of environment associated outliers -----------------------

rda.loadings <- scores(crab.env.rda, choices=c(1:3), display="species")
hist(rda.loadings[,1]) #normal distribution

# function from online tutorial for selecting statistical outliers based on loadings and standard deviations
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# apply function to RDA axes of interest

out1 <- outliers(rda.loadings[,1], 3) #542
out2 <- outliers(rda.loadings[,2], 3) #392
out3 <- outliers(rda.loadings[,3], 3) #452

out1 <- cbind.data.frame(rep(1,times=length(out1)), names(out1), unname(out1))
out2 <- cbind.data.frame(rep(2,times=length(out2)), names(out2), unname(out2))
out3 <- cbind.data.frame(rep(3,times=length(out3)), names(out3), unname(out3))

colnames(out1) <- colnames(out2) <- colnames(out3) <- c("axis","snp","loading")

outliers <- rbind(out1, out2, out3)
outliers$snp <- as.character(outliers$snp)
ncand <- length(out1$snp) + length(out2$snp) + length(out3$snp)

#Add correlations with SNP loading and environment to see which loci load against which vars

env.cor <- matrix(nrow=ncand, ncol=7)  # 7 columns for 7 predictors
colnames(env.cor) <- c("GLORYS.depth..m.", "summer_sal", "winter_temp_min", "winter_temp_mean", "spring_temp_mean", "summer_temp_max", "fall_temp_min")

i=NULL
for (i in 1:length(outliers$snp)) {
  nam <- outliers[i,2]
  snp.gen <- exon.imp[,nam]
  env.cor[i,] <- apply(env.dat.inds[,5:length(colnames(env.dat.inds))], 2, function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(outliers,env.cor)  
head(cand)

#check for and remove duplicates 
length(cand$snp[duplicated(cand$snp)])  # 59 duplicate SNPs
cand <- cand[!duplicated(cand$snp), ]

#Now determine which predictor each SNP is correlated with

i=NULL
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,11] <- names(which.max(abs(bar[4:10]))) # gives the variable
  cand[i,12] <- max(abs(bar[4:10]))            # gives the correlation
}

colnames(cand)[11] <- "predictor"
colnames(cand)[12] <- "correlation"

table(cand$predictor) 

# Plot the SNPs and their environmental correlates

sel <- cand$snp
env <- cand$predictor
env[env=="GLORYS.depth..m."] <- '#1f78b4'
env[env=="summer_sal"] <- '#a6cee3'
env[env=="winter_temp_min"] <- '#6a3d9a'
env[env=="winter_temp_mean"] <- '#e31a1c'
env[env=="spring_temp_mean"] <- '#33a02c'
env[env=="summer_temp_max"] <- '#ffff33'
env[env=="fall_temp_min"] <- '#fb9a99'


# color by predictor:
col.pred <- rownames(crab.env.rda$CCA$v) # pull the SNP names

i=NULL
for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("Locus",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99')

# axes 1 & 2
plot(crab.env.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(crab.env.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(crab.env.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(crab.env.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("GLORYS.depth..m.", "summer_sal", "winter_temp_min", "winter_temp_mean", "spring_temp_mean", "summer_temp_max", "fall_temp_min"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


# axes 1 & 3
plot(crab.env.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(crab.env.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(crab.env.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(crab.env.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend= c("GLORYS.depth..m.", "summer_sal", "winter_temp_min", "winter_temp_mean", "spring_temp_mean", "summer_temp_max", "fall_temp_min"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 2 & 3
plot(crab.env.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3))
points(crab.env.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(2,3))
points(crab.env.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(2,3))
text(crab.env.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))
legend("bottomright", legend= c("GLORYS.depth..m.", "summer_sal", "winter_temp_min", "winter_temp_mean", "spring_temp_mean", "summer_temp_max", "fall_temp_min"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# Save data ---------------------------------------------------------------

save.image(file = "Crab_Enviro_RDA.RData")
