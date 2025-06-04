# Compare SNPs between poolseq and lcwgs 

pool.snps <- read.table("/mnt/sdb/Snow_Crab_PoolSeq/Trimmed/PoolSeq_SNPinfo.txt", header = T)
head(pool.snps)
pool.snps$SNP_POS <- paste0(pool.snps$crabpools.snp.info.Chromosome, sep="_", pool.snps$crabpools.snp.info.Position)

ind.snps <- read.table("/mnt/sdb/SnowCrab_LCWGS/indiv_snps.tsv", sep = "\t", header = F)
head(ind.snps)
ind.snps$SNP_POS <- paste0(ind.snps$V1, sep="_", ind.snps$V2)

#find intersecting snps
int.snps <- dplyr::intersect(pool.snps$SNP_POS, ind.snps$SNP_POS)
dim(pool.snps) #2829820 5
dim(ind.snps) #8384961 3
length(int.snps) #1834608

write.table(int.snps, file="/home/mcrg/Documents/GitHub/SnowCrabGenomics/data/Poolseq_LCWGS_IntersectingSNPS.txt", quote = F, row.names = F)
