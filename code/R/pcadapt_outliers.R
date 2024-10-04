library(pcadapt)
library(qvalue)

setwd("/mnt/sdb/SnowCrab_LCWGS/")
all.snp <-read.pcadapt(input = "snowcrab.recode.bed", type = "bed")

#choose k Principal components
y<-pcadapt(all.snp, K=1)
plot(y, option="manhattan")
x <- pcadapt(all.snp, K=10)
plot(x,option="screeplot")
plot(x, option="scores")
plot(x, option="scores", i=1, j=3)
plot(x, option="manhattan")
#reduce K=3 now
xx <- pcadapt(all.snp, K=3)
summary(xx)
plot(xx, option="qqplot")
plot(xx, option="stat.distribution")

par(mfrow=c(2,2))
for (i in 1:3)
  plot(xx$loadings[,i], pch=19, cex=.3, ylab=paste0("Loadings PC",i))
