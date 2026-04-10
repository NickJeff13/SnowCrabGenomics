#Fst between clusters

plink \
  --bfile snowcrab.maffiltered \
  --fst \
  --within MAFfilteredClusters.txt \
  --out fst_clusters

# Cluster 1
plink \
  --bfile snowcrab.maffiltered \
  --keep cluster1_keep.txt \
  --het \
  --out het_cluster1

# Cluster 2
plink \
  --bfile snowcrab.maffiltered \
  --keep cluster2_keep.txt \
  --het \
  --out het_cluster2
