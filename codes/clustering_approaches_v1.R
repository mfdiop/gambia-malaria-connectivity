# 1️⃣ Hierarchical clustering (Ward linkage) on IBD distance

# Load packages
library(vcfR)
library(adegenet)     # for genlight object
library(SNPRelate)    # efficient IBD/IBS computation
library(ape)          # for dendrogram plotting
library(igraph)
library(mclust)   # for ARI
library(aricode)  # for NMI
library(tidyverse)

# Step 1: Read VCF
vcf <- read.vcfR("../01_objective1/01_data/raw_data/gambia.filtered.vcf.gz")
genlight_obj <- vcfR2genlight(vcf)

# Convert to GDS for SNPRelate
snpgdsVCF2GDS("../01_objective1/01_data/raw_data/gambia.filtered.vcf.gz", "filtered_data.gds")
genofile <- snpgdsOpen("filtered_data.gds")

# Step 2: Compute pairwise relatedness
# IBS (identity-by-state) or IBD (kinship-like)
ibs <- snpgdsIBS(genofile, num.thread = 4, autosome.only = FALSE)
ibd <- snpgdsIBDMLE(genofile, maf = 0.05, autosome.only = FALSE, kinship = TRUE)  # MLE-based IBD estimates

# Use IBD matrix
ibd_matrix <- ibd$kinship  # kinship as proxy for relatedness
rownames(ibd_matrix) <- ibd$sample.id
colnames(ibd_matrix) <- ibd$sample.id

# Step 3: Convert to distance matrix for clustering
dist_mat <- as.dist(1 - ibd_matrix)

# Step 4: Hierarchical clustering (Ward linkage)
hc <- hclust(dist_mat, method = "ward.D2")

# Plot dendrogram
plot(hc, cex=0.6, hang=-1, main="Hierarchical clustering (Ward) on IBD")


# 2️⃣ Louvain clustering on an IBD network
# Assume genofile already open from above

# Step 2: Build a network
# Apply threshold (e.g. kinship > 0.1 to define edges)
threshold <- 0.1
adj_matrix <- (ibd_matrix > threshold) * 1  # binary adjacency matrix
diag(adj_matrix) <- 0

# Step 3: Convert to igraph object
g <- graph_from_adjacency_matrix(adj_matrix, mode="undirected")

# Step 4: Louvain community detection
cl <- cluster_louvain(g)

# Step 5: Save / inspect results
membership <- membership(cl)
table(membership)

# Plot network with communities
plot(cl, g, vertex.size = 5, vertex.label = NA,
     main = "Louvain clustering on IBD network")

# Load metadata
metadata <- readxl::read_xlsx("../01_objective1/01_data/meta/GamMetadata_Final_imputemissingdate.xlsx") 
head(metadata)

# Make sure order matches VCF samples
metadata <- metadata[match(ibd$sample.id, metadata$SampleID), ]


# 1️⃣ Hierarchical Clustering + Metadata Concordance
# Cut tree into clusters (choose k by elbow/silhouette)
k <- 5
clusters_hc <- cutree(hc, k = k)

# Add to metadata
metadata$Cluster_HC <- clusters_hc

# ---- Concordance Tests ----
# Adjusted Rand Index
ari_village <- adjustedRandIndex(metadata$Cluster_HC, metadata$Village)
ari_year    <- adjustedRandIndex(metadata$Cluster_HC, metadata$Year)

# Normalized Mutual Information
nmi_village <- NMI(metadata$Cluster_HC, metadata$Village)
nmi_year    <- NMI(metadata$Cluster_HC, metadata$Year)

# Chi-square test for independence (Cluster vs Village)
chisq_village <- chisq.test(table(metadata$Cluster_HC, metadata$Village))
chisq_year    <- chisq.test(table(metadata$Cluster_HC, metadata$Year))

# ---- Report Results ----
cat("Hierarchical Clustering Results:\n")
cat("ARI (Village):", ari_village, "\n")
cat("ARI (Year):", ari_year, "\n")
cat("NMI (Village):", nmi_village, "\n")
cat("NMI (Year):", nmi_year, "\n")
print(chisq_village)
print(chisq_year)

# ---- Visualize clusters by metadata ----
ggplot(metadata, aes(x=Year, y=VillageCode, color=factor(Cluster_HC))) +
   geom_jitter(width=0.3, height=0.3, size=2) +
   theme_minimal() +
   labs(title="Hierarchical Clustering vs Space-Time Metadata",
        color="Cluster")

# 2️⃣ Louvain Clustering + Metadata Concordance

# Louvain clustering already done:
metadata$Cluster_Louvain <- membership(cl)

# ---- Concordance Tests ----
# Adjusted Rand Index
ari_village_louv <- adjustedRandIndex(metadata$Cluster_Louvain, metadata$VillageCode)
ari_year_louv    <- adjustedRandIndex(metadata$Cluster_Louvain, metadata$Year)

# NMI
nmi_village_louv <- NMI(metadata$Cluster_Louvain, metadata$VillageCode)
nmi_year_louv    <- NMI(metadata$Cluster_Louvain, metadata$Year)

# Chi-square test
chisq_village_louv <- chisq.test(table(metadata$Cluster_Louvain, metadata$VillageCode))
chisq_year_louv    <- chisq.test(table(metadata$Cluster_Louvain, metadata$Year))

# ---- Report Results ----
cat("Louvain Clustering Results:\n")
cat("ARI (Village):", ari_village_louv, "\n")
cat("ARI (Year):", ari_year_louv, "\n")
cat("NMI (Village):", nmi_village_louv, "\n")
cat("NMI (Year):", nmi_year_louv, "\n")
print(chisq_village_louv)
print(chisq_year_louv)

# ---- Plot network with cluster & metadata ----
plot(cl, g, vertex.size=5, vertex.label=NA,
     vertex.color=as.factor(metadata$VillageCode),
     main="Louvain Clusters with Village Metadata")

# K-means on PCA of IBS/IBD + Metadata Concordance
library(factoextra)   # for PCA + clustering visualization

# Convert to genlight for distance calculation
gl <- vcfR2genind(vcf)

### 2. Compute Pairwise IBS/IBD (IBD would need external caller results)
# IBS: simple proportion of allele sharing
ibs_mat <- propShared(gl)   # from adegenet (IBS proxy)

# Assume you also have IBD matrix already estimated, e.g. hmmIBD output
# ibd_mat <- as.matrix(read.csv("ibd_matrix.csv", row.names = 1))

### 3. PCA on IBS (or IBD)
# IBS first
pca_ibs <- prcomp(ibs_mat, scale. = TRUE)

# Keep first 10 PCs (adjust after variance inspection)
pcs <- as.data.frame(pca_ibs$x[,1:10])
pcs$sample <- rownames(pcs)

# Merge metadata
pcs <- left_join(pcs, metadata, by = c("sample" = "SampleID"))

### 4. K-means clustering
set.seed(123)

# You should decide K based on elbow/silhouette, not randomly
fviz_nbclust(pcs[,1:10], kmeans, method = "wss") + theme_minimal()
fviz_nbclust(pcs[,1:10], kmeans, method = "silhouette") + theme_minimal()

# Example: run K=3
kmeans_res <- kmeans(pcs[,1:10], centers = 5, nstart = 25)

pcs$cluster <- factor(kmeans_res$cluster)

### 5. Compare Clusters vs Metadata
# Adjusted Rand Index (ARI) for concordance
library(mclust)
ari_village <- adjustedRandIndex(pcs$cluster, pcs$VillageCode)
ari_year    <- adjustedRandIndex(pcs$cluster, pcs$Year)

cat("ARI (cluster vs village):", ari_village, "\n")
cat("ARI (cluster vs year):", ari_year, "\n")

# Cross-tabulation
table_cluster_village <- table(pcs$cluster, pcs$VillageCode)
table_cluster_year    <- table(pcs$cluster, pcs$Year)

print(table_cluster_village)
print(table_cluster_year)

### 6. Visualization

shapes <- seq(1:10)

# PCA plot with clusters
ggplot(pcs, aes(x = PC1, y = PC2, color = cluster, shape = VillageCode)) +
   geom_point(size = 3) +
   labs(title = "K-means clustering on PCA (IBS)") +
   theme_minimal()


