
# Clustering + Evaluation Metrics
# 
# Clustering evaluation metrics fall into two categories: internal metrics, which assess cluster quality (compactness and separation) 
# without ground truth labels (e.g., Silhouette Score, Davies-Bouldin Index, Dunn Index, Inertia), and external metrics, 
# which compare clustering results to known labels to measure agreement (e.g., Adjusted Rand Index, V-Measure, Normalized Mutual Information). 
# Key internal metrics include the Silhouette Score (higher is better) and Davies-Bouldin Index (lower is better), 
# while Adjusted Rand Index and V-Measure are popular external metrics
# 
# Determining The Optimal Number Of Clusters: 3 Must Know Methods (elbow, silhouette and gap statistic)
# Determining the optimal number of clusters (k) in a data set is a fundamental issue in partitioning clustering, 
# such as k-means clustering, which requires the user to specify the number of clusters k to be generated.
# 
# Computing the number of clusters using R
# 1. fviz_nbclust() function [in factoextra R package]: It can be used to compute the three different methods 
# [elbow, silhouette and gap statistic] for any partitioning clustering methods [K-means, K-medoids (PAM), CLARA, HCUT (factoextra package)]
# 
# 2. NbClust() function [ in NbClust R package] (Charrad et al. 2014): It provides 30 indices for determining the relevant number of clusters 
# and proposes to users the best clustering scheme from the different results obtained by varying all combinations of number of clusters, 
# distance measures, and clustering methods.

# ==============================
#        Load libraries
# ==============================
library(vcfR)
library(tidyverse)
library(NbClust)      # Compute the number of clusters
library(cluster)      # silhouette, PAM, clustering
library(factoextra)   # visualization
library(ape)          # cophenetic correlation
library(kernlab)      # spectral clustering (specc)
library(proxy)        # distance calculations
library(stats)
library(mclust)       # adjusted rand index
library(aricode)      # NMI
library(clusterSim)
library(patchwork)


# ===============================
#         Functions
# ==============================        
ibs_pair <- function(g1, g2) {
   valid <- !is.na(g1) & !is.na(g2)
   g1 <- g1[valid]; g2 <- g2[valid]
   if(length(g1) == 0) return(NA)
   sum(g1 == g2, na.rm = TRUE)/length(g1)
}

get_ibs_matrix <- function(vcf, ignore_het = FALSE, impute = FALSE){
   
   # gt <- extract.gt(vcf, element = "GT")
   # samples <- colnames(gt)
   
   # Extract the Read Count and the Coverage
   ad <- vcfR::extract.gt(vcf, element = "AD")
   dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = T)
   
   # Estimate Within Sample Allele Frequencies (WSAFs)
   altad <- vcfR::masplit(ad, record = 2, sort = F)
   wsaf <- altad/dp
   
   # gt <- extract_gt_tidy(vcf)
   # 
   # names <- c("Key", "Indiv", "gt_AD", "gt_DP", "gt_GT", "gt_GT_alleles")
   # 
   # if(names %in% colnames(gt)) 
   #    gt <- gt %>% select(names)
   
   # geno_numeric <- apply(gt, 2, function(x) {
   #    gsub("\\|", "/", x)
   #    g <- suppressWarnings(as.numeric(gsub("/", "", x)))
   #    g[is.na(g)] <- NA
   #    return(g)
   # })
   # 
   # geno_numeric <- t(geno_numeric)
   # rownames(geno_numeric) <- colnames(gt)
   
   
   # ==============================
   #   deal with heterozygous sites
   # ==============================
   if (ignore_het) {
      wsaf[wsaf != 0 & wsaf != 1] <- NA
   } else {
      wsaf <- round(wsaf)
   }
   
   # ==============================
   #     Impute missing values
   # ==============================
   if(impute) {
      locus_impute <- apply(wsaf, 2, median, na.rm = TRUE)
      locus_impute <- outer(rep(1, nrow(wsaf)), locus_impute)
      wsaf[is.na(wsaf)] <- locus_impute[is.na(wsaf)]
      
      # Transpose the matrix to have individuals on rows
      wsaf_impute <- as.matrix(t(wsaf))
      
      rm(locus_impute, wsaf)
   } else {
      # Transpose the matrix to have individuals on rows
      wsaf_impute <- as.matrix(t(wsaf))
   }
   
   samples <- colnames(ad)
   n <- nrow(wsaf_impute)
   ibs_mat <- matrix(NA, n, n)
   rownames(ibs_mat) = colnames(ibs_mat) <- samples
   
   pb <- txtProgressBar(min = 0, max = n, style = 3) # Initialize progress bar
   
   for (i in 1:n) {
      for (j in i:n) {
         ibs_val <- round(ibs_pair(wsaf_impute[i,], wsaf_impute[j,]), 3)
         ibs_mat[i,j] <- ibs_val
         ibs_mat[j,i] <- ibs_val
      }
      Sys.sleep(0.1) # Simulate task (replace with your code)
      setTxtProgressBar(pb, i) # Update progress bar
   }
   
   close(pb) # Close the progress bar
   
   return(list(ibs = ibs_mat, wsaf = wsaf_impute))
}


# =================================
# 1. Load data (VCF) and metadata
# =================================
vcf <- read.vcfR("../01_objective1/01_data/raw_data/gambia.filtered.vcf.gz")
meta <- readxl::read_xlsx("../01_objective1/01_data/meta/GamMetadata_Final_imputemissingdate.xlsx")   # columns: sample, village, year

# ========================================
# 2. Build genetic distance or IBS matrix
# ========================================
# Option 1: use existing distance matrix (e.g., SNP distances)
# dist_mat <- as.dist(read.csv("genetic_distance_matrix.csv", row.names = 1))

# Option 2: compute IBS as fallback
results <- get_ibs_matrix(vcf, impute = TRUE)

# Convert IBS similarity to distance
dist_mat <- as.dist(1 - results$ibs)

samples <- colnames(results$ibs)

meta <- meta %>% 
   filter(SampleID %in% samples)

# ======================================
# 3. Computing the number of clusters
# ======================================
# a. fviz_nbclust() function: Elbow, Silhouhette and Gap statistic methods

# Elbow method
p1 <- fviz_nbclust(as.data.frame(results$ibs), kmeans, method = "wss") +
   geom_vline(xintercept = 4, linetype = 2)+
   labs(subtitle = "Elbow kmeans method")

p2 <- fviz_nbclust(as.data.frame(results$ibs), hcut, method = "wss") +
   labs(subtitle = "Elbow hcut method")

p3 <- fviz_nbclust(as.data.frame(results$ibs), pam, method = "wss") +
   labs(subtitle = "Elbow pam method")

# Shows k = 3 or 4

# Silhouette method
p4 <- fviz_nbclust(as.data.frame(results$ibs), kmeans, method = "silhouette") +
   labs(subtitle = "Silhouette kmeans method")

p5 <- fviz_nbclust(as.data.frame(results$ibs), hcut, method = "silhouette") +
   labs(subtitle = "Silhouette hcut method")

p6 <- fviz_nbclust(as.data.frame(results$ibs), pam, method = "silhouette") +
   labs(subtitle = "Silhouette pam method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
p7 <- fviz_nbclust(as.data.frame(results$ibs), kmeans, nstart = 25, k.max = 20, method = "gap_stat", nboot = 100)+
   labs(subtitle = "Gap statistic kmeans method")

p8 <- fviz_nbclust(as.data.frame(results$ibs), hcut, nstart = 25, k.max = 20, method = "gap_stat", nboot = 100)+
   labs(subtitle = "Gap statistic hcut method")

p9 <- fviz_nbclust(as.data.frame(results$ibs), pam, nstart = 25, k.max = 20, method = "gap_stat", nboot = 100)+
   labs(subtitle = "Gap statistic pam method")

# b. NbClust() function: 30 indices for choosing the best number of clusters
# nb <- NbClust(results$wsaf, distance = "euclidean", min.nc = 2,
#               max.nc = 10, method = "complete", index = "all" )
# 
# # Visualize the result
# fviz_nbclust(nb) + theme_minimal()

combined <- gridExtra::grid.arrange(p1, p4, p7, 
                        p2, p5, p8, 
                        p3, p6, p9, nrow = 3)

# Save optimal number of clusters plot
ggsave("../benchmarking_transmission_methods/results/optimal_k.pdf", 
       plot = combined, width = 12, height = 10, dpi = 600)

# ======================
# Run PCA on IBS matrix
# ======================
pca <- prcomp(results$ibs, scale. = TRUE)
pcs <- as.data.frame(pca$x[,1:10]) # first 10 PCs

# =========================================
# Remove constant or zero-variance columns
# ===========================================
wsaf_clean <- results$wsaf[, apply(results$wsaf, 2, function(x) sd(x) != 0)]

# =========================================================
#    Clustering using eclust function (factorextra package)
# =========================================================
# A. Within-Sample Allele Frequencies (WSAF)
# K-means clustering
res.km <- eclust(wsaf_clean, "kmeans", k = k,
                 nstart = 25, graph = FALSE)

plot_km <- fviz_cluster(res.km, geom = "point", ellipse = FALSE, pointsize = 3, shape = 19,
                        show.clust.cent = TRUE, main = "K-Means Clustering (WSAF)", legend = 'none')

# Enhanced hierarchical clustering
res.hc <- eclust(wsaf_clean, "hclust", k = k,
                 nstart = 25, graph = FALSE)

plot_hc <- fviz_cluster(res.hc, geom = "point", ellipse = FALSE, pointsize = 3, shape = 19,
                        show.clust.cent = TRUE, main = "hierarchical Clustering (WSAF)", legend = 'none')

cluster <- as.data.frame(xx) %>% rownames_to_column("SampleID") %>% 
   inner_join(., meta, by = "SampleID")

# Dendrogram
fviz_dend(res.hc, rect = TRUE, show_labels = TRUE) 

# PAM clustering
res.pam <- eclust(wsaf_clean, "pam", k = k,
                  nstart = 25, graph = FALSE)

plot_pam <- fviz_cluster(res.pam, geom = "point", ellipse = FALSE, pointsize = 3, shape = 19,
                         show.clust.cent = TRUE, main = "PAM Clustering (WSAF)", legend = 'none')

# B. Principal Component Analysis on IBS (PCS)
# K-means clustering
km.res <- eclust(pcs, "kmeans", k = k,
                 nstart = 25, graph = FALSE)

km_pcs <- fviz_cluster(km.res, geom = "point", ellipse = FALSE, pointsize = 3, shape = 19,
                        show.clust.cent = TRUE, main = "K-Means Clustering (PCs)", legend = 'none')

# Enhanced hierarchical clustering
hc.res <- eclust(pcs, "hclust", k = k,
                 nstart = 25, graph = FALSE)

hc_pcs <- fviz_cluster(hc.res, geom = "point", ellipse = FALSE, pointsize = 3, shape = 19,
             show.clust.cent = TRUE, main = "Hierarchical Clustering (PCs)", legend = 'none')

# Dendrogram
fviz_dend(hc.res, rect = TRUE, show_labels = TRUE) 

# PAM clustering
pam.res <- eclust(pcs, "pam", k = k,
                  nstart = 25, graph = FALSE)

pam_pcs <- fviz_cluster(pam.res, geom = "point", ellipse = FALSE, pointsize = 3, shape = 19,
             show.clust.cent = TRUE, main = "PAM Clustering (PCs)") +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal", 
         legend.title = element_blank(),
         legend.text = element_text(size = 12, face = 'bold'))


# Function to extract legend from a ggplot object
get_legend <- function(a.gplot){
   tmp <- ggplot_gtable(ggplot_build(a.gplot))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   return(legend)
}

# Extract legend from one of the plots
shared_legend <- get_legend(pam_pcs)

pam_pcs <- pam_pcs + theme(legend.position = "none")

cl_combined <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_km, plot_hc, plot_pam, 
                                       km_pcs, hc_pcs, pam_pcs, ncol = 3),
                                       shared_legend,
                                       nrow = 2,
                                       heights = c(10, 1)  # adjust legend height
)

# Save optimal number of clusters plot
ggsave("../benchmarking_transmission_methods/results/clustering.pdf", 
       plot = cl_combined, width = 12, height = 10, dpi = 600)

# ==================================
# 6. Clustering Evaluate/Validation 
# ==================================
# a. silhouette
sil_km <- silhouette(kmeans_res$cluster, dist(pcs))
sil_pam <- silhouette(pam_res$clustering, dist(pcs))
cat("\nK-means avg silhouette:", round(mean(sil_km[,3]), 2), "\n")
cat("PAM avg silhouette:", round(mean(sil_pam[,3]), 2), "\n")

# b. Davies Bouldin's Index (DBI) using centroids
db_index_centroids <- index.DB(iris_data, kmeans_result$cluster, centrotypes = "centroids")
print(db_index_centroids)

# c. Dunn Index (cluster.stats function)
clustering_indices <- cluster.stats(dist(kmeans_data), kmeans_model$cluster)
# Extract the Dunn index from the clustering indices list
dunn_index <- clustering_indices$dunn
cat("Dunn Index:", dunn_index, "\n")

# ==============================
# 4. Hierarchical clustering
# ==============================
linkages <- c("ward.D2", "complete", "average", "single")
hc_results <- list()

for (l in linkages) {
   hc <- hclust(dist_mat, method = l)
   hc_results[[l]] <- hc
   
   # Cophenetic correlation
   coph_cor <- cor(cophenetic(hc), dist_mat)
   
   # Silhouette score (cut tree at k clusters, e.g., k=3)
   k <- 4
   clusters <- cutree(hc, k=k)
   sil <- silhouette(clusters, dist_mat)
   avg_sil <- mean(sil[, 3])
   
   cat("\nLinkage:", l, "\nCophenetic corr:", round(coph_cor,3),
       "\nAvg silhouette (k=", k, "):", round(avg_sil,3), "\n")
   
   # Evaluate clustering against metadata (With Ground Truth)
   ARI <- adjustedRandIndex(clusters, as.factor(meta$VillageCode))
   NMI_val <- NMI(clusters, as.factor(meta$VillageCode))
   
   cat("\nLinkage:", l, "\nAdjusted Rand Index for (k=", k, "):", round(ARI, 3), "\n")
   cat("\nLinkage:", l, "\nNormalised Mutual Information (k=", k, "):", round(NMI_val, 3), "\n")
   
   plot(sil, main="Silhouette Analysis for IBS Clustering")
   
   # Plot dendrogram
   plot(hc, main=paste("Hierarchical Clustering -", l), cex=0.6, hang=-1) # hang=-1
}

# The Dunn Index can be calculated in R using packages like clValid or fpc.
# Assuming 'dist_matrix' is your distance matrix and 'cluster_assignments' is a vector of cluster memberships
# from a clustering algorithm (e.g., k-means, hierarchical clustering).

# Using clValid (if you've run clValid already)
# results <- clValid(data, nClust = 2:5, clMethods = "kmeans", validation = "internal")
# dunn_index_values <- results@Dunn

# Using fpc
# library(fpc)
# stats <- cluster.stats(d = dist_matrix, clustering = cluster_assignments)
# dunn_index <- stats$dunn

# ==============================
# 5. K-means / K-medoids
# ==============================
# Dimensionality reduction (PCA on IBS matrix)
pca <- prcomp(results$ibs, scale. = TRUE)
pcs <- as.data.frame(pca$x[,1:10]) # first 10 PCs

# a. K-means clustering
kmeans_res <- kmeans(pcs, centers=k, nstart=25)
fviz_cluster(kmeans_res, data=pcs)


kmeans_wsaf <- kmeans(wsaf_clean, centers=k, nstart=25)
fviz_cluster(kmeans_wsaf, data=wsaf_clean, geom = "point", repel = TRUE, 
             ellipse = FALSE, pointsize = 3)

# b. K-medoids / PAM
pam_res <- pam(pcs, k=k)
fviz_cluster(pam_res, geom="point", data=pcs, repel = TRUE, 
             ellipse = FALSE, pointsize = 3)

pam_wsaf <- pam(wsaf_clean, k=k)
fviz_cluster(pam_wsaf, geom="point", data=wsaf_clean, repel = TRUE, 
             ellipse = FALSE, pointsize = 3)

# ==============================
# 5. Spectral clustering
# ==============================
# specc requires similarity matrix (use IBS)
# Normalize IBS to 0-1 if needed
ibs_sim <- results$ibs
diag(ibs_sim) <- 1  # ensure diagonal = 1
spectral_res <- specc(as.matrix(ibs_sim), centers=k)

# Evaluate silhouette
sil_spec <- silhouette(spectral_res@.Data, dist_mat)
cat("\nSpectral clustering avg silhouette:", round(mean(sil_spec[,3]), 2), "\n")

# Visualization
pcs$hc_cluster <- cutree(hc_results[["ward.D2"]], k=k)
pcs$kmeans_cluster <- kmeans_res$cluster
pcs$pam_cluster <- pam_res$clustering
pcs$spectral_cluster <- spectral_res@.Data

ggplot(pcs, aes(x=PC1, y=PC2, color=factor(hc_cluster))) +
   geom_point(size=2) + theme_minimal() + labs(title="Hierarchical Clustering (Ward)")

ggplot(pcs, aes(x=PC1, y=PC2, color=factor(kmeans_cluster))) +
   geom_point(size=2) + theme_minimal() + labs(title="K-means Clustering")

ggplot(pcs, aes(x=PC1, y=PC2, color=factor(pam_cluster))) +
   geom_point(size=2) + theme_minimal() + labs(title="PAM Clustering")

ggplot(pcs, aes(x=PC1, y=PC2, color=factor(spectral_cluster))) +
   geom_point(size=2) + theme_minimal() + labs(title="Spectral Clustering")





# Use the function eclust() [in factoextra] which provides several advantages as described in the previous chapter: Visual Enhancement of Clustering Analysis.
# 
# eclust() stands for enhanced clustering. It simplifies the workflow of clustering analysis and, it can be used to compute hierarchical clustering and partititioning clustering in a single line function call.


library(clValid)
# Iris data set:
# - Remove Species column and scale
df <- scale(iris[, -5])

# Compute clValid
clmethods <- c("hierarchical","kmeans","pam")
intern <- clValid(df, nClust = 2:6, 
                  clMethods = clmethods, validation = "internal")
# Summary
summary(intern)















