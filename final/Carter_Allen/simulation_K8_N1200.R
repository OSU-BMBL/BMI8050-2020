# simulation of spatially organized scRNA data clustering
library(ggplot2)
library(umapr)
library(scran)
library(igraph)
library(mvtnorm)
library(mlsbm)
library(Seurat)
library(pheatmap)
library(sbmlhelpers)
library(viridis)
library(Giotto)

rm(list = ls())
#setwd("~/Documents/Research/Code/spatial")

source("switch_labels.R")
source("f_dist.R")
source("cluster_sm.R")

# define parameters 
n <- 1200 # number of cells
G <- 800 # number of genes
M <- 400 # number of genes that are markers for cell types
K <- 8 # number of cell types

# sample cluster labels
pi <- rep(1/K,K)
z <- sample(1:K,
            size = n,
            replace = TRUE,
            prob = pi)
save(z,file = "z_true.RData")

# sample the marker genes for each cell type
g <- rep(0,G)
# these are the M of G genes that are markers
gm <- sample(1:G,
             size = M,
             replace = FALSE)
g[gm] <- sample(1:K,
                size = M,
                replace = TRUE,
                prob = pi)


# set up scRNA expression matrix
counts <- matrix(0,nrow = n,ncol = G)
for(i in 1:n) # loop through cells
{
    for(j in 1:G) # loop through genes
    {
        if(z[i] == g[j]) # for the genes that are markers for cluster z_i
        {
            # sample 0 count with low probability
            # sample signal count with high probability
            pz1 = 0.02 # probability of zero count
            # sample from negative binomial
            # mu = mean
            # size = dispersion
            # higher size -> lower dispersion/lower variance
            # sample off-diagonal mixing, i.e., blurring of block structure
            n_o <- rpois(1,2)
            js <- (max(1,j-n_o)):(min(j+n_o,G))
            counts[i,js] = ifelse(rbinom(1,1,pz1) == 0,
                                  rnbinom(n = 1, size = 0.5, mu = 55),
                                  0)
        }
        else
        {
            pz2 = 0.05
            counts[i,j] = ifelse(rbinom(1,1,pz2) == 0,
                                 rnbinom(n = 1, size = 0.5, mu = 15),
                                 0)
        }
    }
}
rownames(counts) <- paste("cell",1:nrow(counts),sep = "_")
colnames(counts) <- paste("gene",1:ncol(counts),sep = "_")
tcounts <- t(counts)

counts_re <- counts[order(z),order(g)]
save(counts, file = "counts.RData")

# add row names and cell names

count_heatmap_true <- pheatmap(counts_re,
                               cluster_rows = FALSE,
                               cluster_cols = FALSE,
                               show_colnames = FALSE,
                               show_rownames = FALSE)
ggsave(plot = count_heatmap_true,
       filename = "count_heatmap_true.pdf",
       height = 10,
       width = 10)
ggsave(plot = count_heatmap_true,
       filename = "count_heatmap_true.jpg",
       height = 10,
       width = 10)

count_heatmap <- pheatmap(counts,
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          show_colnames = FALSE,
                          show_rownames = FALSE,
                          grey.colors(10))
ggsave(plot = count_heatmap,
       filename = "count_heatmap.pdf",
       height = 10,
       width = 10)

# Seurat analysis
tcounts <- t(counts)
rownames(tcounts) <- paste0("g",1:G)
colnames(tcounts) <- paste0("c",1:n)
# seur <- CreateSeuratObject(counts = tcounts)
# seur <- NormalizeData(seur)
# seur <- FindVariableFeatures(seur)
# seur <- ScaleData(seur)
# seur <- RunPCA(seur)
# seur <- FindNeighbors(seur)
#seur <- FindClusters(seur,algorithm = 4)
# seur <- RunUMAP(seur, features = 1:100)
# DimPlot(seur, reduction = "umap")
# 
# z_seur <- as.numeric(seur$seurat_clusters)
# z_seur_s <- switch_labels(z,z_seur)
# table(z,z_seur_s)
# 
# A1s <- as.matrix(seur@graphs$RNA_nn)
# fit1s <- fit_sbm(A1s,K)
# fit1s_zs <- switch_labels(z,fit1s$z)
# table(z,fit1s_zs)



emb <- umap(counts)
UMAP <- emb[,c("UMAP1","UMAP2")]

umap_plot <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z))) + 
    geom_point() + 
    theme(legend.position = "none")
ggsave(plot = umap_plot,
       filename = "umap_plot.pdf",
       width = 11,
       height = 8)

X1 <- as.matrix(UMAP)
D1 <- as.matrix(dist(X1))
G1 <- buildKNNGraph(D1,floor(sqrt(n)))
A1 <- as_adjacency_matrix(G1,sparse = FALSE)
G1s <- graph_from_adjacency_matrix(A1,mode = "undirected")

choose_K_svt(A1)

# louv1_R <- NetworkToolbox::louvain(A1)
# louv1 <- cluster_louvain(G1s)
# z1 <- louv1$membership
# z1 <- cluster_walktrap(G1) %>%
#     cut_at(no = K)
# z1 <- switch_labels(z,z1)
# table(z,z1)

# simulate cell locations using a bivariate normal distribution
MU <- matrix(c(1,3,
               2,3,
               3,3,
               4,3,
               1,2,
               2,2,
               3,2,
               4,2),
             nrow = 8,
             ncol = 2,
             byrow = TRUE)
SIG <- matrix(c(0.1,0,
                0,0.1),
              nrow = 2,
              byrow = TRUE)

coords <- matrix(0,nrow = n,ncol = 2)
for(i in 1:n)
{
    zi = z[i]
    mui = MU[zi,]
    xy = rmvnorm(1,mean = mui,sigma = SIG)
    coords[i,] = xy
}
coords <- as.data.frame(coords)
colnames(coords) <- c("X","Y")
rownames(coords) <- rownames(counts)
save(coords, file = "coords.RData")

spatial_plot <- ggplot(data = coords, aes(x = X, y = Y, color = as.factor(z))) + 
    geom_point() + 
    theme_void() +
    theme(legend.position = "none") 
ggsave(plot = spatial_plot,
       filename = "spatial_plot.pdf",
       width = 10,
       height = 10)
ggsave(plot = spatial_plot,
       filename = "spatial_plot.jpg",
       width = 10,
       height = 10)

save(counts,file = "counts_K8_N1200.RData")
save(coords,file = "coords_K8_N1200.RData")

write.csv(counts,
          file = "counts_K8_N1200.csv",
          row.names = FALSE,
          quote = FALSE)
write.csv(coords,
          file = "coords_K8_N1200.csv",
          row.names = FALSE,
          quote = FALSE)


# Giotto analysis
# giot <- createGiottoObject(raw_exprs = tcounts,
#                            spatial_locs = coords)
# giot <- normalizeGiotto(giot)
# giot_sn <- createSpatialNetwork(giot)
# save(giot_sn, file = "giot_K6_N600.RData")
# giot_sgs <- binSpect(giot_sn)
# giot_hmrf <- doHMRF(giot_sn, spatial_genes = giot_sgs$genes, k =6)

#giot_sn <- doLeidenCluster(giot_sn)
#giot_sn <- addHMRF(giot_sn,HMRFoutput = giot_hmrf,k = 6)
#spatPlot(giot_sn,)
#viewHMRFresults(giot_sn,giot_hmrf,k = 6)

# 1) PCA reduction

# 2) UMAP reduction clustering

X2 <- as.matrix(coords)
D2 <- as.matrix(dist(X2))
G2 <- buildKNNGraph(D2,floor(sqrt(n)))

G3 <- graph_from_adjacency_matrix(A1 * f_dist(D2), 
                                  mode = "undirected",
                                  weighted = TRUE,
                                  diag = FALSE)
z_umap <- cluster_louvain(G1)$membership
z_umap_s <- switch_labels(z,z_umap)
z_umap_sm <- cluster_sm(G3,D2,c_init = z)
z_umap_sm_s <- switch_labels(z,z_sm)

cluster_accuracy(z,z_umap_s)
cluster_accuracy(z,z_umap_sm_s)

ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_umap_sm))) + 
    geom_point() + 
    theme(legend.position = "none")

z_mat <- cluster_sm(G1,D2)
p1 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[1,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 1")
p2 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[2,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 2")
p3 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[3,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none")+
    ggtitle("Iteration 3")
p4 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[4,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 4")
p5 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[5,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 5")
p6 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[6,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 6")
p7 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[7,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 7")
p8 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[8,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 8")
p9 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[9,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 9")
p10 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[10,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 10")
p11 <- ggplot(UMAP, aes(x = UMAP1, y = UMAP2, color = as.factor(z_mat[11,]))) + 
    geom_point() + 
    theme_void() + 
    theme(legend.position = "none") +
    ggtitle("Iteration 11")
ggsave(p1,filename = "demonstration_plots/p1.jpg",height = 5,width = 5)
ggsave(p2,filename = "demonstration_plots/p2.jpg",height = 5,width = 5)
ggsave(p3,filename = "demonstration_plots/p3.jpg",height = 5,width = 5)
ggsave(p4,filename = "demonstration_plots/p4.jpg",height = 5,width = 5)
ggsave(p5,filename = "demonstration_plots/p5.jpg",height = 5,width = 5)
ggsave(p6,filename = "demonstration_plots/p6.jpg",height = 5,width = 5)
ggsave(p7,filename = "demonstration_plots/p7.jpg",height = 5,width = 5)
ggsave(p8,filename = "demonstration_plots/p8.jpg",height = 5,width = 5)
ggsave(p9,filename = "demonstration_plots/p9.jpg",height = 5,width = 5)
ggsave(p10,filename = "demonstration_plots/p10.jpg",height = 5,width = 5)
ggsave(p11,filename = "demonstration_plots/p11.jpg",height = 5,width = 5)

