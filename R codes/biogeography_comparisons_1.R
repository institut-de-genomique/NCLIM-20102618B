library('mclust')
library('RColorBrewer')
library('gplots')
library('tidyr')
library('dplyr')
library('pvclust')
library('dendextend')
library('vegan')
#library('readxl')
library('grDevices')

cambria <- readRDS('cambria.rds')
grd_bio <- readRDS('biogeographies.rds')
dom_woa <- readRDS('dominant_communities_woa.rds')
dom_woa <- dom_woa[,seq(1,12,2)]
clusts_phate_medoids <- readRDS('phate_cluster_medoids_prob_7_woa.rds')
clusts_phate_medoids0 <- readRDS('phate_cluster_medoids_prob_4_woa.rds')

clusts_phate_0690 <- readRDS('phate_cluster_medoids_prob_13_0690.rds')
clusts_tsne_0690 <- readRDS('tsne_cluster_medoids_prob_10_0690.rds')

seq1 <- seq(1, length(clusts_tsne_0690),2)
seq2 <-seq(2, length(clusts_tsne_0690),2)
phate_06 <- clusts_phate_0690[seq1]
phate_90 <- clusts_phate_0690[seq2]
tsne_06 <- clusts_tsne_0690[seq1]
tsne_90 <- clusts_tsne_0690[seq2]


bray_curtis <- readRDS('model-mean_bray_curtis.rds')
mapping <- match(bray_curtis$cell,grd_bio$cell)
grd_bio <- grd_bio[mapping,]
# to_remove <- !is.na(grd_bio$biome)
# grd_bio <- grd_bio[to_remove,]
# bray_curtis <- bray_curtis[to_remove,]
# dom_woa<- dom_woa[to_remove,]
fractions <- c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
colnames(dom_woa)<-fractions

to_compare <- cbind(dom_woa,clusts_phate_medoids0, clusts_phate_medoids,grd_bio[,c(3,5,6,7)], 
                    phate_06, phate_90,tsne_06, tsne_90)
colnames(to_compare)[7]<-'all_phate_4'
colnames(to_compare)[8]<-'all_phate_7'

partition_test <- function(x){
  u <- apply(to_compare,2,adjustedRandIndex, y=to_compare[,x])
}

index_mat <- t(sapply(1:16, partition_test))
rownames(index_mat)<-colnames(to_compare)

pdf(family="Helvetica",'randIndex_biogeographies.pdf', height= 7, width=10)
br <- c(seq(0, 1, length.out = 101))
set = colorRampPalette(c('blue','white', 'red'))(100)
index_mat0 <- index_mat[1:12, 1:12]
index_mat0[upper.tri(index_mat0, diag = T)] <- NA
heatmap.2(index_mat0, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1, col=set ,margins = c(13.5,13.5),
          breaks=br, symkey = F, cexRow=1, 
          cellnote = ifelse(is.na(index_mat0), NA, format(round(index_mat0, 2), nsmall = 2)),
                            notecol = 'black', denscol = 'black')
dev.off()
pdf(family="Helvetica",'randIndex_biogeographies_1.pdf', height= 7, width=10)
br <- c(seq(0, 1, length.out = 101))
set = colorRampPalette(c('blue','white', 'red'))(100)
index_mat0 <- index_mat[1:11, 1:11]
index_mat0[upper.tri(index_mat0, diag = T)] <- NA
heatmap.2(index_mat0, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1, col=set ,margins = c(13.5,13.5),
          breaks=br, symkey = F, cexRow=1, 
          cellnote = ifelse(is.na(index_mat0), NA, format(round(index_mat0, 2), nsmall = 2)),
          notecol = 'black', denscol = 'black')
dev.off()

pdf(family="Helvetica",'randIndex_biogeographies_all.pdf', height= 7, width=10)
br <- c(seq(0, 1, length.out = 101))
set = colorRampPalette(c('blue','white', 'red'))(100)
index_mat[upper.tri(index_mat, diag = T)] <- NA
heatmap.2(index_mat, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1, col=set ,margins = c(13.5,13.5),
          breaks=br, symkey = F, cexRow=1, 
          cellnote = ifelse(is.na(index_mat), NA, format(round(index_mat, 2), nsmall = 2)),
          notecol = 'black', denscol = 'black')
dev.off()
