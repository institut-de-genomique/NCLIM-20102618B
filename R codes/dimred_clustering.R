library('ClusterR')
library('parallel')
cambria <- readRDS('cambria.rds')
type = commandArgs(trailingOnly = T)[1]
mod = commandArgs(trailingOnly = T)[2]
period = commandArgs(trailingOnly = T)[3]
phate_fit <- readRDS(paste(mod,'_fit_',period,'_',type,'.rds', sep=''))

if (mod=='phate'){
  data <- phate_fit$embedding
  if (type=='discrete'){
    mx_clust =150
  } else{
    mx_clust <-15
  }

} else if (mod=='tsne'){
  data <- phate_fit$Y
  if (type=='discrete'){
    mx_clust =150
  } else{
    mx_clust <-20
  }
}

# subset <- sample(1:dim(data)[1], 100)
# tu <- Optimal_Clusters_Medoids(data = data[subset,], plot_clusters=F,seed=2,threads = 1,
#                                criterion = 'silhouette', max_clusters = 15, distance_metric = "euclidean")


tu <- Optimal_Clusters_Medoids(data = data, plot_clusters=F,seed=2,threads = 5,
                               criterion = 'silhouette', max_clusters = mx_clust, 
                               distance_metric = "euclidean")


av_sil <- NULL
av_sim <- NULL
for (j in 2:mx_clust){
  cl <- tu[[j]]
  avg_sil <- cl$avg_width_silhouette
  avg_disim <- cl$avg_intra_clust_dissimilarity
  av_sil <-append(av_sil, avg_sil)
  av_sim <- append(av_sim, avg_disim)
}

opt_n <- which.max(av_sil[4:length(av_sil)])+4
pdf(family="Helvetica",paste('cluster_silhouette_',type,'_',mod,'_',period,'.pdf', sep=''))
plot(2:mx_clust, av_sim, col='red', xlab='n_clust', ylab='Dissimilarity')
if (mod=='phate'){
  points(2:mx_clust, av_sil/100, col='blue', pch=21)
} else if (mod=='tsne'){
  points(2:mx_clust, av_sil*10, col='blue', pch=21)
}
plot(2:mx_clust, av_sil, col='blue', pch=21,xlab='n_clust', ylab='avg_width_silhouette', 
     main=paste('n_opt=',opt_n))
dev.off()
# n_opt is found to be the optimal number of clusters


