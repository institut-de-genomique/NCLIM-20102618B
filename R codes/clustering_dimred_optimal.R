library('ClusterR')

type = commandArgs(trailingOnly = T)[1]
mod = commandArgs(trailingOnly = T)[2]
opt_n <- as.integer(commandArgs(trailingOnly = T)[3])
year <- commandArgs(trailingOnly = T)[4]
phate_fit <- readRDS(paste(mod,'_fit_',year,'_',type,'.rds', sep=''))
if (mod=='phate'){
  data <- phate_fit$embedding
} else if (mod=='tsne'){
  data <- phate_fit$Y
}
                    
clust_res <- Cluster_Medoids(data = data, 
                             clusters = opt_n, 
                             distance_metric = 'euclidean',
                             seed = 2)
#saveRDS(clust_res, file = paste(mod,'_cluster_medoids_result_',type,'_',year,'.rds', sep=''))
saveRDS(clust_res$clusters, paste(mod,'_cluster_medoids_',type,'_',opt_n,'_',year,'.rds', sep=''))
