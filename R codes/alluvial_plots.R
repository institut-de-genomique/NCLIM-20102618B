library('ggalluvial')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
dom_com06 <- readRDS('dominant_communities_2006.rds')
dom_com90 <- readRDS('dominant_communities_2090.rds')
cambria <- readRDS('cambria.rds')
coslat <- readRDS('model-mean_weights_cos.rds')
col_sets <- readRDS('colors_provinces.rds')
fractions <- c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
code_fraction <-c('F', 'E','D','C','B','A')
plot_list <- rep(list(NA), 6)
for (i in 1:6){
  pairs <- NULL
  letter <- code_fraction[i]
  uni <- unique(dom_com06[,i*2-1])
  uni <-append(uni, 0)
  for (u in uni){
    for (j in uni){
        pairs <- append(pairs, paste(u,j, sep='_'))
    }
  }
  data_alluvial <- rep(list(0), length(pairs))
  for (j in 1:length(dom_com06[,1])){
    dc06 <- dom_com06[j,i*2-1]
    dc90 <- dom_com90[j,i*2-1]
    pair = paste(dc06,dc90, sep='_')
    diff_proba = dom_com06[j,i*2]- dom_com90[j,i*2]
    if (diff_proba>0){
      pair1 = paste(dc06, 0, sep='_')
    } else{
      pair1 = paste(0, dc90, sep='_')
    }
    pair2 = paste(0,0, sep='_')
    data_alluvial[[match(pair, pairs)]] = data_alluvial[[match(pair, pairs)]] + 
      abs((dom_com06[j,i*2]-abs(diff_proba)))*111*111*coslat[j]
    data_alluvial[[match(pair1, pairs)]] = data_alluvial[[match(pair1, pairs)]] +
      abs(diff_proba)*111*111*coslat[j]
    data_alluvial[[match(pair2, pairs)]] = data_alluvial[[match(pair2, pairs)]] +
      (1-abs(diff_proba)-abs((dom_com06[j,i*2]-abs(diff_proba))))*111*111*coslat[j]
  }
  data_alluvial0 <- NULL
  for (j in 1:length(pairs)){
    pair = pairs[[j]]
    dc06 <- strsplit(pair, split = '_')[[1]][1]
    dc90 <- strsplit(pair, split = '_')[[1]][2]
    n = round(data_alluvial[[j]][[1]])
    data_alluvial0 <- rbind(data_alluvial0, c(dc06, '2006',pair, n))
    data_alluvial0 <- rbind(data_alluvial0, c(dc90, '2090',pair, n))
  }
  #data_alluvial0 <- as.data.frame(data_alluvial0)
  #cols <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(max(dom_com06[,i*2-1]))[sort(unique(dom_com06[,i*2-1]))]
  cols <- col_sets[[i]][!is.na(col_sets[[i]])]
  #print(cols)
  # cols <- rev(cols)
  cols <- append('white', cols)
  #cols <- match(unique(dom_com06[,i*2-1]), unique(data_alluvial0$dc))
  colnames(data_alluvial0) <- c('dc', 'year','pair', 'Freq')
  data_alluvial0 <-as.data.frame(data_alluvial0)
  data_alluvial0$Freq <- as.numeric(levels(data_alluvial0$Freq))[data_alluvial0$Freq]
  data_alluvial0$dc <- factor(data_alluvial0$dc, 
                              sort(as.numeric(levels(data_alluvial0$dc)), decreasing = F))
  data_alluvial0$dc0 <- paste(letter , data_alluvial0$dc, sep='')
  #data_alluvial0$dc<- as.numeric(levels(data_alluvial0$dc))[data_alluvial0$dc]
  #ploti <- 
  pdf(family="Helvetica",paste('alluvial_plot_areas_', fractions[i],'.pdf', sep=''))
  ploti <- ggplot(data_alluvial0,
         aes(x = year, stratum = dc, alluvium = pair,
             y = Freq,
             fill = dc, label = dc)) +
    scale_x_discrete(expand = c(.2, .2)) +
    geom_flow() +
    geom_stratum(alpha = 0.9) +
    scale_fill_manual(values = cols)+
    #geom_text(stat = "stratum", size = 3) +
    theme(#axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  print(ploti)
  dev.off()
  #plot_list[[i]] <- ploti
}

# for (i in 1:6){
#   pdf(family="Helvetica",paste('alluvial_plot_areas_', fractions[i],'.pdf', sep=''))
#   print(plot_list[[i]])
#   dev.off()
# }
