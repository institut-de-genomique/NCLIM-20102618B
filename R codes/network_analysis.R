#!bin/usr/bin/env Rscript
library("gbm")
library('randomForest')
library('mgcv')
library('nnet')
library("dismo")
library("FactoMineR")
# library("factoextra")
#library("readxl")
# library("ggplot2")
# library("reshape2")
library("gplots")
# library("plotly")
library("stringr")
# library("caret")
library('mapproj')
library('mapplots')
library('maptools')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library("CDFt")
library('plotrix')
library('png')
library('grid')
library('VennDiagram')
#library('tidygraph')
#library('ggraph')
library('igraph')
library('animation')

setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
source('axis_map0.R')
source('hide_arctic.R')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
lati<-seq(-89.5, 60, 1)
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
lati_sorted<-seq(-89.5, 60, 1)
longi_sorted =sort(longi)
cambria <- readRDS('cambria.rds')

dom_2006 <- readRDS('dominant_communities_2006.rds')
dom_2090 <- readRDS('dominant_communities_2090.rds')
bray_curtis <- readRDS('model-mean_bray_curtis.rds')
number_of_changes <- readRDS('number_of_changes.rds')
rev_mapping <- readRDS('rev_mapping_lo_lt.rds')
mapping <- readRDS('mapping_lo_lt.rds')
bc_new <- readRDS('model-mean_bray-curtis_dom.rds')
bc_new <- bc_new[rev_mapping]
bray_curtis_fish_eez <- readRDS('model-mean_bray_curtis_cond_fish_eez.rds')
bray_curtis_fish_eez$cond_fish[is.na(bray_curtis_fish_eez$cond_fish)] <- F
bray_curtis_fish_eez$cond_eez[is.na(bray_curtis_fish_eez$cond_eez)] <- F

fracs <- c("180-2000", "20-180", "5-20", '0.8-5', "0.22-3", '0-0.2' )
niches_annotation <- NULL
annotations0 <- c('polar', 'temperate','tropico-equatorial',
                 'polar', 'tropico-equatorial', 'temperate',
                 'temperate', 'equatorial',  'tropical',
                 'polar', 'cryptic', 'subtropical', 'temperate', 'equatorial', 'cryptic', 'tropical',
                 'subtropical', 'equatorial', 'tropical','temperate1', 'temperate', 'polar',
                 'tropical', 'subtropical', 'equatorial', 'cryptic', 'temperate')
color_vec = c('deepskyblue', 'forestgreen', 'orangered', 'darkred', 'red','grey','darkgreen','darkblue' )
annotations <- c('p',  'a', 'te',
                 'p',  'te', 'a',
                 'a', 'e', 't',
                 'p', 'c', 'st', 'a', 'e', 'c', 't',
                 'st', 'e', 't','a1', 'a', 'p',
                 't', 'st', 'e', 'c', 'a')
count <- 1
for (i in 1:6){
  
  fraction <- fracs[i]
  for (k in sort(unique(dom_2006[,2*i-1]))){
    niches_annotation <- rbind(niches_annotation, c(fraction, k, annotations[count], annotations0[count]))
    count<- count+1
  }
}
write.table(niches_annotation, 'niches_annotation.txt', row.names = F, col.names = F)

coffs =c(-0.01,-0.01,-0.01, 1/6, 1/6, 1/6)
names <- c('no_cut_off','no_cut_off_fish','no_cut_off_eez', 
           'cut_off_0_16', 'cut_off_0_16_fish', 'cut_off_0_16_eez')
for (az in 1:length(coffs)){
  cof <- coffs[az]
  name <- names[az]
  if (az %in% c(1,4)){
    conditional <- bc_new>cof
  } else if (az %in% c(2,5)){
    conditional <- bc_new>cof & bray_curtis_fish_eez$cond_fish==T
  } else if (az %in% c(3,6)){
    conditional <- bc_new>cof & bray_curtis_fish_eez$cond_eez==T
  }
  
  number_of_changes_sel <- number_of_changes[conditional]
  dc_06 <- dom_2006[conditional,seq(1,12,2)]
  dc_90 <- dom_2090[conditional,seq(1,12,2)]
  
  func_com <- function(u){
    com <- u[1]
    for (i in u[2:length(u)]){
      com <- paste(com, i, sep='_')
    }
    return(com)
  }
  
  com06 <- apply(dc_06,1, FUN=func_com)
  com90 <- apply(dc_90,1, FUN=func_com)
  if (az==1){
    com06_all <- unique(com06)
    com90_all <- unique(com90)
  }
  
  dc_06_annot <- NULL
  dc_90_annot <- NULL
  for (i in 1:6){
    selec <- which(niches_annotation[,1]==fracs[i])
    choices <- niches_annotation[selec,3]
    vec <- choices[match(dc_06[,i], niches_annotation[selec,2])]
    vec1 <- choices[match(dc_90[,i], niches_annotation[selec,2])]
    dc_06_annot <- cbind(dc_06_annot, vec)
    dc_90_annot <- cbind(dc_90_annot, vec1)
  }
  
  com06_an <- apply(dc_06_annot,1, FUN=func_com)
  com90_an <- apply(dc_90_annot,1, FUN=func_com)
  
  data_network <- as.numeric(com06==com90)
  for (i in which(data_network==0)){
    if ( (com90[i] %in% com06_all) & (com06[i] %in% com90_all) ){
      data_network[i]=2
    } else if (!(com90[i] %in% com06_all) & (com06[i] %in% com90_all)){
      data_network[i]=4 # 2090 specific
    } else if (!(com06[i] %in% com90_all) & (com90[i] %in% com06_all)){
      data_network[i]=3 # 2006 specific
    } else if (!(com06[i] %in% com90_all) & !(com90[i] %in% com06_all)){
      data_network[i]=5 # 2006 +2090 specific
    }
  }
  
  bray_curtis_co <- bray_curtis[conditional,]
  data_contour<- rep(list(matrix(NA, ncol=150, nrow=360)),5)
  data_contour1 <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(bray_curtis_co)[1]){
    dc=data_network[i]
    data_contour[[dc]][which(longi_sorted==bray_curtis_co$Long[i]),which(lati==bray_curtis_co$Lat[i])]= dc
    data_contour1[which(longi_sorted==bray_curtis_co$Long[i]),which(lati==bray_curtis_co$Lat[i])]= 1
  }
  data_contour1[is.na(data_contour1)]<-0
  for (i in 1:5){
    data_contour[[i]][is.na(data_contour[[i]])]<-0
  }
  color_pal = c('white','lightblue', 'darkblue', 'red', 'darkorchid1', 'blue')
  pdf(family="Helvetica",file=paste('changes_network_',name,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  set_col=color_pal
  ec=scales::alpha('white',0)
  legs = c('No significant change','Same network','Change of network','2006 specific network','2090 specific network', '2006 & 2090 specific')
  for (i in 1:5){
    .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  
                    col=c(ec,set_col[i+1], set_col[i+1], set_col[i+1]),levels=c(0,0.1,i, i+1))
  }
  # contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
  # contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
  # contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
  # legend(x = 59, y = 66, 
  #        legend = legs,
  #        col=color_pal,pch=15,cex=8, ncol=1, bty = 'n')
  # contour(x=longi_sorted,y=lati_sorted, z=data_contour,col='white',add=T, lwd=0.4,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
  pdf(family="Helvetica",file=paste('changes_network_legend.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  legend(x = 59, y = 66, 
         legend = legs,
         col=color_pal,pch=15,cex=8, ncol=1, bty = 'n')
  hide_arctic()
  axis_map0()
  dev.off()
  
  cos_lat<-readRDS('model-mean_weights_cos.rds')
  cos_lat_selec <- cos_lat[conditional]
  cos_lat <- cos_lat[mapping]
  stats_net <- NULL
  for (i in 1:5){
    selection <- data_network == i
    area_ <- sum(cos_lat_selec[selection]*111*111)/sum(cos_lat*111*111)
    stats_net <- rbind(stats_net, c(legs[i+1], area_))
  }
  stats_net <- rbind(c('No significant change', 
                       1-sum(as.numeric(stats_net[,2]), na.rm=T)),
                     stats_net)
  write.table(stats_net, paste('stats_area_network_',name,'.txt', sep=''))
  pdf(family="Helvetica",'barplot_network_',name,'.pdf', width = 8)
  barplot(as.numeric(stats_net[,2]), col=color_pal, names.arg=stats_net[,1], cex.names=0.6)
  dev.off()


  stats_net0 <- NULL
  stats_net0 <- rbind(stats_net0, c('2006', length(unique(com06))))
  stats_net0 <- rbind(stats_net0, c('2090', length(unique(com90))))
  unique_90 <- length(unique(com90)[!(unique(com90) %in% unique(com06))])
  unique_06 <- length(unique(com06)[!(unique(com06) %in% unique(com90))])
  stats_net0 <- rbind(stats_net0, c('2006 specific', unique_06))
  stats_net0 <- rbind(stats_net0, c('2090 specific', unique_90))
  stats_net0 <- rbind(stats_net0, c('common', length(intersect(com06, com90))))
  write.table(stats_net0, paste('stats_2006_vs_2090_network',name,'.txt', sep=''))
  pdf(family="Helvetica",paste('venn_diagram_network_',name,'.pdf', sep=''))
  draw.pairwise.venn(length(unique(com06)), length(unique(com90)), length(intersect(com06, com90)), 
                     col = c('red', 'darkorchid1' ), fill = c('red', 'darkorchid1' ), category=c('2006', '2090'))
  dev.off()


  data_fractions <- list()
  data_frac <- NULL
  changes_frac <- NULL
  data_fractions_an <- list()
  data_frac_an <- NULL
  changes_frac_an <- NULL
  for (i in 1:6){
    vec <- which(dom_2006[,2*i-1] == dom_2090[,2*i-1])
    vec1 <- as.numeric(dom_2006[,2*i-1] == dom_2090[,2*i-1])
    vec2 <- paste(dom_2006[,2*i-1],dom_2090[,2*i-1], sep='_')
    vec3 <- which(dc_06_annot[,i] == dc_90_annot[,i])
    vec4 <- as.numeric(dc_06_annot[,i] == dc_90_annot[,i])
    vec5 <- paste(dc_06_annot[,i],dc_90_annot[,i], sep='_')
    print(length(vec))
    data_fractions[[i]]=vec
    data_frac <- cbind(data_frac, vec1)
    changes_frac  <- cbind(changes_frac, vec2)
    data_fractions_an[[i]]=vec3
    data_frac_an <- cbind(data_frac_an, vec4)
    changes_frac_an  <- cbind(changes_frac_an, vec5)
  }
  
  percent_common <- function(u, mat){
    percent_common = apply(mat, 2, function(w) sum(u==w)/length(u))
  }
  percent_mat = apply(data_frac, 2, percent_common, data_frac)
  percent_mat <- as.matrix(percent_mat)
  vector_shift <- NULL
  vector_shift_an <- NULL
  for (i in 1:6){
    print(length(unique(changes_frac[,i]))-length(unique(dom_2006[,2*i-1])))
    percent_mat[i,i]=NA
    vector_shift <- paste(vector_shift, changes_frac[,i])
    vector_shift_an <- paste(vector_shift_an, changes_frac_an[,i])
  }
  #venn.diagram(data_fractions, filename ='test.tiff', category.names = c("180-2000", "20-180", "5-20", '0.8-5', "0.22-3" ))
  pdf(family="Helvetica",paste('common_change_fractions_network_',name,'.pdf', sep=''))
  heatmap.2(percent_mat, trace="none",symm=TRUE, 
            labRow = c("180-2000", "20-180", "5-20", '0.8-5', "0.22-3", '0-0.2' ),
            labCol = c("180-2000", "20-180", "5-20", '0.8-5', "0.22-3", '0-0.2' ),
            keysize=1,margins=c(10,10), col= greenred(10), 
            symkey = F, hclustfun = hclust)
  dev.off()
  
  data_frac_changes_only <- data_frac[number_of_changes!=0 & conditional,]
  changes_frac_changes_only <- changes_frac[number_of_changes!=0 & conditional,]
  vector_shift_changes_only <- vector_shift[number_of_changes!=0 & conditional]
  
  data_frac_changes_only_an <- data_frac_an[number_of_changes_sel!=0,]
  changes_frac_changes_only_an <- changes_frac_an[number_of_changes_sel!=0,]
  vector_shift_changes_only_an <- vector_shift_an[number_of_changes_sel!=0]
  
  test <- as.vector(table(vector_shift_changes_only))
  test_an <- as.vector(table(vector_shift_changes_only_an))
  top_changes <- vector_shift_changes_only[order(test, decreasing = T)][1:10]
  unique_changes_frac <- list()
  unique_changes_frac_an <- list()
  glob_stats_changes <- NULL
  for (i in 1:6){
    unique_changes_frac[[i]]<- data.frame(sort(table(changes_frac_changes_only[data_frac_changes_only[,i]==0,i]), decreasing = T))
    unique_changes_frac_an[[i]]<- data.frame(sort(table(changes_frac_changes_only_an[data_frac_changes_only[,i]==0,i]), decreasing = T))
    vec <- c(fracs[i])
    for (j in 1:2){
      first <- strsplit(as.character(unique_changes_frac_an[[i]]$Var1[j]), '_')[[1]][1]
      second <- strsplit(as.character(unique_changes_frac_an[[i]]$Var1[j]), '_')[[1]][2]
      c1 <- annotations0[match(first, annotations)]
      c2 <- annotations0[match(second, annotations)]
      vec<-append(vec,c(as.character(unique_changes_frac[[i]]$Var1[j]),
                        as.character(unique_changes_frac_an[[i]]$Var1[j]),
                        paste(c1,c2,sep='->'),
                        as.integer(unique_changes_frac[[i]]$Freq[j])/sum(unique_changes_frac[[i]]$Freq)))
    }
    glob_stats_changes <- rbind(glob_stats_changes, vec)
  }
  write.table(glob_stats_changes, paste('glob_stats_changes_frac_network_',name,'.txt', sep=''), 
              col.names = F, row.names = F)
  
  coms06_annot_distrib <- data.frame(table(com06_an))
  coms90_annot_distrib <- data.frame(table(com90_an))
  coms06_vs_90 <- data.frame('06'=com06_an, '90'=com90_an)
  number_of_changes_sel <- number_of_changes[conditional]
  coms06_vs_90 <- coms06_vs_90[number_of_changes_sel!=0,]
  coms06_vs_90 <- cbind(coms06_vs_90, cos_lat_selec[as.numeric(row.names(coms06_vs_90))]*111*111)
  connection <- paste(coms06_vs_90$X06, coms06_vs_90$X90)
  coms06_vs_90 <- cbind(coms06_vs_90,connection)
  areas_vec <- NULL
  for (co in unique(connection)){
    area <- sum(coms06_vs_90[coms06_vs_90$connection==co,3])
    areas_vec <- append(areas_vec, area)
  }
  #graph_0690 <- data.frame(table(connection))
  graph_0690 <- data.frame(unique(connection) ,areas_vec)
  graph_0690[,1]<- as.character(graph_0690[,1])
  nodes06_90 <- t(as.data.frame(strsplit(graph_0690[,1], ' ')))
  nodes06_90 <- cbind(nodes06_90, graph_0690[,2])
  row.names(nodes06_90)<-NULL
  nodes06_90 <- as.data.frame(nodes06_90)
  nodes06_90$V3 <- as.numeric(as.character(nodes06_90$V3))
  nodes06_90 <- nodes06_90[!is.na(nodes06_90$V3),]
  nodes06_90$V4 <- nodes06_90$V3/sum(nodes06_90$V3, na.rm = T)
  cumul_area <- 0
  vec_cum_area <-NULL
  for (u in sort(nodes06_90$V4, decreasing = T)){
    cumul_area = cumul_area + u
    vec_cum_area <- append(vec_cum_area, cumul_area)
  }
  selecs <- order(nodes06_90$V3, decreasing = T)
  nodes06_90<-nodes06_90[selecs,]
  nodes06_90$V5 <- vec_cum_area
  numbs <- number_of_changes_sel[number_of_changes_sel!=0][match(graph_0690[,1], connection)]
  numbs <- numbs[selecs]
  #quartiles <- quantile(nodes06_90$V4, probs = seq(0,1,0.25))
  cut_offs <- c(0.25, 0.5)
  nodes06_90_2 <- nodes06_90[nodes06_90$V5<0.75,]
  vi<-log10(nodes06_90_2$V3)
  for (co in cut_offs){
    nodes06_90_1 <- nodes06_90[nodes06_90$V5<co,]
    numbs_1 <- numbs[nodes06_90$V5<co]
    
    net <- graph_from_data_frame(nodes06_90_1)
    annot_nodes <- NULL
    col_code<-c('red', 'darkorchid1' , 'darkblue')
    col_nodes <- NULL
    area_2006 <- NULL
    area_2090 <- NULL
    for (u in V(net)$name){
      area06 <- sum(cos_lat_selec[com06_an==u]*111*111)
      area90 <- sum(cos_lat_selec[com90_an==u]*111*111)
      if (area06==0){
        area_2006 <- append(area_2006, 0)
      } else{
        area_2006 <- append(area_2006, log10(area06))
      }
      if (area90==0){
        area_2090 <- append(area_2090, 0)
      } else{
        area_2090 <- append(area_2090, log10(area90))
      }
      if (u %in% unique(com06_an) & u %in% unique(com90_an)){
        annot_nodes <- append(annot_nodes, 'common_network')
        col_nodes <- append(col_nodes, col_code[3])
      } else if ( !(u %in% unique(com06_an)) & u %in% unique(com90_an)){
        annot_nodes <- append(annot_nodes, 'specific_network_2090')
        col_nodes <- append(col_nodes, col_code[2])
      } else if ( u %in% unique(com06_an) & !(u %in% unique(com90_an))){
        annot_nodes <- append(annot_nodes, 'specific_network_2006')
        col_nodes <- append(col_nodes, col_code[1])
      }
    }
    deg <- degree(net)
    E(net)$log10_number<-log10(nodes06_90_1$V3)
    E(net)$number <- numbs_1
    V(net)$Network_type <- annot_nodes
    #V(net)$nodes_colour <- col_nodes
    #75 h=7, w=10
    #50 h=11, w=16
    #25 h=20, w=30
    #0 h=60, w=90
    color_frame <-cbind(color_vec, unique(annotations0), unique(annotations))
    
    values_color <- NULL
    for (v in V(net)$name){
      u =strsplit(v, '_')
      for (com in u[[1]]){
        values_color<-append(values_color,color_frame[color_frame[,3]==com,1])
      }
    }
    values_color <- matrix(values_color, ncol=6, byrow = T)
    col_list <- color_vec[match(annotations0, color_frame[,2])]
    selec<-readRDS('fractions0.rds')
    dec <- c(0,3,6,9,16,21)
    values<-apply(values_color,1,function(x){
      vec = rep(0,length(col_list))
      i=1
      for (coli in x){
        vec[match(coli, col_list[selec[[i]]])+dec[i]]=vec[match(coli, col_list[selec[[i]]])+dec[i]]+1
        i=i+1
      }
      return(vec)
      #sapply(col_list,function(y){print(y);sum(x==y)})
    })
    values<-as.list(as.data.frame(values))
    set.seed(1)
    l<-layout_nicely(net)
    pnt =cbind(x =c(1,1.05,1.05,1), y =c(1,1,1.2,1.2))
    color_edges <- colorRampPalette(c('blue', 'red'))(100)
    u=log10(nodes06_90_1$V3)
    vec <- round((u-min(vi))*100/(max(u)-min(vi)))
    pdf(family="Helvetica",paste('major_changes_all_frac_pie_', co,'_2006_network_',name,'.pdf', sep=''), width = 10, height = 10)
    plot.igraph(net, vertex.shape="pie", vertex.pie=values,
                vertex.pie.color=list(col_list),
                vertex.size=area_2006**3/30, layout=l, edge.arrow.size=0.9, 
                vertex.label=rep(NA,length(deg)), edge.width=numbs_1*2, 
                edge.color=color_edges[vec])
    
    legend.gradient(pnt, cols=color_edges, limits=c(round(min(nodes06_90_2$V3)), round(max(nodes06_90_1$V3))))
    legend(x = -1.2, y=-0.9, legend = unique(annotations0), fill = unique(color_vec), bty='n')
    text(x=1.1,y=-1.5+max(numbs_1)/10+0.1, 'Number of changes')
    for (i in sort(unique(numbs_1))){
      arrows(x0=1,y0=-1.5+i/10,x1=1.1,y1=-1.5+i/10, lwd=i*1.5)
      text(x=1.2, y=-1.5+i/10, i)
    }
    dev.off()
    pdf(family="Helvetica",paste('major_changes_all_frac_pie_', co,'_2090_network',name,'.pdf', sep=''), width = 10, height = 10)
    plot.igraph(net, vertex.shape="pie", vertex.pie=values,
                vertex.pie.color=list(col_list),
                vertex.size=area_2090**3/30, layout=l, edge.arrow.size=0.9, 
                vertex.label=rep(NA,length(deg)), edge.width=numbs_1*2, 
                edge.color=color_edges[vec])
    legend.gradient(pnt, cols=color_edges, limits=c(round(min(nodes06_90_2$V3)), round(max(nodes06_90_1$V3))))
    legend(x = -1.2, y=-0.9, legend = unique(annotations0), fill = unique(color_vec), bty='n')
    text(x=1.1,y=-1.5+max(numbs_1)/10+0.1, 'Number of changes')
    for (i in sort(unique(numbs_1))){
      arrows(x0=1,y0=-1.5+i/10,x1=1.1,y1=-1.5+i/10, lwd=i*1.5)
      text(x=1.2, y=-1.5+i/10, i)
    }
    dev.off()
    
#     saveGIF(
#       { plot.igraph(net, vertex.shape="pie", vertex.pie=values,
#                     vertex.pie.color=list(col_list),
#                     vertex.size=area_2006**3/30, layout=l, edge.arrow.size=0.9,
#                     vertex.label=rep(NA,length(deg)), edge.width=0,
#                     edge.color='white')
#         legend.gradient(pnt, cols=color_edges, limits=c(round(min(nodes06_90_2$V3)), round(max(nodes06_90_1$V3))))
#         legend(x = -1.2, y=-0.9, legend = unique(annotations0), fill = unique(color_vec), bty='n')
#         text(x=1.1,y=-1.5+max(numbs_1)/10+0.1, 'Number of changes')
#         for (i in sort(unique(numbs_1))){
#           arrows(x0=1,y0=-1.5+i/10,x1=1.1,y1=-1.5+i/10, lwd=i*1.5)
#           text(x=1.2, y=-1.5+i/10, i)
#         }
#         title('2006')
#   
#         plot.igraph(net, vertex.shape="pie", vertex.pie=values,
#                     vertex.pie.color=list(col_list),
#                     vertex.size=area_2006**3/30, layout=l, edge.arrow.size=0.9,
#                     vertex.label=rep(NA,length(deg)), edge.width=numbs_1*2,
#                     edge.color=color_edges[vec])
#         legend.gradient(pnt, cols=color_edges, limits=c(round(min(nodes06_90_2$V3)), round(max(nodes06_90_1$V3))))
#         legend(x = -1.2, y=-0.9, legend = unique(annotations0), fill = unique(color_vec), bty='n')
#         text(x=1.1,y=-1.5+max(numbs_1)/10+0.1, 'Number of changes')
#         for (i in sort(unique(numbs_1))){
#           arrows(x0=1,y0=-1.5+i/10,x1=1.1,y1=-1.5+i/10, lwd=i*1.5)
#           text(x=1.2, y=-1.5+i/10, i)
#         }
#         title('2006')
#   
#         plot.igraph(net, vertex.shape="pie", vertex.pie=values,
#                     vertex.pie.color=list(col_list),
#                     vertex.size=area_2090**3/30, layout=l, edge.arrow.size=0.9,
#                     vertex.label=rep(NA,length(deg)), edge.width=numbs_1*2,
#                     edge.color=color_edges[vec])
#         legend.gradient(pnt, cols=color_edges, limits=c(round(min(nodes06_90_2$V3)), round(max(nodes06_90_1$V3))))
#         legend(x = -1.2, y=-0.9, legend = unique(annotations0), fill = unique(color_vec), bty='n')
#         text(x=1.1,y=-1.5+max(numbs_1)/10+0.1, 'Number of changes')
#         for (i in sort(unique(numbs_1))){
#           arrows(x0=1,y0=-1.5+i/10,x1=1.1,y1=-1.5+i/10, lwd=i*1.5)
#           text(x=1.2, y=-1.5+i/10, i)
#         }
#         title('2090')
#         plot.igraph(net, vertex.shape="pie", vertex.pie=values,
#                     vertex.pie.color=list(col_list),
#                     vertex.size=area_2090**3/30, layout=l, edge.arrow.size=0.9,
#                     vertex.label=rep(NA,length(deg)), edge.width=0,
#                     edge.color='white')
#         legend.gradient(pnt, cols=color_edges, limits=c(round(min(nodes06_90_2$V3)), round(max(nodes06_90_1$V3))))
#         legend(x = -1.2, y=-0.9, legend = unique(annotations0), fill = unique(color_vec), bty='n')
#         text(x=1.1,y=-1.5+max(numbs_1)/10+0.1, 'Number of changes')
#         for (i in sort(unique(numbs_1))){
#           arrows(x0=1,y0=-1.5+i/10,x1=1.1,y1=-1.5+i/10, lwd=i*1.5)
#           text(x=1.2, y=-1.5+i/10, i)
#         }
#         title('2090')
#         }, interval = 3, movie.name=paste('major_changes_all_frac_pie_', co,'_network_',name,'.gif', sep='')
#     )
  }
}
# ggraph(lay, layout = 'circlepack', weight="size")+
#   geom_node_point(colour=rep('green',length(V(net))), alpha=1, size=1)+
#   geom_node_point(colour=rep('red',length(V(net))), size=2, alpha=0.5)+
#   geom_node_point(colour=rep('blue',length(V(net))), size=3, alpha=0.25)


# pdf(family="Helvetica",paste('major_changes_all_frac_', co,'.pdf', sep=''), height = 7, width = 10)
# set.seed(3)
# lay <- create_layout(net, layout = "nicely")
# ggraph(lay) +
#   geom_edge_fan(aes(colour=log10_number, width=factor(number)), alpha = 1, 
#                 arrow = arrow(length = unit(2, "mm")),
#                 start_cap = circle(7, "mm"),
#                 end_cap = circle(7, "mm"),  
#                 angle_calc = 'along',
#                 label_dodge = unit(2.5, 'mm')) +
#   geom_node_point(size = deg/1.5)+
#   geom_node_label(aes(label=name, color=Network_type), size=2)+
#   #geom_node_text(aes(label=name), size=2.5, vjust = -1.5)+
#   theme_graph(base_family = "sans")+
#   theme(legend.margin = margin(1,1))+
#   scale_edge_color_continuous(name='log10(number of changes)',low='blue',high = 'red', limits=c(1, max(log10(nodes06_90_1$V3))))+
#   scale_edge_width_manual(name='Number of\ncommunity change',values=c(0.5,1,1.5,2, 2.5, 3))+
#   scale_color_manual(name="Type of Network",
#                      values=c(common_network='darkblue', specific_network_2090='darkorchid1', specific_network_2006='red' ),
#                      labels=c('Common network', '2006 specific', '2090 specific'))
# dev.off()
