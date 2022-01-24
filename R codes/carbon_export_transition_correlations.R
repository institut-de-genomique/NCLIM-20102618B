library('FactoMineR')
library('RColorBrewer')
library('gplots')

trans <- readRDS('transitions.rds')
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
longi_sorted =sort(longi)
lati<-seq(-89.5, 60, 1)
for (i in 1:6){
  trans[[i]][trans[[i]]!='FALSE'] <- sapply(trans[[i]][trans[[i]]!='FALSE'], 
                                  FUN = function(x){strsplit(x, ' ')[[1]][1]})
}
frcs <- c('0-0.2', '0.22-3','0.8-5', '5-20', '20-180', '180-2000')
names(trans)<- rev(frcs)
sigs <- readRDS('significant_transitions_wilcox_sizes_bis.rds')
func <- function(x){
  v <-strsplit(x, ' ')[[1]]
  u <- paste(v[1:(length(v)-1)], collapse = ' ')
  return(u)
}
sigs$type0 <- sapply(sigs$type, FUN = func)
sigs$fulltype <- paste(sigs$type0, sigs$frac)
type <- '_extrapolated_henson'
delta_flux <- readRDS(paste('delta_fluxes_all',type,'.rds', sep=''))

rev_mapping <- readRDS('rev_mapping_lo_lt.rds')
bc_new <- readRDS('model-mean_bray-curtis_dom.rds')
bc_new <- bc_new[rev_mapping]
cond_bc <- bc_new>1/6
dat <- NULL
for (j in unique(sigs$fulltype)){
  fr <- tail(strsplit(j, ' ')[[1]], n=1)
  v <- trans[[fr]]
  for (tr in unique(v)){
    if (tr %in% sigs$Trans[sigs$fulltype==j] ){
      v[v==tr]<- sigs$value[sigs$frac==fr & sigs$Trans==tr & sigs$fulltype==j]
    } else{
      v[v==tr]<-0
    }
  }
  dat <- cbind(dat, as.numeric(v))
}
colnames(dat)<-unique(sigs$fulltype)
dat <- as.data.frame(dat)
dat <- cbind(delta_flux,dat)
cond <- !is.na(dat$delta_flux) & cond_bc
coords <- readRDS('new_stations_1deg.rds')
coords <- coords[cond,]
dat<- dat[cond,] 



qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

dat <- dat[, c(1,order(colnames(dat)[2:dim(dat)[2]])+1)]
sel <- !is.na(dat$delta_flux)
sel0 <- apply(dat, 1, FUN = function(x){sum(x, na.rm=T)})
sel1 <- apply(dat, 2, FUN = function(x){sum(x, na.rm=T)})
dat <- dat[sel & sel0!=0  ,]
coords <- coords[sel & sel0!=0,]
dat <- dat[, sel1 !=0 & !is.nan(sel1) & !is.infinite(sel1)]
u <- PCA(dat, graph = F)
res.km <- kmeans(u$ind$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
u$ind$coord[,1] <- (u$ind$coord[,1]-min(u$ind$coord[,1]))*2/(max(u$ind$coord[,1])-min(u$ind$coord[,1]))-1
u$ind$coord[,2] <- (u$ind$coord[,2]-min(u$ind$coord[,2]))*2/(max(u$ind$coord[,2])-min(u$ind$coord[,2]))-1
acronyms <- c(paste('Delta', 'Export', sep=' '), 'AlMicro 20-180','AlNano 0.22-3', 'AlNano 0.8-5', 'AlNano 180-2000', 'AlNano 20-180', 'AlNano 5-20',
              'AlPico 0.22-3', 'AlPico 0.8-5', 'CyaMicro 20-180', 'CyaMicro 5-20', 'CyaNano 0.22-3', 'CyaNano 0.8-5', 'CyaNano 180-2000', 'CyaNano 20-180', 'CyaNano 5-20',
              'CyaPico 0.22-3','CyaPico 0.8-5','CyaPico 5-20','DiaNano 0.8-5', 
              'DiaNano 20-180', 'DiaNano 5-20', 'Diaz 0.8-5', 'Diaz 5-20', 'MesoZoo 180-2000', 'MicroZoo 180-2000', 'MicroZoo 20-180')
acronyms <- acronyms[sel1 !=0]
ku <- length(colnames(dat))
pdf('PCA_carbon_export.pdf')
plot.PCA(u, axes=c(1, 2), choix="ind", habillage="none", col.ind=scales::alpha('white', 0), 
         ylim =c(1.15*min(u$var$coord[,2]), 1.15*max(u$var$coord[,2])), 
         xlim= c(1.15*min(u$var$coord[,1]), 1.15*max(u$var$coord[,1])))#, col.ind = cols[grp]
x <- u$var$coord[,1]
y<-u$var$coord[,2]
selec <- abs(x)>0.1 | abs(y)>0.1
arrows(x0=rep(0, ku),y0= rep(0, ku),x1=x, y1=y, col =col_vector[1:ku], lwd = 2 )
text(x[selec]+0.04*sign(x[selec]) , y[selec]+0.04*sign(y[selec]),cex=0.5,
     labels=acronyms[selec])
plot.PCA(u, axes=c(3, 4), choix="ind", habillage="none", col.ind=scales::alpha('white', 0), 
         ylim =c(1.15*min(u$var$coord[,4]), 1.15*max(u$var$coord[,4])), 
         xlim= c(1.15*min(u$var$coord[,3]), 1.15*max(u$var$coord[,3])))#, col.ind = cols[grp]
x <- u$var$coord[,3]
y<-u$var$coord[,4]
selec <- abs(x)>0.1 | abs(y)>0.1
arrows(x0=rep(0, ku),y0= rep(0, ku),x1=x, y1=y, col =col_vector[1:ku], lwd = 2 )
text(x[selec]+0.04*sign(x[selec]) , y[selec]+0.04*sign(y[selec]),cex=0.5,
     labels=acronyms[selec])
plot(0,0, col='white',axes=FALSE, frame.plot=F, xlab = '', ylab='')
legend('topright', bty='n', legend =  acronyms, fill = col_vector[1:ku], cex=0.85 ,ncol = 1)
dev.off()

library('randomForest')
colnames(dat)[2:dim(dat)[2]]<-acronyms[2:length(acronyms)]
saveRDS(dat, 'training_data_carbon_export_all.rds')
saveRDS(coords, 'coords_carbon_export_all.rds')
#set.seed(42)
#rf <- randomForest(x = dat[, 2:dim(dat)[2]],y=dat$delta_flux, ntree = 1000, mtry = 10, importance = T)
#imps <- rf$importance[,1]*100/sum(rf$importance[,1])
#names(imps)<-acronyms[2:length(acronyms)]
#od <- order(imps, decreasing = T)
#imps <- imps[od]
#pdf('random_forest_caron_export_variable_importance.pdf')
#par(mar= c(10.1, 5.1, 4.1, 2.1))
#barplot(imps, las=2, col = col_vector[1:ku][od])
#par(mfrow=c(2,2), mar= c(4.1, 5.1, 4.1, 2.1))
#dt <- NULL
#for (v in colnames(dat)[2:dim(dat)[2]][od]){
#  u <- randomForest::partialPlot(x = rf, pred.data = dat[, 2:dim(dat)[2]], x.var = c(v), plot = F, xlab = v, ylab='Delta flux (gC/m2/year)', 
#                            main=paste('Partial Dependence on', v), n.pt = 10)
#  dt <- rbind(dt, u$y)
#}
#dev.off()

#heatmap.2(dt, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
#          dendrogram = "none", keysize=1,margins=c(10,22), 
#          labRow=colnames(dat)[2:dim(dat)[2]][od], col= redgreen(15),  symkey = F, cexRow=1.2)




