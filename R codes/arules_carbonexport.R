library('arules')
library('gplots')
library('RColorBrewer')
#library('SDMTools')

nom <- commandArgs(trailingOnly = T)[1]
supp <- commandArgs(trailingOnly = T)[2]
launch <- commandArgs(trailingOnly = T)[3]
dt <- readRDS('training_data_carbon_export_all.rds')
coords <- readRDS('coords_carbon_export_all.rds')
fl <- dt[,1]
# dt$delta_flux <- as.factor(sign(dt$delta_flux))
for (i in names(dt)[1:dim(dt)[2]]){
  dt[[i]] <- as.factor(sign(dt[[i]]))
}
cond_all <- rep(T, length(coords[,1]))
cond_eq <- abs(coords[,1]) <20 
cond_subtr_n <- coords[,1] >20 & coords[,1]<40
cond_subtr_s <- coords[,1] < -20 & coords[,1] > -40
cond_subp <- abs(coords[,1]) >40
list_conds <- list(cond_eq, cond_subtr_n,cond_subtr_s, cond_subp)

apriori_launch <- function(dt0, cond, num, sp){
  dt0 <- dt0[cond,]
  test <- arules::apriori(dt0, parameter =list(supp=sp, conf=0.9, target="rules", minlen=2, maxlen=num),
                          control = list(filter=0.1, tree=T, heap=T, memopt=T, load=T, sort=2, verbose=T))
  
  rules <- sort(test, by ="lift")
  
  v <- as(rules, "data.frame")
  v$r1 <- sapply(v$rules, FUN = function(x){strsplit(as.character(x), '=>')[[1]][1]})
  v$r2 <- sapply(v$rules, FUN = function(x){strsplit(as.character(x), '=>')[[1]][2]})
  u <- v[grepl("delta_flux",v$r2),]
  u$length <- sapply(u$r1, FUN=function(x){length(strsplit(x, ',')[[1]])})
  t <- u[!grepl('=0', u$r1),]
  #saveRDS(t, 'rules.rds')
  tnew <- NULL
  for (j in unique(t$count)){
    h <- t[t$count==j,]
    k <- h[which.max(h$length),]
    tnew <- rbind(tnew, k)
  }
  return(tnew)
}

if (launch=='1'){
  res <- NULL
  i=1
  for (cds in list_conds){
    resu <- apriori_launch(dt0 = dt, cond = cds, num=as.numeric(nom), sp= as.numeric(supp))
    if (!is.null(resu)){
      res[[i]] <- resu
    } else{
      res[[i]] <- NA
    }
    i=i+1
  }
  names(res)<- c( 'equatorial', 'subtropical_north','subtropical_south', 'subpolar')
  saveRDS(res, paste('arules_carbon_export_',nom,'_',supp,'.rds', sep=''))
} else{
  res <- readRDS(paste('arules_carbon_export_',nom,'_',supp,'.rds', sep=''))
}


for (i in 1:length(res)){
  dt0 <- dt[list_conds[[i]],]
  flu <- fl[list_conds[[i]]]
  dot <- res[[i]]
  dat <- matrix(NA, nrow=dim(dot)[1], ncol = length(colnames(dt))+1)
  index_mat0 <- matrix(NA, nrow=dim(dot)[1], ncol = length(colnames(dt))+1)
  flux_mat <- matrix(NA, nrow=dim(dot)[1], ncol = length(colnames(dt))+1)
  colnames(dat) <- c(colnames(dt)[1], NA, colnames(dt)[2:length(colnames(dt))])
  for (j in 1:dim(dot)[1]){
    flux <- ifelse(grepl('=1', dot$r2[j]), 2, -2)
    index_mat0[j,1] <- dot$lift[j]
    dat[j, 1] <- flux
    trans <- strsplit(dot$r1[j], '\\{')[[1]][2]
    trans <- strsplit(trans, '\\}')[[1]][1]
    trans <- strsplit(trans, ',')[[1]]
    cd <- rep(T, length(flu))
    for (tr in trans){
      tx <- strsplit(tr, '=')[[1]][1]
      sign <- as.numeric(strsplit(tr, '=')[[1]][2])
      dat[j, which(colnames(dt)==tx)+1] <- sign
      cd <- cd & dt0[,which(colnames(dt0)==tx)]==sign
    }
    flux_mat[j,1] <- mean(flu[cd], na.rm=T)
  }
  n <-20
  pdf(paste('associations_carbon_export_',nom,'_', names(res)[i],'_',supp, '.pdf', sep=''))
  par(mfrow=c(1,2))
  dat[is.na(dat)] <- 0
  dat[,2] <- NA
  heatmap.2(dat, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
            dendrogram = "none", keysize=1,margins=c(10,5),
            labCol=c('sign of delta flux (lift)', NA, colnames(dt)[2:length(colnames(dt))]),las=2,
            labRow=NA, col= c('goldenrod1','mediumorchid1', 'gray90', 'royalblue', 'deepskyblue'), 
            symkey = F, cexRow=1.5,
            cellnote = ifelse(is.na(index_mat0), NA, format(round(index_mat0, 1))), 
            notecol = 'black' )
  mx <- 5
  bks <- seq(-mx, mx, 2*mx/100)
  col_scale=colorRampPalette(c('red', 'orangered', 'orange',
                               'grey92','lightblue1' ,'blue', 'darkblue'))(100)
  heatmap.2(flux_mat, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
            dendrogram = "none", keysize=1,margins=c(10,5),
            labCol = c('delta flux (gC/m2/year)', NA),
            labRow=NA, breaks = bks,
            col= col_scale, 
            symkey = F, cexRow=1.5 , key.xtickfun = function() {
              breaks = pretty(parent.frame()$breaks)
              breaks = breaks[c(1,length(breaks))]
              list(at = parent.frame()$scale01(breaks),
                   labels = breaks)
            })
  pnt <- cbind(x =c(0,17,17,0), y =c(0,50,50,0))
  xmx <- 25
 # plot(0,0, col='white', xlim=c(0,xmx), ylim=c(0, 60), axes=FALSE, frame.plot=F, xlab = '', ylab='')
 # SDMTools::legend.gradient(pnt, col_scale, limits=c(-round(mx,-10),round(mx,-10)), title ='Delta flux (gC/m2/year)' , cex=1)
  dev.off()
}

