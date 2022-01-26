source('vioplot1.R')

data_prov <- readRDS('Provinces_env_parameters_woa_scaled.rds')
data_prov$Fraction[data_prov$Fraction=='43952'] <- '5-20'
data_prov$lab <- paste(data_prov$Fraction, data_prov$Province, sep='_')
to_remove <- readRDS('excluded_niches.rds')
data_prov <- data_prov[!(data_prov$lab %in% to_remove),]
valids <- readRDS('valids.rds')
data_prov <- data_prov[data_prov$lab %in% valids,]
stats_prov <- unique(data_prov$Station)
frcs <- unique(data_prov$Fraction)
letters<- c('A', 'C', 'E', 'F', 'B', 'D')
data_prov$prov <- paste(letters[match(data_prov$Fraction, frcs)], data_prov$Province, sep='')
data_prov$stat <- as.numeric(sapply(data_prov$Station, FUN = function(x){strsplit(x, '_')[[1]][1]}))
col_prov <- read.table('color_provinces.txt', header = T)

hydro <- read.table('hydrofull_for_margaux.txt', header = T, sep='\t')

nas <- c('POC','PIC','flux500m', 'flux200m' )
for (n in nas){
  data_prov[[n]] <- hydro[[n]][match(data_prov$stat, hydro$Station)]
}

init_vec <- function(){
  comps_ca <- list(rep(list(list()), 6))
  names(comps_ca)<- 'carb'
  names(comps_ca$carb) <- frcs
  for (frac in frcs){
    for (p in sort(unique(data_prov$Province[data_prov$Fraction==frac]))){
      pr <- paste(letters[match(frac, frcs)], p,sep='')
      comps_ca[['carb']][[frac]][[pr]] <- rep(NA, length(nas))
    }
  }
  return(comps_ca)
}
comps_ca <- init_vec()
comps_cam <- init_vec()

save_dt <- rep(list(NULL), length(nas))
sig_carb <- NULL
pdf('carbon_export_violin.pdf', height = 12, width = 12)
for (frac in frcs){
  dt <- data_prov[data_prov$Fraction==frac,]
  par(mfrow=c(2,2))
  co=1
  for (n in nas){
    dt0 <- dt[!is.na(dt[[n]]),]
    an <- aov(dt0[[n]]~dt0$prov)
    #print(summary(an))
    u<-pairwise.wilcox.test(dt0[[n]], dt0$prov, p.adjust.method = 'BH')$p.value
    dta <- list()
    for (p in sort(unique(dt0$Province))){
      pr <- paste(letters[match(frac, frcs)], p,sep='')
      print(pr)
      save_dt[[co]] <- append(save_dt[[co]], median(dt0[[n]][dt0$prov==pr]))
      comps_ca[['carb']][[frac]][[pr]][co] <- median(dt0[[n]][dt0$prov==pr])
      comps_cam[['carb']][[frac]][[pr]][co] <- mean(dt0[[n]][dt0$prov==pr])
      print(median(dt0[[n]][dt0$prov==pr]))
      dta[[pr]] <- dt0[[n]][dt0$prov==pr]
    }
    mx <- max(unlist(dta))
    xs0 <- NULL
    xs1 <- NULL
    ys0 <- NULL
    ys1 <-NULL
    sigs <- NULL
    x_sigs <- NULL
    c=1
    for (i in 1:dim(u)[1]){
      for (j in 1:dim(u)[1]){
        pv <- u[i,j]
        if (!is.na(pv)){
          if (pv<0.05){
            xu <- which(names(dta)==rownames(u)[i])
            xi <- which(names(dta)==colnames(u)[j])
            xs0 <- c(xs0,xu)
            xs1 <- c(xs1,xi)
            ys0 <- c(ys0, mx+0.05*c*mx)
            ys1 <- c(ys1, mx+0.05*c*mx)
            sig_carb <- rbind(sig_carb , c(paste(rownames(u)[i], colnames(u)[j], sep='->'), pv, n))
            sig_carb <- rbind(sig_carb , c(paste(colnames(u)[j], rownames(u)[i] , sep='->'), pv, n))
            if (pv<0.0001){
              sigs <- append(sigs ,'***')
              x_sigs <- append(x_sigs, mean(c(xu, xi)))
            } else if (pv <0.001 & pv >0.0001){
              sigs <- append(sigs ,'**')
              x_sigs <- append(x_sigs, mean(c(xu, xi)))
            } else if (pv <0.05 & pv >0.001){
              sigs <- append(sigs ,'*')
              x_sigs <- append(x_sigs, mean(c(xu, xi)))
            }
            c=c+1
          }
        }
      }
    }
    mix <- mx+0.05*(c+1)*mx
    vioplot1(dta, col=as.character(col_prov$Color[match(names(dta), col_prov$Province)]), names=names(dta), ylim = c(0, mix), cex = 1.7)
    text(n, x = ceiling(length(dta)/2), y=mix, cex=2)
    if (length(xs0)>0){
      segments(x0=xs0, y0=ys0, x1=xs1, y1=ys1, lwd=2)
      text(x = x_sigs, y = ys0+0.025*mx, labels = sigs, cex=2)
    }
    
    stripchart(dta, add = T, pch=19, vertical=T, method='jitter')
    co=co+1
  }
}
dev.off()

comps_ca1 <- comps_ca
saveRDS(comps_ca1, 'provinces_carbon_characteristics_absolute.rds')
saveRDS(comps_cam, 'provinces_carbon_characteristics_absolute_mean.rds')
norms <- NULL
for (i in 1:4){norms <- append(norms, max(save_dt[[i]]))}
for (fr in frcs){
  for (pr in names(comps_ca$carb[[fr]])){
    comps_ca$carb[[fr]][[pr]] <- comps_ca$carb[[fr]][[pr]]/norms
    comps_ca$carb[[fr]][[pr]][is.na(comps_ca$carb[[fr]][[pr]])] <- 0
  }
}
colnames(sig_carb)<- c('Trans', 'p.val', 'Clade')
sig_carb <- as.data.frame(sig_carb)
saveRDS(comps_ca, 'provinces_carbon_characteristics.rds')
saveRDS(sig_carb, 'significant_differences_carbon_wilcox.rds')
