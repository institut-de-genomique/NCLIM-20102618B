library('matlab')
source('axis_map0.R')
source('hide_arctic.R')
source('vioplot1.R')
setwd('/env/cns/scratch_TaraOcean/Niches_models_Neogen/final_model_5/')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')

letters <- c('A', 'B', 'C', 'D', 'E', 'F')
frcs <- c('0-0.2', '0.22-3','0.8-5', '5-20', '20-180', '180-2000')

cols <-  colorRampPalette(c('red','darkred','white','darkgreen','green'))(100)
coords<- readRDS('new_stations_1deg.rds')
coords <- as.data.frame(coords)
colnames(coords) <- c('lat', 'long')
coords$cell <- paste(coords$lat, coords$long, sep='_')

longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
longi_sorted =sort(longi)
lati<-seq(-89.5, 60, 1)
dom_2006 <- readRDS('dominant_communities_2006.rds')
dom_2090 <- readRDS('dominant_communities_2090.rds')
dom_woa<- readRDS('dominant_communities_woa.rds')

sigs <- readRDS('significant_transitions_wilcox_sizes_bis.rds')
unique_tr<-unique(sigs$type)

preds <- readRDS('predictions_new_stations_model-mean_1deg.rds')
ann <- read.table('niches_annotation.txt')
ann$prov <- paste(letters[match(ann$V1, frcs)], ann$V2, sep='')
colnames(preds)<- ann$prov
for (l in letters){
  su <- apply(preds[,grepl(l, colnames(preds))], 1, sum)
  preds[,grepl(l, colnames(preds))] <- preds[,grepl(l, colnames(preds))]/su
}

bray_curtis_fish_eez <- readRDS('model-mean_bray_curtis_cond_fish_eez.rds')
bray_curtis_fish_eez$cond_fish[is.na(bray_curtis_fish_eez$cond_fish)] <- F
bray_curtis_fish_eez$cond_eez[is.na(bray_curtis_fish_eez$cond_eez)] <- F
cond_eez <- bray_curtis_fish_eez$cond_eez
cond_fish <- bray_curtis_fish_eez$cond_fish


filled_concentric_circles <- function(x,y,r,nsteps=100,off=NULL,...){
  rs <- seq(0,2*pi,len=nsteps)
  xc <- x-(r+off)*cos(rs)
  yc <- y-(r+off)*sin(rs)
  yc <- c(yc,  y+(r)*sin(rs))
  xc <- c(xc,  x-(r)*cos(rs))
  graphics::polygon(xc,yc,border = NA, ...)
}

add_taxo <- function(name, dat){
  cols <-  c('darkgoldenrod4', 'darkgreen', 'cyan', 'coral4', 'blueviolet', 'deepskyblue', 'darkorchid1', 'forestgreen', 'deeppink',
             'dodgerblue','darkblue' , 'coral3', 'turquoise1','green2' , 'chartreuse4', 'brown3', 'cyan3', 'firebrick1', 'darkturquoise','gold')
  
  trans <- readRDS('transitions.rds')
  if (name=='hex'){
    sel=grepl('zooplankton', sigs$Clade)
    sel0=grepl('zooplankton', unique_tr)
  } else if (name=='diaz'){
    sel=grepl('Diazotrophic Cyano', sigs$Clade)
    sel0=grepl('Diazotrophic Cyano', unique_tr)
  } else if (name=='phot_algae'){
    sel=grepl('Algae', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans)) 
    sel0=grepl('Algae', unique_tr)
  } else if (name=='phot_diatom'){
    sel=grepl('Diatom', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans))  
    sel0=grepl('Diatom', unique_tr)
  } else if (name=='phot_cyano'){
    sel=grepl('Cyano', sigs$Clade) & !grepl('Diazotrophic Cyano', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans))
    sel0=grepl('Cyano', unique_tr) & !grepl('Diazotrophic Cyano', unique_tr)
  } else if (name=='phot_algae_big'){
    sel=grepl('Algae', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans))
    sel0=grepl('Algae', unique_tr)
  } else if (name=='phot_diatom_big'){
    sel=grepl('Diatom', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans))
    sel0=grepl('Diatom', unique_tr)
  } else if (name=='phot_cyano_big'){
    sel=grepl('Cyano', sigs$Clade) & !grepl('Diazotrophic Cyano', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans))
    sel0=grepl('Cyano', unique_tr) & !grepl('Diazotrophic Cyano', unique_tr)
  }
  rev_mapping <- readRDS('rev_mapping_lo_lt.rds')
  bc_new <- readRDS('model-mean_bray-curtis_dom.rds')
  bc_new <- bc_new[rev_mapping]
  cond_bc <- bc_new>1/6
  
  unique_tr0 <- unique_tr[sel0]
  #sizes_c <- c(1:length(unique_tr0)-0)
  #write(match(sigs$Trans[sel], unique_tr0) , 'test.txt')
  if (name != 'diaz'){
     if (name != 'hex'){
       sigs$sizes_c[sel][grepl('Pico', sigs$type[sel]) & grepl('-', sigs$type[sel])] <- 1
       sigs$sizes_c[sel][grepl('Pico', sigs$type[sel]) & !grepl('-', sigs$type[sel])] <- 2
       sigs$sizes_c[sel][grepl('Nano', sigs$type[sel]) & grepl('-', sigs$type[sel])] <- 3
       sigs$sizes_c[sel][grepl('Nano', sigs$type[sel]) & !grepl('-', sigs$type[sel])] <- 4
       sigs$sizes_c[sel][grepl('Micro', sigs$type[sel]) & grepl('-', sigs$type[sel])] <- 5
       sigs$sizes_c[sel][grepl('Micro', sigs$type[sel]) & !grepl('-', sigs$type[sel])] <- 6
     } else{
       sigs$sizes_c[sel][grepl('Micro', sigs$type[sel]) & grepl('-', sigs$type[sel])] <- 3
       sigs$sizes_c[sel][grepl('Micro', sigs$type[sel]) & !grepl('-', sigs$type[sel])] <- 4
       sigs$sizes_c[sel][grepl('Meso', sigs$type[sel]) & grepl('-', sigs$type[sel])] <- 5
       sigs$sizes_c[sel][grepl('Meso', sigs$type[sel]) & !grepl('-', sigs$type[sel])] <- 6
     }
     #sizes_c <- unique(sigs$sizes_c[sel])   
  } else{
    sizes_c <- c(1:length(unique_tr0))+1
    #print(sizes_c)
    sigs$sizes_c[sel] <- sizes_c[match(sigs$type[sel], unique_tr0)]
  }
  #sigs$sizes_c[sel] <- sizes_c[match(sigs$type[sel], unique_tr0)]
  cols0 <- cols[sel0]
  j=1
  passed <- NULL
  for (t in sigs$Trans[sel]){
    frc <- frcs[sapply(letters, FUN = function(x){grepl(x, t)})]
    cond_oor <- trans[[which(rev(frcs)==frc)]]
    cond_oor[!grepl(t, cond_oor)] <- F 
    cond_oor[grepl(t, cond_oor)] <-T
    cond_oor<- as.logical(cond_oor)
    cond_oor[!cond_bc | is.na(dat)]<- F
    co <- 1
    #if (t %in% passed){
    #  h <- 2
    #  k <- 4
    #} else{
    #  h<- 1
    #  k<-3
    #}
    h<- 1
    k<- 4
    n <- 6
    for (i in unique(coords$long)){
      q <- sum(coords$long==i)
      if (co%%n==h ){
        cond_oor[coords$long==i][c(1:q)[!(seq(1,q,1) %in% seq(h, q,n))]] <- F
      } else if(co%%n==k){
        cond_oor[coords$long==i][c(1:q)[!(seq(1,q,1) %in% seq(k, q,n))]] <- F
      } else {
        cond_oor[coords$long==i] <- F
      } 
      co <- co+1
    }
    #points(x=coords$long[cond_oor], y =coords$lat[cond_oor], pch=sigs$pch[sel][j],
    #       col=scales::alpha(cols0[which(unique_tr0==sigs$type[sel][j])], 0.7),cex=0.2)
    #points(x=coords$long[cond_oor], y =coords$lat[cond_oor],pch=1,lwd=0.3,
    #       col=scales::alpha(cols0[which(unique_tr0==sigs$type[sel][j])], 1),cex=(sigs$pch[sel][j])/30 )
    for (l in 1:length(coords$long[cond_oor])){
      filled_concentric_circles(x=coords$long[cond_oor][l], y=coords$lat[cond_oor][l], r=sigs$sizes_c[sel][j]/4, off=0.25, col=cols0[which(unique_tr0==sigs$type[sel][j])])
    }
    passed <- c(passed, t)
    j=j+1
  }
  ordre <- order(unique_tr0)
  if (name != 'diaz'){
    sizes_c <- rep(NA, length(unique_tr0))
    if (name !='hex'){
      sizes_c[grepl('Pico', unique_tr0) & grepl('-', unique_tr0)] <- 1
      sizes_c[grepl('Pico', unique_tr0) & !grepl('-', unique_tr0)] <- 2
      sizes_c[grepl('Nano', unique_tr0) & grepl('-', unique_tr0)] <- 3
      sizes_c[grepl('Nano', unique_tr0) & !grepl('-', unique_tr0)] <- 4
      sizes_c[grepl('Micro', unique_tr0) & grepl('-', unique_tr0)] <- 5
      sizes_c[grepl('Micro', unique_tr0) & !grepl('-', unique_tr0)] <- 6
    } else {
      sizes_c[grepl('Micro', unique_tr0) & grepl('-', unique_tr0)] <- 3
      sizes_c[grepl('Micro', unique_tr0) & !grepl('-', unique_tr0)] <- 4
      sizes_c[grepl('Meso', unique_tr0) & grepl('-', unique_tr0)] <- 5
      sizes_c[grepl('Meso', unique_tr0) & !grepl('-', unique_tr0)] <- 6
    }
    print(unique_tr0)
    print(grepl('Pico', unique_tr0) & grepl('-', unique_tr0))
    print(grepl('Pico', unique_tr0) & !grepl('-', unique_tr0))
    print(sizes_c)
  }
  plot(0,0, col='white',  ylim=c(0,length(cols0)+3),xlim=c(-3, length(cols0)-1),axes=FALSE, frame.plot=F, xlab = '', ylab='')
  text(x=rep(0.5, length(cols0)), y = 1:length(cols0),  labels = unique_tr0[ordre], cex=2)
  for (l in 1:length(cols0)){
    filled_concentric_circles(x=-3, y=l, col=cols0[ordre][l],r=sizes_c[ordre][l]/64, off=0.015625)
  }
  #points(x=rep(0, length(cols0)), y = 1:length(cols0),pch=which(sel0==T)[ordre],col=scales::alpha(cols0[ordre]), cex=2 )
}


make_carbon_bilan <- function(type){
  sel <- seq(11,1,-2)
  sels <- list(sel[c(1,2,3,4,5,6)], sel[c(1,2,3,5,6)])
  sel0s <- list(c(1,2,3,4,5,6), c(1,2,3,5,6))
  versions <- c('_all', '_no5-20')
  dtosave <- NULL
  for (i in c(1, 2)){
    sel <- sels[[i]]
    sel0<-sel0s[[i]]
    ver <- versions[i]
    ass <- function(x, sel, sel0){
      a <- paste0(letters[sel0],x[sel] )
      #print(a)
      paste(a, collapse='_')
    }
    
    
    ass_woa <- apply(dom_woa, 1, FUN = ass, sel=sel, sel0=sel0)
    ass_2006 <- apply(dom_2006, 1, FUN = ass, sel=sel, sel0=sel0)
    ass_2090 <- apply(dom_2090, 1, FUN = ass, sel=sel, sel0=sel0)
    
    if (type==''){
      cexp_data <- read.table('carbon_flux_Henson_2019.txt', sep='\t', header = T)
    } else if (type=='_guidi'){
      cexp_data <- read.table('hydrofull_for_margaux.txt', header = T, sep='\t')
      colnames(cexp_data)[2] <- 'Latitude'
      colnames(cexp_data)[3] <- 'Longitude'
      colnames(cexp_data)[67] <- 'Export.flux'
    }else if (type=='_extrapolated_henson'){
      cexp_data <- readRDS('Henson_2012_export.rds')
    } else if (type=='_extrapolated_eppley'){
      cexp_data <- readRDS('eppley_export.rds')
    } else if (type=='_extrapolated_laws'){
      cexp_data <- readRDS('laws_export.rds')
    } else if (type=='_extrapolated_schlitzer'){
      cexp_data <- readRDS('schlitzer_export.rds')
    }
    name <- 'Export.flux'
    #cexp_data <- read.table('LeMoigne_2013_carbonflux.tab', sep='\t', header = T)
    # name <- TOC.flux..mg.m..2.day.
    cexp_data <- cexp_data[cexp_data$Latitude<60,]
    cexp_data$lt_cell <- sapply(cexp_data$Latitude, FUN = function(x){lati[which.min(abs(x-lati))]})
    cexp_data$lo_cell <- sapply(cexp_data$Longitude, FUN = function(x){longi[which.min(abs(x-longi))]})
    cexp_data$cell <- paste(cexp_data$lt_cell, cexp_data$lo_cell, sep='_')
    cexp_data$assemblage <- ass_2006[match(cexp_data$cell, coords$cell)]
    cexp_data$assemblage90 <- ass_2090[match(cexp_data$cell, coords$cell)]
    
    #wilc_test <-pairwise.wilcox.test(cexp_data$TOC.flux..mg.m..2.day., cexp_data$assemblage)
    counts <- as.data.frame(table(cexp_data$assemblage))
    cexp_data<- cexp_data[cexp_data$assemblage %in% counts$Var1[counts$Freq>2],]
    #wilc_test <-pairwise.wilcox.test(cexp_data$TOC.flux..mg.m..2.day., cexp_data$assemblage, p.adjust.method = 'fdr')
    
    
    
    dt <- list()
    asse <- NULL
    for (a in unique(cexp_data$assemblage)){
      carbs <- cexp_data[[name]][cexp_data$assemblage==a]
      #print(carbs)
      if (length(carbs[!is.na(carbs)])>2){
        dt[[a]]<- carbs[!is.na(carbs)]
        asse <- append(asse, a)
      }
    }
    cexp_data<- cexp_data[cexp_data$assemblage %in% asse,]
    
    cexp_data1 <- cexp_data[cexp_data$assemblage90 %in% unique(cexp_data$assemblage),]
    cexp_data1 <- cexp_data1[!is.na(cexp_data1$Export.flux),]
    w_coso <- cos(cexp_data1$Latitude*2*pi/360)
    if (type=='_extrapolated_henson'){
      data_expo <- sum(as.numeric(111*111*w_coso*cexp_data1$Export.flux*1.e6))
    } else{
      data_expo <- sum(as.numeric(111*111*w_coso*cexp_data1$Export.flux*365*1.e6/1.e3))
    }
    write(data_expo, paste('carbon_export_bilan_restricted',type,'_',ver,'.txt', sep=''))
    
    wilc_test <-pairwise.wilcox.test(cexp_data[[name]], cexp_data$assemblage, p.adjust.method = 'holm')
    
    aova_ass <- aov(cexp_data[[name]]~cexp_data$assemblage)
    print(AIC(aova_ass))
    if (i==1){
     jk <- 6
    } else{
     jk <- 5
    }
    dot <- NULL 
    for (i in 1:jk){
      vo <- sapply(cexp_data$assemblage, FUN = function(x){strsplit(x, split = '_')[[1]][i]})
      dot <- cbind(dot, vo)
    }
    for (koi in 1:jk){
      aova <- aov(cexp_data[[name]]~dot[,koi])
      print(AIC(aova))
    }
    
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    u <- wilc_test$p.value
    print(type)
    print(sum(sum(u<0.05, na.rm=T)))
    print(sum(sum(!is.na(u))))
    print(sum(sum(u<0.05, na.rm=T))/sum(sum(!is.na(u))))

    mx <- max(unlist(dt))
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
            xu <- which(names(dt)==rownames(u)[i])
            xi <- which(names(dt)==colnames(u)[j])
            xs0 <- c(xs0,xu)
            xs1 <- c(xs1,xi)
            ys0 <- c(ys0, mx+0.05*c*mx)
            ys1 <- c(ys1, mx+0.05*c*mx)
    #         sig_carb <- rbind(sig_carb , c(paste(rownames(u)[i], colnames(u)[j], sep='->'), pv, n))
    #         sig_carb <- rbind(sig_carb , c(paste(colnames(u)[j], rownames(u)[i] , sep='->'), pv, n))
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
    
    ass<- NULL
    for (a in unique(cexp_data$assemblage)){
      ass <- append(ass, c(a, median(cexp_data[[name]][cexp_data$assemblage==a], na.rm=T),
                    mean(cexp_data[[name]][cexp_data$assemblage==a], na.rm=T)))
    }
    ass <- matrix(ass, ncol=3, byrow = T)
    ass <-as.data.frame(ass, stringsAsFactors = F)
    saveRDS(ass, paste('assemblages_carbon_export',ver,type,'.rds', sep=''))
    
    pdf(paste('carbon_export_violin_asssemblages',ver,type,'.pdf', sep=''), height = 12, width = 12)
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    vioplot1(dt, col=col_vector[1:length(dt)], names=NA,ylim = c(0, mix), cex = 1.7, xaxt='n', ylab='TOC flux (mgC/m2/d)')
    axis(1, at = 1:length(dt), las=2, names(dt))
    if (length(xs0)>0){
      segments(x0=xs0, y0=ys0, x1=xs1, y1=ys1, lwd=2)
      text(x = x_sigs, y = ys0+0.025*mx, labels = sigs, cex=2)
    }
    
    stripchart(dt, add = T, pch=19, vertical=T, method='jitter')
    dev.off()
    
    vals <- as.numeric(ass$V2[match(ass_2006, ass$V1)])
    vals_m <- as.numeric(ass$V3[match(ass_2006, ass$V1)])
    
    vals_f <- as.numeric(ass$V2[match(ass_2090, ass$V1)])
    vals_m_f <- as.numeric(ass$V3[match(ass_2090, ass$V1)])
    
    carb_map <- function(data_exp0,data_exp2=NULL, name, col_scale, sym, leg_title, cond, type, taxo=NULL){
      
      data_exp0 <- as.numeric(data_exp0)
      w_cos <- cos(coords$lat*2*pi/360)
      if (type=='_extrapolated_henson'){
        data_export <- as.numeric(111*111*w_cos*data_exp0*as.numeric(cond)*1.e6)
        data_expi <- as.numeric(data_exp0*as.numeric(cond)) 
      } else{
        data_export <- as.numeric(111*111*w_cos*data_exp0*365*as.numeric(cond)*1.e6/1.e3)
        data_expi <- as.numeric(data_exp0*as.numeric(cond)*365/1.e3)
      }
      data_exp0[!cond]<-NA
      if (sym==F){
        if (type=='_extrapolated_henson'){
          data_exp1 <- log10(data_exp0)
          data_exp3 <- log10(data_exp2)
        } else {
          data_exp1 <- log10(data_exp0*365/1.e3)
          data_exp3 <- log10(data_exp2*365/1.e3)
        }
      } else{
        if (type=='_extrapolated_henson'){
          data_exp1 <- data_exp0
        } else {
          data_exp1 <- data_exp0*365/1.e3
        }
      }
      
      if (sym!=T){
        mi <- min(c(data_exp1, data_exp3), na.rm = T)
        mx <- max(c(data_exp1, data_exp3), na.rm = T)
        data_norm <- round(99*(data_exp1-mi)/(mx-mi)+1)
      } else{
        mox <- max(c(abs(max(data_exp1, na.rm = T)), abs(min(data_exp1, na.rm = T))))
        data_norm <- round(99*(data_exp1+mox)/(2*mox)+1)
        mi <- -mox
        mx <- mox
      }
      
      data_contour <- matrix(NA, ncol=150, nrow=360)
      for (j in 1:length(coords$lat)){
        data_contour[which(longi_sorted==coords$long[j]),which(lati==coords$lat[j])]= as.numeric(data_norm[j])
      }
      pdf(family="Helvetica",file=paste(name, sep=''),
          width=10,height=4.065)
      par(mar=c(0,0,0,0))
      maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
      plot(coastline,lwd=0.0475, col='black', add=T, cex=0.284)
      .filled.contour(x=longi_sorted,y=lati, z=data_contour,  col=col_scale, levels=c(1:100))
      #print(data_contour)
      hide_arctic()
      axis_map0()
      if (!is.null(taxo)){
        add_taxo(taxo,dat=data_exp0)
      }
      pnt <- cbind(x =c(0,2,2,0), y =c(0,50,50,0))
      xmx <- 17
      plot(0,0, col='white', xlim=c(0,xmx), ylim=c(0, 60), axes=FALSE, frame.plot=F, xlab = '', ylab='')
      SDMTools::legend.gradient(pnt, col_scale, limits=c(round(mi,2),round(mx,2)), title =leg_title , cex=1)
      dev.off()
      return(list(data_export, data_expi))
    }
    rev_mapping <- readRDS('rev_mapping_lo_lt.rds')
    bc_new <- readRDS('model-mean_bray-curtis_dom.rds')
    bc_new <- bc_new[rev_mapping]
    cond_all <- rep(T, length(vals))
    cond_bc <- bc_new>1/6
    col_scale <- colorRampPalette(c('darkorchid4', 'blue', 'green', 'orange', 'red'))(100)
    col_scale0 <- colorRampPalette(c('red', 'orangered', 'orange','grey92','lightblue1' ,'blue', 'darkblue'))(100)
    dte_m_bc_cyano <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_cyanobacteria',ver,type,'.pdf', sep=''),
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type,
                              taxo='phot_cyano')

    te <- carb_map(data_exp0 = vals,data_exp2=vals_f,name = paste('carbon_flux_present_assemblages',ver,type,'.pdf', sep=''), 
                   col_scale = col_scale, sym=F,leg_title = 'Carbon flux log10(gC/m2/year)', cond = cond_all, type=type)
    te_m <- carb_map(data_exp0 = vals_m,data_exp2=vals_m_f,name = paste('carbon_flux_present_assemblages_mean',ver,type,'.pdf', sep=''), 
                     col_scale = col_scale, sym=F,leg_title = 'Carbon flux log10(gC/m2/year)', cond = cond_all, type=type)
    te_f <- carb_map(data_exp0 = vals_f,data_exp2=vals,name = paste('carbon_flux_futur_assemblages',ver,type,'.pdf', sep=''),
                     col_scale = col_scale, sym=F,leg_title = 'Carbon flux log10(gC/m2/year)',cond= cond_all, type=type)
    te_m_f <- carb_map(data_exp0 = vals_m_f,data_exp2=vals_m,name = paste('carbon_flux_futur_assemblages_mean',ver,type,'.pdf', sep=''),
                       col_scale = col_scale, sym=F,leg_title = 'Carbon flux log10(gC/m2/year)',cond =  cond_all, type=type)
                     
    
    dte <- carb_map(data_exp0 = vals_f-vals,name = paste('carbon_delta_flux_assemblages',ver,type,'.pdf', sep=''), col_scale=col_scale0, sym=T, 
                    leg_title = 'Delta Carbon flux (gC/m2/year)', cond = cond_all, type=type)
    dte_m <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean',ver,type,'.pdf', sep=''), 
                      col_scale = col_scale0, sym=T, 
                      leg_title = 'Delta Carbon flux (gC/m2/year)', cond=cond_all, type=type)
    dte_bc <- carb_map(data_exp0 = vals_f-vals,name = paste('carbon_delta_flux_assemblages_bc',ver,type,'.pdf', sep=''), col_scale=col_scale0, 
                       sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)', cond=cond_bc, type=type)   
    dte_m_bc <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc',ver,type,'.pdf', sep=''), 
                         col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type)
    saveRDS(dte_m_bc[[2]], paste('delta_fluxes', ver,type,'.rds', sep=''))
    dte_m_bc_diaz <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_cyano-diazotrophs',ver,type,'.pdf', sep=''), 
                         col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type, 
                         taxo='diaz')
    dte_m_bc_hex <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_copepods',ver,type,'.pdf', sep=''), 
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type,
                             taxo='hex')
    dte_m_bc_algae <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_algae',ver,type,'.pdf', sep=''), 
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type, 
                              taxo='phot_algae')
    dte_m_bc_diat <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_diatoms',ver,type,'.pdf', sep=''), 
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type, 
                              taxo='phot_diatom')
 #   dte_m_bc_cyano <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_cyanobacteria',ver,type,'.pdf', sep=''), 
  #                            col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type, 
   #                           taxo='phot_cyano')
    dte_m_bc_algae_b <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_algae_big',ver,type,'.pdf', sep=''),
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type,
                              taxo='phot_algae_big')
    dte_m_bc_diat_b <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_diatoms_big',ver,type,'.pdf', sep=''),
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type,
                              taxo='phot_diatom_big')
    dte_m_bc_cyano_b <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_cyanobacteria_big',ver,type,'.pdf', sep=''),
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc, type=type,
                              taxo='phot_cyano_big') 
    
    dte_fish <- carb_map(data_exp0 = vals_f-vals,name = paste('carbon_delta_flux_assemblages_fish',ver,type,'.pdf', sep=''), col_scale=col_scale0, sym=T, 
                    leg_title = 'Delta Carbon flux (gC/m2/year)', cond = cond_fish, type=type)
    dte_m_fish <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_fish',ver,type,'.pdf', sep=''), col_scale = col_scale0, sym=T, 
                      leg_title = 'Delta Carbon flux (gC/m2/year)', cond=cond_fish, type=type)
    dte_bc_fish <- carb_map(data_exp0 = vals_f-vals,name = paste('carbon_delta_flux_assemblages_bc_fish',ver,type,'.pdf', sep=''), col_scale=col_scale0, 
                       sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)', cond=cond_bc&cond_fish, type=type)
    dte_m_bc_fish <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_fish',ver,type,'.pdf', sep=''), 
                         col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc&cond_fish, type=type)
    
    dte_eez <- carb_map(data_exp0 = vals_f-vals,name = paste('carbon_delta_flux_assemblages_eez',ver,type,'.pdf', sep=''), col_scale=col_scale0, sym=T, 
                         leg_title = 'Delta Carbon flux (gC/m2/year)', cond = cond_eez, type=type)
    dte_m_eez <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_eez',ver,type,'.pdf', sep=''), col_scale = col_scale0, sym=T, 
                           leg_title = 'Delta Carbon flux (gC/m2/year)', cond=cond_eez, type=type)
    dte_bc_eez <- carb_map(data_exp0 = vals_f-vals,name = paste('carbon_delta_flux_assemblages_bc_eez',ver,type,'.pdf', sep=''), col_scale=col_scale0, 
                            sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)', cond=cond_bc&cond_eez, type=type)
    dte_m_bc_eez <- carb_map(data_exp0 = vals_m_f-vals_m,name = paste('carbon_delta_flux_assemblages_mean_bc_eez',ver,type,'.pdf', sep=''), 
                              col_scale = col_scale0, sym=T, leg_title = 'Delta Carbon flux (gC/m2/year)',cond =  cond_bc&cond_eez, type=type)
                      
    cexp06 <- sum(te[[1]], na.rm=T)
    cexp90 <- sum(te_f[[1]], na.rm=T)
    cexp06_m <- sum(te_m[[1]], na.rm=T)
    cexp90_m <- sum(te_m_f[[1]], na.rm=T)
    cexp06_r <- sum(te[[1]][!is.na(dte[[1]])], na.rm=T)
    cexp90_r <- sum(te_f[[1]][!is.na(dte[[1]])], na.rm=T)
    cexp06_m_r <- sum(te_m[[1]][!is.na(dte_m[[1]])], na.rm=T)
    cexp90_m_r <- sum(te_m_f[[1]][!is.na(dte_m[[1]])], na.rm=T)
    de <- sum(dte[[1]], na.rm=T)
    de_m <- sum(dte_m[[1]], na.rm=T)
    delt <- sum(dte[[1]], na.rm=T)*100/sum(te[[1]][!is.na(dte[[1]])], na.rm=T)
    delt_m <- sum(dte_m[[1]], na.rm=T)*100/sum(te_m[[1]][!is.na(dte_m[[1]])], na.rm=T)
    de_bc <- sum(dte_bc[[1]], na.rm=T)
    de_m_bc <- sum(dte_m_bc[[1]], na.rm=T)
    delt_bc <- sum(dte_bc[[1]], na.rm=T)*100/sum(te[[1]][!is.na(dte_bc[[1]])], na.rm=T)
    delt_m_bc <- sum(dte_m_bc[[1]], na.rm=T)*100/sum(te_m[[1]][!is.na(dte_m_bc[[1]])], na.rm=T)
    dtosave <- append(dtosave,c(ver, cexp06, cexp90, cexp06_m, cexp90_m,cexp06_r, cexp90_r, cexp06_m_r, cexp90_m_r, de, de_m, delt, delt_m,
                                de_bc, de_m_bc,delt_bc, delt_m_bc))
  }
  dtosave <- matrix(dtosave, nrow = 2, byrow = T)
  colnames(dtosave)<- c('version', 'export06_median(g/year)', 'export90_median(g/year)','export06_mean(g/year)', 'export90_mean(g/year)',
                        'export06_median(g/year)_restricted', 'export90_median(g/year)_restricted','export06_mean(g/year)_restricted', 
                        'export90_mean(g/year)_restricted', 'delta_median(gC/year)', 'delta_mean(gC/year)', 'delta_median(%)', 'delta_mean(%)',
                        'delta_median_bc(gC/year)', 'delta_mean_bc(gC/year)', 'delta_median_bc(%)', 'delta_mean_bc(%)')
  write.table(dtosave, paste('carbon_export_statistics',type,'.txt', sep=''), row.names = FALSE)
}


make_carbon_bilan(type='_extrapolated_henson')
make_carbon_bilan(type = '')
make_carbon_bilan(type='_guidi')
make_carbon_bilan(type = '_extrapolated_eppley')
make_carbon_bilan(type='_extrapolated_laws')
make_carbon_bilan(type='_extrapolated_schlitzer')
