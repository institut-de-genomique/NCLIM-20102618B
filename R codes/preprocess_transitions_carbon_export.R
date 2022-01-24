source('axis_map0.R')
source('hide_arctic.R')
trans <- readRDS('transitions.rds')
sig_hex <- read.table('significant_differences_hexanauplia_wilcox.txt', header = T)
sig_hex <- sig_hex[grepl('Micro',sig_hex$Clade) | grepl('Meso',sig_hex$Clade),]
comp_all <- readRDS('SMAGs_provinces_compositional_functions_genomics_trophic.rds')
comp_hex <- comp_all$`Abundance Hex_bis`
tax_s <- readRDS('SMAGs_Statistics_02.rds')
Hex_list_bis <- unique(tax_s$Hex_bis)
Hex_list_bis <- Hex_list_bis[!is.na(Hex_list_bis) & Hex_list_bis!='0']
Hex_list_bis <- Hex_list_bis[c(4, 3, 2, 1, 5)]
letters<- c('A', 'B', 'C', 'D', 'E', 'F')
fracs <- c('0-0.2', '0.22-3', '0.8-5', '5-20', '20-180', '180-2000')
dt <- NULL
vals <- NULL
frcs <- NULL
i=1
for (t in sig_hex$Trans){
  fr <- fracs[sapply(letters, FUN = function(x){grepl(x, t)})]
  p1 <- strsplit(t, '->')[[1]][1]
  p2 <- strsplit(t, '->')[[1]][2]
  if (comp_hex[[fr]][[p1]][which(Hex_list_bis==sig_hex$Clade[i])] > comp_hex[[fr]][[p2]][which(Hex_list_bis==sig_hex$Clade[i])]){
    v <- '-'
  } else{
    v <- '+'
  }
  vals <- c(vals,log10(comp_hex[[fr]][[p2]][which(Hex_list_bis==sig_hex$Clade[i])]+1)-log10(comp_hex[[fr]][[p1]][which(Hex_list_bis==sig_hex$Clade[i])]+1) )
  dt <-c(dt, v)
  frcs <- append(frcs, fr)
  i=i+1
}
sig_hex$sign <- dt
sig_hex$value <- vals
sig_hex$frac <- frcs
saveRDS(sig_hex, 'transitions_hexanauplia_sizes.rds')

sig_diaz <- read.table('significant_differences_diazotrophy_wilcox.txt', header = T)
sig_diaz <- sig_diaz[grepl('Cyanobacteria',sig_diaz$Clade) ,]
comp_all_d <- readRDS('MAGprok_provinces_compositional_functions_diazotrophy.rds')
comp_diaz <- comp_all_d$`Abundance Diazotrophy`
tax <- read.table('MAGprok_taxo_1.txt', sep='\t', header = T)
tax$N_Fixation_Phylum <- as.character(levels(tax$N_Fixation_Phylum))[tax$N_Fixation_Phylum]
Diaz_list <- unique(tax$N_Fixation_Phylum)
Diaz_list <- Diaz_list[Diaz_list!=""]
letters<- c('A', 'B', 'C', 'D', 'E', 'F')
fracs <- c('0-0.2', '0.22-3', '0.8-5', '5-20', '20-180', '180-2000')
dt <- NULL
vals <- NULL
frcs <- NULL
i=1
for (t in sig_diaz$Trans){
  fr <- fracs[sapply(letters, FUN = function(x){grepl(x, t)})]
  p1 <- strsplit(t, '->')[[1]][1]
  p2 <- strsplit(t, '->')[[1]][2]
  if (comp_diaz[[fr]][[p1]][which(Diaz_list==sig_diaz$Clade[i])] > comp_diaz[[fr]][[p2]][which(Diaz_list==sig_diaz$Clade[i])]){
    v <- '-'
  } else{
    v <- '+'
  }
  vals <- c(vals,
            log10(comp_diaz[[fr]][[p2]][which(Diaz_list==sig_diaz$Clade[i])]+1)-log10(comp_diaz[[fr]][[p1]][which(Diaz_list==sig_diaz$Clade[i])]+1) )
  dt <-c(dt, v)
  frcs <- append(frcs, fr)
  i=i+1
}
sig_diaz$sign <- dt
sig_diaz$value <- vals
sig_diaz$Clade <- paste('Diazotrophic', sig_diaz$Clade)
sig_diaz$frac <- frcs
saveRDS(sig_diaz, 'transitions_diazotrophy_cyanobacteria.rds')

sig_phot <- read.table('significant_differences_phototrophy_wilcox.txt', header = T)
sig_phot <- sig_phot[grepl('Pico',sig_phot$Clade) | grepl('Nano',sig_phot$Clade ) ,]
comp_all_p <- readRDS('all_MAGs_provinces_compcompositional_functions_phyto.rds')
comp_phot <- comp_all_p
tax_a <- readRDS('all_MAGs_statistics.rds')
tax_a <- data.frame(lapply(tax_a, as.character), stringsAsFactors=FALSE)
Phyto_list =unique(tax_a$Phyto_bis)
Phyto_list=Phyto_list[Phyto_list!='0']
#ordre=c(2, 4, 1, 3, 5, 6, 8, 11, 12, 9, 7, 10)
phot_list <- Phyto_list
letters<- c('A', 'B', 'C', 'D', 'E', 'F')
fracs <- c('0-0.2', '0.22-3', '0.8-5', '5-20', '20-180', '180-2000')
dt <- NULL
vals <- NULL
frcs <- NULL
i=1
for (t in sig_phot$Trans){
  fr <- fracs[sapply(letters, FUN = function(x){grepl(x, t)})]
  p1 <- strsplit(t, '->')[[1]][1]
  p2 <- strsplit(t, '->')[[1]][2]
  if (comp_phot[[fr]][[p1]][which(phot_list==sig_phot$Clade[i])] > comp_phot[[fr]][[p2]][which(phot_list==sig_phot$Clade[i])]){
    v <- '-'
  } else{
    v <- '+'
  }
  vals <- c(vals,
            log10(comp_phot[[fr]][[p2]][which(phot_list==sig_phot$Clade[i])]+1)- log10(comp_phot[[fr]][[p1]][which(phot_list==sig_phot$Clade[i])]+1))
  dt <-c(dt, v)
  frcs <- c(frcs, fr)
  i=i+1
}
sig_phot$sign <- dt
sig_phot$value <- vals
sig_phot$frac <- frcs
sig_phot <- sig_phot[sig_phot$frac %in% c('0.22-3', '0.8-5', '5-20'),]
saveRDS(sig_phot, 'transitions_phototrophy.rds')

sig_phot_bis <- read.table('significant_differences_phototrophy_bis_wilcox.txt', header = T)
sig_phot_bis <- sig_phot_bis[grepl('Pico',sig_phot_bis$Clade) | grepl('Nano',sig_phot_bis$Clade ) | grepl('Micro',sig_phot_bis$Clade ) ,]
comp_all_p <- readRDS('all_MAGs_provinces_compcompositional_functions_phyto_bis.rds')
comp_phot <- comp_all_p
tax_a <- readRDS('all_MAGs_statistics.rds')
tax_a <- data.frame(lapply(tax_a, as.character), stringsAsFactors=FALSE)
Phyto_list =unique(tax_a$Phyto_ter)
Phyto_list=Phyto_list[Phyto_list!='0']
#ordre=c(2, 4, 1, 3, 5, 6, 8, 11, 12, 9, 7, 10)
phot_list <- Phyto_list
letters<- c('A', 'B', 'C', 'D', 'E', 'F')
fracs <- c('0-0.2', '0.22-3', '0.8-5', '5-20', '20-180', '180-2000')
dt <- NULL
vals <- NULL
frcs <- NULL
i=1
for (t in sig_phot_bis$Trans){
  fr <- fracs[sapply(letters, FUN = function(x){grepl(x, t)})]
  p1 <- strsplit(t, '->')[[1]][1]
  p2 <- strsplit(t, '->')[[1]][2]
  if (comp_phot[[fr]][[p1]][which(phot_list==sig_phot_bis$Clade[i])] > comp_phot[[fr]][[p2]][which(phot_list==sig_phot_bis$Clade[i])]){
    v <- '-'
  } else{
    v <- '+'
  }
  vals <- c(vals,
            log10(comp_phot[[fr]][[p2]][which(phot_list==sig_phot_bis$Clade[i])]+1)- log10(comp_phot[[fr]][[p1]][which(phot_list==sig_phot_bis$Clade[i])]+1))
  dt <-c(dt, v)
  frcs <- c(frcs, fr)
  i=i+1
}
sig_phot_bis$sign <- dt
sig_phot_bis$value <- vals
sig_phot_bis$frac <- frcs
sig_phot_bis <- sig_phot_bis[sig_phot_bis$frac %in% c('0.22-3', '0.8-5', '5-20', '20-180', '180-2000'),]
saveRDS(sig_phot_bis, 'transitions_phototrophy_bis.rds')

sigs <- rbind(sig_hex, sig_diaz, sig_phot_bis)

sigs$Clade<- as.character(levels(sigs$Clade))[sigs$Clade]
sigs$Clade[grepl('Mesozooplankton', sigs$Clade)] <- 'Mesozooplankton'
sigs$Clade[grepl('Microzooplankton', sigs$Clade)] <- 'Microzooplankton'
unique_tr <- unique(paste(sigs$Clade, sigs$sign))
sigs$type <- paste(sigs$Clade, sigs$sign)
sigs$pch <- sapply(sigs$type, FUN = function(x){which(unique_tr==x)})
saveRDS(sigs,'significant_transitions_wilcox_sizes_bis.rds')

sigsb <- rbind(sig_hex, sig_diaz, sig_phot)

sigsb$Clade<- as.character(levels(sigsb$Clade))[sigsb$Clade]
sigsb$Clade[grepl('Mesozooplankton', sigsb$Clade)] <- 'Mesozooplankton'
sigsb$Clade[grepl('Microzooplankton', sigsb$Clade)] <- 'Microzooplankton'
unique_trb <- unique(paste(sigsb$Clade, sigsb$sign))
sigsb$type <- paste(sigsb$Clade, sigsb$sign)
sigsb$pch <- sapply(sigsb$type, FUN = function(x){which(unique_trb==x)})
saveRDS(sigsb,'significant_transitions_wilcox_sizes.rds')

coords<- readRDS('new_stations_1deg.rds')
coords <- as.data.frame(coords)
colnames(coords) <- c('lat', 'long')
cols <- c('darkgoldenrod4', 'darkgreen', 'cyan', 'coral4', 'blueviolet', 'deepskyblue', 'darkorchid1', 'forestgreen', 'deeppink',
          'dodgerblue','darkblue' , 'coral3', 'darkmagenta','darksalmon', 'chartreuse4', 'brown3', 'cyan2', 'firebrick1', 'darkturquoise','gold')
rev_mapping <- readRDS('rev_mapping_lo_lt.rds')
bc_new <- readRDS('model-mean_bray-curtis_dom.rds')
bc_new <- bc_new[rev_mapping]
cond_bc <- bc_new>1/6
make_map <- function(name, sel, sel0, bc){
  unique_tr0 <- unique_tr[sel0]
  cols0 <- cols[sel0]
  pdf(family="Helvetica",file=paste('transitions_',name,'.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T, cex=0.284)
  j=1
  passed <- NULL
  for (t in sigs$Trans[sel]){
    frc <- fracs[sapply(letters, FUN = function(x){grepl(x, t)})]
    cond_oor <- trans[[which(rev(fracs)==frc)]]
    cond_oor[!grepl(t, cond_oor)] <- F 
    cond_oor[grepl(t, cond_oor)] <-T
    cond_oor<- as.logical(cond_oor)
    if (bc==T){
      cond_oor[!cond_bc]<- F
    }
    co <- 1
    if (t %in% passed){
      h <- 2
      k <- 4
    } else{
      h<- 1
      k<-3
    }
    for (i in unique(coords$long)){
      q <- sum(coords$long==i)
      if (co%%4==h ){
        cond_oor[coords$long==i][c(1:q)[!(seq(1,q,1) %in% seq(1, q,4))]] <- F
      } else if(co%%4==k){
        cond_oor[coords$long==i][c(1:q)[!(seq(1,q,1) %in% seq(3, q,4))]] <- F
      } else {
        cond_oor[coords$long==i] <- F
      } 
      co <- co+1
    }
    points(x=coords$long[cond_oor], y =coords$lat[cond_oor], pch=sigs$pch[sel][j],
           col=scales::alpha(cols0[which(unique_tr0==sigs$type[sel][j])], 0.7),cex=0.2)
    passed <- c(passed, t)
    j=j+1
  }
  hide_arctic()
  axis_map0()
  ordre <- order(unique_tr0)
  plot(0,0, col='white',  ylim=c(0,length(cols0)),axes=FALSE, frame.plot=F, xlab = '', ylab='')
  text(x=rep(0.5, length(cols0)), y = 1:length(cols0),  labels = unique_tr0[ordre])
  points(x=rep(0, length(cols0)), y = 1:length(cols0),pch=which(sel0==T)[ordre],col=scales::alpha(cols0[ordre]) )
  dev.off()
}
make_map(name = 'diaz', sel=grepl('Diazotrophic Cyano', sigs$Clade), sel0=grepl('Diazotrophic Cyano', unique_tr), bc=F)

make_map(name = 'hex', sel=grepl('zooplankton', sigs$Clade), sel0=grepl('zooplankton', unique_tr), bc=F)

make_map(name = 'photo_cyano',
         sel=grepl('Cyano', sigs$Clade) & !grepl('Diazotrophic Cyano', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans)), 
         sel0=grepl('Cyano', unique_tr) & !grepl('Diazotrophic Cyano', unique_tr), bc=F)

make_map(name = 'photo_algae', sel=grepl('Algae', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans)) ,
         sel0=grepl('Algae', unique_tr), bc=F)

make_map(name= 'photo_cyano_big',
         sel=grepl('Cyano', sigs$Clade) & !grepl('Diazotrophic Cyano', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans)),
         sel0=grepl('Cyano', unique_tr) & !grepl('Diazotrophic Cyano', unique_tr), bc=F)

make_map(name = 'photo_algae_big', sel=grepl('Algae', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans))  ,
         sel0=grepl('Algae', unique_tr), bc=F)

make_map(name = 'photo_diatom_big', sel=grepl('Diatom', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans))  ,
         sel0=grepl('Diatom', unique_tr), bc=F)


make_map(name = 'diaz_bc', sel=grepl('Diazotrophic Cyano', sigs$Clade), sel0=grepl('Diazotrophic Cyano', unique_tr), bc=T)

make_map(name = 'hex_bc', sel=grepl('zooplankton', sigs$Clade), sel0=grepl('zooplankton', unique_tr), bc=T)

make_map(name = 'photo_cyano_bc',
         sel=grepl('Cyano', sigs$Clade) & !grepl('Diazotrophic Cyano', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans)), 
         sel0=grepl('Cyano', unique_tr) & !grepl('Diazotrophic Cyano', unique_tr), bc=T)

make_map(name = 'photo_algae_bc', sel=grepl('Algae', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans)) ,
         sel0=grepl('Algae', unique_tr), bc=T)

make_map(name = 'photo_diatom_bc', sel=grepl('Diatom', sigs$Clade) & (grepl('B', sigs$Trans) | grepl('C', sigs$Trans) | grepl('D', sigs$Trans)) , 
         sel0=grepl('Diatom', unique_tr), bc=T)

make_map(name= 'photo_cyano_big_bc',
         sel=grepl('Cyano', sigs$Clade) & !grepl('Diazotrophic Cyano', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans)),
         sel0=grepl('Cyano', unique_tr) & !grepl('Diazotrophic Cyano', unique_tr), bc=T)

make_map(name = 'photo_algae_big_bc', sel=grepl('Algae', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans))  ,
         sel0=grepl('Algae', unique_tr), bc=T)

make_map(name = 'photo_diatom_big_bc', sel=grepl('Diatom', sigs$Clade) & (grepl('E', sigs$Trans) | grepl('F', sigs$Trans))  ,
         sel0=grepl('Diatom', unique_tr), bc=T)
#make_map(name = 'all', sel=rep(T, length(sigs$Clade)))

