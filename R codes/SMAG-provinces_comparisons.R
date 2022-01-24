library('gplots')
#library('readxl')
library('matlab')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_6")
cambria <- readRDS('cambria.rds')
type <- commandArgs(trailingOnly = T)[1]
version <- commandArgs(trailingOnly = T)[2]
data_SMAGs <- read.table(paste('Occurence_',type,'.txt',sep=''), header = T)
# if (type=='SMAGs'){
perc_smags <- read.table(paste(type,'_percent_mapped_read.txt', sep=''), header = T)
perc_smags<- perc_smags[,match(colnames(data_SMAGs), colnames(perc_smags))]
  #perc_smags <- perc_smags[,-1]
#}
if (type=='MAGprok'){
  frc_c <- c('GGZZ', 'GGMM', 'CCKK', 'QQSS', 'SSUU', 'MMQQ', 'CCII', 'KKQQ', 'QQRR')
  frc_s <- c('0.8-2000', '0.8-5', '0.22-3', '20-180', '180-2000', '5-20', '0.2-1.6','3-20', '20-RR')
  sizes <- c('Filter_size', rep(NA, dim(data_SMAGs)[2]-1))
  for (frc in frc_c){
    sizes[grepl(frc, colnames(data_SMAGs))] <- frc_s[which(frc_c==frc)]
  }
  data_SMAGs <- rbind(sizes, data_SMAGs)
}
nst <- dim(data_SMAGs)[2]-1
nmag <- dim(data_SMAGs)[1]-1
stat_SMAGs <- colnames(data_SMAGs)[2:(nst+1)]
smags <- as.character(levels(data_SMAGs[2:(nmag+1),1]))[data_SMAGs[2:(nmag+1),1]]
if (type=='MAGprok' & version=='old'){
  smags <- smags[!grepl('[[:digit:]]_MAG', smags)] # retrieve only published MAGs
  cond <- data_SMAGs[,1] %in% smags
  cond[1]<-T
  data_SMAGs <- data_SMAGs[cond,]
  nst <- dim(data_SMAGs)[2]-1
  nmag <- dim(data_SMAGs)[1]-1
}


extract_stat <- function(x){
  if (grepl('DCM', x)){
    v <- strsplit(x, split='DCM')[[1]][1]
    sepo <- 'DCM'
  } else{
    v <- strsplit(x, split='SUR')[[1]][1]
    sepo <- 'SUR'
  }
  j <- substr(v, 2, nchar(v))
  st <- paste(j, sepo, sep='_')
  return(st)
}

stat_SMAGs_all <- mapply(stat_SMAGs,FUN = extract_stat)
colnames(data_SMAGs) <- c('SMAG', stat_SMAGs_all)
colnames(perc_smags) <- c('SMAG', stat_SMAGs_all)

data_prov <- readRDS('Genocenoses_env_parameters_woa_scaled.rds')
data_prov$Fraction[data_prov$Fraction=='43952'] <- '5-20'
data_prov$lab <- paste(data_prov$Fraction, data_prov$Genocenose, sep='_')
to_remove <- readRDS('excluded_niches.rds')
data_prov <- data_prov[!(data_prov$lab %in% to_remove),]
valids <- readRDS('valids.rds')
data_prov <- data_prov[data_prov$lab %in% valids,]
stats_prov <- unique(data_prov$Station)

stat_SMAGs_prov <- stat_SMAGs_all[stat_SMAGs_all %in% stats_prov]
data_SMAGs0 <- data_SMAGs[,colnames(data_SMAGs) %in% stats_prov]
perc_smags0 <- perc_smags[,colnames(perc_smags) %in% stats_prov]
#perc_smags0 <- rbind(as.character(), matrix(perc_smags0))
colnames(data_SMAGs0) <-stat_SMAGs_prov 
rownames(data_SMAGs0)<- c('Filter_size', smags)
extract_frac <- function(x){
  return(as.character(x[[1]]))
}
test <- sapply(data_SMAGs[1,], FUN=extract_frac)
fracs <- unique(test)
if (type=='SMAGs'){
  fracs <- fracs[c(3,5,6,7,10)]
  letters <- c('B', 'C', 'F', 'E', 'D')
} else if (type=='MAGprok'){
  fracs <- fracs[c(2, 3, 5,6,7)]
  letters <- c('B', 'C', 'E', 'F', 'D')
}


jaccard_index <- function(x, occ){
  m11 <- sum(x==1 & occ==1)
  #m00 <- sum(x==1 & occ==1)
  m01 <- sum(x==0 & occ==1)
  m10 <- sum(x==1 & occ==0)
  jac <- m11/(m11+m01+m10)
  return(jac)
}

perc_stations <- function(x, occ){
  m11 <- sum(x==1 & occ==1)
  m10 <- sum(x==0 & occ==1)
  perc <- m11/(m11+m10)
  return(perc)
}

ab_stations <- function(x, occ){
  ab_score <- mean(x[occ==1])
  return(ab_score)
}
jaccard_list <- list()
percent_list <- list()
abundance_list <- list()
#pdf(family="Helvetica",paste(type,'_vs_provinces.pdf', sep=''), width = 15, height = 15)
for (fra in fracs){
  df_frac <- data_prov[data_prov$Fraction==fra,]
  dsm <- as.numeric(as.matrix(data_SMAGs0[2:(nmag+1),data_SMAGs0[1,]==fra]))
  
  stats <- stat_SMAGs_prov[data_SMAGs0[1,]==fra]
  dsm <-matrix(dsm, nrow = nmag, ncol=length(dsm)/nmag, byrow = F)
  # if (type=='SMAGs'){
  dsm_ab <- as.numeric(as.matrix(perc_smags0[1:nmag,data_SMAGs0[1,]==fra]))
  dsm_ab <- matrix(dsm_ab, nrow = nmag, ncol=length(dsm)/nmag, byrow = F)
  colnames(dsm_ab)<- stats
  # }
  
  colnames(dsm) <- stats
  
  dsm <- dsm[,colnames(dsm) %in% df_frac$Station]
  
  df_frac <- df_frac[df_frac$Station %in% stats,]
  sums <- apply(dsm, 1, sum)
  dsm <- dsm[sums>4,]
  to_take <- match(df_frac$Station,colnames(dsm))
  to_take <- to_take[!is.na(to_take)]
  dsm0 <- dsm[,to_take]
  Smgs <- smags[sums>4]
  let <- letters[fracs==fra]
  # if (type=='SMAGs'){
  dsm_ab <- dsm_ab[,colnames(dsm_ab) %in% df_frac$Station]
  dsm_ab <- dsm_ab[sums>4,]
  dsm0_ab <- dsm_ab[,to_take]
  # }
  jacc_vecs <- NULL
  perc_vecs <- NULL
  ab_vecs <- NULL
  labs <- NULL
  for (gp in unique(df_frac$lab)){
    occ_vec <- as.numeric(df_frac$lab==gp)
    num <- df_frac$Genocenose[df_frac$lab==gp][1]
    jacc_vec <- apply(dsm0, 1,FUN =  jaccard_index, occ=occ_vec)
    perc_vec <- apply(dsm0, 1,FUN =  perc_stations, occ=occ_vec)
    # if (type=='SMAGs'){
    ab_vec <- apply(dsm0_ab, 1,FUN =  ab_stations, occ=occ_vec)
    # }
    jacc_vecs <- cbind(jacc_vecs, jacc_vec)
    perc_vecs <- cbind(perc_vecs, perc_vec)
    # if (type=='SMAGs'){
    ab_vecs <- cbind(ab_vecs, ab_vec)
    # }
    labs <- append(labs, paste(let, num, sep=''))
  }
 
  colnames(jacc_vecs) <- labs
  rownames(jacc_vecs) <- Smgs
  colnames(perc_vecs) <- labs
  rownames(perc_vecs) <- Smgs
  # if (type=='SMAGs'){
  colnames(ab_vecs) <- labs
  rownames(ab_vecs) <- Smgs
  # }
  if (type=='SMAGs'){
    br <- seq(0,1, length.out = 100)
  } else if (type=='MAGprok'){
    br <- seq(0,1, length.out = 100)
  }
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  lwid = c(1.5,4)
  lhei = c(1.5,4,1)
  jaccard_list[[fra]]<- jacc_vecs
  percent_list[[fra]]<- perc_vecs
  # if (type=='SMAGs'){
  abundance_list[[fra]]<- ab_vecs
  # }
  if ( (fra == '0.8-5' & type=='MAGprok') | (fra == '0.22-3' & type=='MAGprok')){
    cx <- 0.1
  } else {
    cx <- 0.2
  }
#   heatmap.2(jacc_vecs, col= jet.colors(99) ,breaks = br, trace="none", margins=c(10,10), cexRow=cx,
#             lmat = lmat, lwid = lwid, lhei = lhei, density.info = 'none')
}
#dev.off()

if (type=='MAGprok'){
  tax <- read.table(paste(type, '_taxo_1.txt',sep=''), sep='\t', header = T)
  tax$N_Fixation_Phylum <- as.character(levels(tax$N_Fixation_Phylum))[tax$N_Fixation_Phylum]
} else{
  #tax <- read.table(paste(type, '_16102020.txt', sep=''), header = T, sep = '\t')
  tax <- readRDS('SMAGs_Statistics_02.rds')
  tax <- tax[!grepl('METDB',tax$Genome_Id),]
#   test <- tax[,c(11, 72:75)]
#   keep_max <- function(x){
#     ind <- which.max(x)
#     x[c(1:5)!=ind]=0
#     return(x)
#   }
#   test0 <- t(apply(test, 1, keep_max))
#   test0 <- as.data.frame(test0)
#   tax$Photosynthesis_prediction <- test0$Photosynthesis_prediction
#   tax$Prototrophy_prediction <- test0$Prototrophy_prediction
#   tax$`Phagocyte-rozellid_prediction` <- test0$`Phagocyte-rozellid_prediction`
#   tax$`Phagocyte-generalist_prediction` <- test0$`Phagocyte-generalist_prediction`
#   tax$`Phagocyte-entamoebid_prediction` <- test0$`Phagocyte-entamoebid_prediction`
    tax$Photosynthesis_prediction <- tax$Photosynthesis_prediction/4
    tax$Prototrophy_prediction <- tax$Prototrophy_prediction/4  
    tax$`Phagocyte-generalist_prediction` <- tax$`Phagocyte-generalist_prediction`/6
    tax$`Phagocyte-rozellid_prediction` <- tax$`Phagocyte-rozellid_prediction`/6
    tax$`Phagocyte-entamoebid_prediction` <- tax$`Phagocyte-entamoebid_prediction`/6
}



extract_funct<-function(weighter, funct_list, type0){
  compositional_functions <- rep(list(NULL), length(fracs))
#   compositional_functions_0 <- rep(list(NULL), length(fracs))
#   compositional_functions_1 <- rep(list(NULL), length(fracs))
  names(compositional_functions)<- fracs
#   names(compositional_functions_0)<- fracs
#   names(compositional_functions_1)<- fracs
  
  for (frac in fracs){
    jacc_ <- weighter[[frac]]
    compositional_functions[[frac]]<- rep(list(NULL), dim(jacc_)[2])
    names(compositional_functions[[frac]]) <- colnames(jacc_)
#     compositional_functions_0[[frac]]<- rep(list(NULL), dim(jacc_)[2])
#     names(compositional_functions_0[[frac]]) <- colnames(jacc_)
#     compositional_functions_1[[frac]]<- rep(list(NULL), dim(jacc_)[2])
#     names(compositional_functions_1[[frac]]) <- colnames(jacc_)
    for (pr in 1:dim(jacc_)[2]){
      parts <- NULL
      for (gr in funct_list){
        
        if (type0=='genomic'){
          find_gr <- function(x, gr){
            if (grepl('TOSAG', x)){
              x <- strsplit(x, split = '_scaffold')[[1]][1]
            }
            v <- tax$Functional_group[tax$Genome_Id==x]==gr
            return(v)
          }
          
          sel <- sapply(rownames(jacc_), FUN = find_gr, gr=gr)
          part <- sum(jacc_[sel,pr])
          parts <- append(parts, part)
        } else if (type0=='trophic'){
          find_fu <- function(x, fu){
            if (grepl('TOSAG', x)){
              x <- strsplit(x, split = '_scaffold')[[1]][1]
            }
            s <- tax[[fu]][tax$Genome_Id==x]
            return(s)
          }
          scores <- sapply(rownames(jacc_), FUN = find_fu, fu=gr)
          part <- sum(scores*jacc_[,pr])
          parts <- append(parts, part)
        } else if (type0=='PFT'){
          find_gr1 <- function(x, gr){
            if (grepl('TOSAG', x)){
              x <- strsplit(x, split = '_scaffold')[[1]][1]
            }
            v <- tax$PFT[tax$Genome_Id==x]==gr
            return(v)
          }
          sel1 <- sapply(rownames(jacc_), FUN = find_gr1, gr=gr)
          part <- sum(jacc_[sel1,pr])
          parts <- append(parts, part)
        } else if (type0=='Diaz'){
          find_gr1 <- function(x, gr){
            if (grepl('TOSAG', x)){
              x <- strsplit(x, split = '_scaffold')[[1]][1]
            }
            v <- tax$N_Fixation_Phylum[tax$Genome_Id==x]==gr
            return(v)
          }
          sel1 <- sapply(rownames(jacc_), FUN = find_gr1, gr=gr)
          part <- sum(jacc_[sel1,pr])
          parts <- append(parts, part)
        }  else if (type0 %in% colnames(tax)){
          find_gr1 <- function(x, gr){
            if (grepl('TOSAG', x)){
              x <- strsplit(x, split = '_scaffold')[[1]][1]
            }
            v <- tax[[type0]][tax$Genome_Id==x]==gr
            return(v)
          }
          sel1 <- sapply(rownames(jacc_), FUN = find_gr1, gr=gr)
          part <- sum(jacc_[sel1,pr])
          parts <- append(parts, part)
        }
      }
      
        
      
      # colnames(parts) <- unique(tax$Functional_group)
      compositional_functions[[frac]][[colnames(jacc_)[pr]]] <- parts
#       compositional_functions_0[[frac]][[colnames(jacc_)[pr]]] <- parts_0
#       compositional_functions_1[[frac]][[colnames(jacc_)[pr]]] <- parts_1
    }
  }
  
  return(compositional_functions)
}
funct_list<- c('Photosynthesis_prediction', 'Prototrophy_prediction',
                 'Phagocyte-generalist_prediction', 'Phagocyte-rozellid_prediction',
                 'Phagocyte-entamoebid_prediction')
if (type=='SMAGs'){
  PFT_list <- unique(tax$PFT)
  Group_list <- unique(tax$Functional_group)
  Hex_list <- c("Marine_Hexanauplia_A", "Marine_Hexanauplia_B")
  Hex_list_bis <- unique(tax$Hex_bis)
  Hex_list_bis <- Hex_list_bis[!is.na(Hex_list_bis) & Hex_list_bis!='0']
  Hex_list_bis <- Hex_list_bis[c(4, 3, 2, 1, 5)]
}
if (type=='MAGprok'){
  Diaz_list <- unique(tax$N_Fixation_Phylum)
  Diaz_list <- Diaz_list[Diaz_list!=""]
}
if (type=='SMAGs'){
#   v1 <- extract_funct(jaccard_list, funct_list)
#   v2 <- extract_funct(percent_list, funct_list)
#   v3 <- extract_funct(abundance_list, funct_list)
  compositional_functions <- extract_funct(jaccard_list, Group_list, 'genomic')
  compositional_functions_0 <- extract_funct(jaccard_list, funct_list , 'trophic')
  compositional_functions_1 <- extract_funct(jaccard_list, PFT_list, 'PFT')
  compositional_functions_hex <- extract_funct(jaccard_list, Hex_list, 'Hex')
  compositional_functions_hex_bis <- extract_funct(jaccard_list, Hex_list_bis, 'Hex_bis')
  compositional_functions_p <- extract_funct(percent_list, Group_list, 'genomic')
  compositional_functions_0_p <- extract_funct(percent_list, funct_list, 'trophic')
  compositional_functions_1_p <- extract_funct(percent_list, PFT_list, 'PFT')
  compositional_functions_hex_p <- extract_funct(percent_list, Hex_list, 'Hex')
  compositional_functions_hex_p_bis <- extract_funct(percent_list, Hex_list_bis, 'Hex_bis')
  compositional_functions_a <- extract_funct(abundance_list, Group_list, 'genomic')
  compositional_functions_0_a <- extract_funct(abundance_list, funct_list, 'trophic')
  compositional_functions_1_a <- extract_funct(abundance_list, PFT_list, 'PFT')
  compositional_functions_hex_a <- extract_funct(abundance_list, Hex_list, 'Hex')
  compositional_functions_hex_a_bis <- extract_funct(abundance_list, Hex_list_bis, 'Hex_bis')
  phylums <- unique(tax$Best_taxonomy_PHYLUM)
  compositional_functions_ph <- extract_funct(jaccard_list, phylums, 'Best_taxonomy_PHYLUM')
  compositional_functions_p_ph <- extract_funct(percent_list, phylums, 'Best_taxonomy_PHYLUM')
  compositional_functions_a_ph <- extract_funct(abundance_list, phylums, 'Best_taxonomy_PHYLUM')
  

  func_groups <- unique(tax$Functional_group)
  
  piechart <- function(slices, p, classif, name, name0){
    sel <- slices!=0
    lbls <- classif[sel]
    pct <- round(slices[sel]/sum(slices[sel])*100)
    lbls <- paste(lbls, pct) # add percents to labels
    lbls <- paste(lbls,"%",sep="") # ad % to labels
    lbls <- paste(lbls, '\n','score= ' ,round(slices[sel], 1), sep='')
    if (name0=='Hexanauplia'){
      colos <- c('darkred', 'midnightblue')
    } else if (name0=='Hexanauplia_bis'){
      colos=c('darkred','red','lightpink' ,'lightblue','blue')[sel]
    }else{
      colos=rainbow(length(classif))[sel]
    }
    pie(slices[sel],labels = lbls, col=colos,
        main=paste("Pie Chart of ",name0," functional groups (",name,"): province ", p, sep=''))
  }
  
  for( frac in fracs){
    pdf(paste('functionnal_composition_SMAGs_provinces_', frac, '.pdf', sep=''), width=22.5, 
        height = 10)
    par(mfrow=c(3,3), mar=c(2.1, 5.1, 4.1, 5.1))
    for (na in names(compositional_functions[[frac]])){
      piechart(compositional_functions[[frac]][[na]], na, func_groups, 
               name='weight: Jaccard index', name0='genomic')
      piechart(compositional_functions_1[[frac]][[na]], na, PFT_list,
               name='weight: Jaccard index', name0='PFT')
      piechart(compositional_functions_0[[frac]][[na]], na, funct_list,
               name='weight: Jaccard index', name0='trophic')
      
      piechart(compositional_functions_p[[frac]][[na]], na, func_groups,
               name='weight: %', name0='genomic')
      piechart(compositional_functions_1_p[[frac]][[na]], na, PFT_list,
               name='weight: %', name0='PFT')
      piechart(compositional_functions_0_p[[frac]][[na]], na, funct_list,
               name='weight: %', name0='trophic')
      
      piechart(compositional_functions_a[[frac]][[na]], na, func_groups,
               name='weight: mean abundance', name0='genomic')
      piechart(compositional_functions_1_a[[frac]][[na]], na, PFT_list,
               name='weight: mean abundance', name0='PFT')
      piechart(compositional_functions_0_a[[frac]][[na]], na, funct_list,
               name='weight: mean abundance', name0='trophic')
    }
    dev.off()
    pdf(paste('functionnal_composition_phylums_SMAGs_provinces_', frac, '.pdf', sep=''), width=7.5, 
        height = 10)
    par(mfrow=c(3,1), mar=c(2.1, 5.1, 4.1, 5.1))
    for (na in names(compositional_functions[[frac]])){
      piechart(compositional_functions_ph[[frac]][[na]], na, phylums, 
               name='weight: Jaccard index', name0='genomic')
      piechart(compositional_functions_p_ph[[frac]][[na]], na, phylums,
               name='weight:  %', name0='PFT')
      piechart(compositional_functions_a_ph[[frac]][[na]], na, phylums,
               name='weight: Abundance', name0='trophic')
    }
    dev.off()
    
    if (frac %in% c("180-2000", "20-180" ,  "5-20")){
      pdf(paste('functionnal_composition_hexanaupia_SMAGs_provinces_', frac, '.pdf', sep=''), width=7.5, 
          height = 10)
      par(mfrow=c(3,1), mar=c(2.1, 5.1, 4.1, 5.1))
      for (na in names(compositional_functions_hex[[frac]])){
        piechart(compositional_functions_hex[[frac]][[na]], na, Hex_list, 
                 name='weight: Jaccard index', name0='Hexanauplia')
        piechart(compositional_functions_hex_p[[frac]][[na]], na, Hex_list,
                 name='weight:  %', name0='Hexanauplia')
        piechart(compositional_functions_hex_a[[frac]][[na]], na, Hex_list,
                 name='weight: Abundance', name0='Hexanauplia')
      }
      dev.off()
      pdf(paste('functionnal_composition_hexanauplia_bis_SMAGs_provinces_', frac, '.pdf', sep=''), width=7.5, 
          height = 10)
      par(mfrow=c(3,1), mar=c(2.1, 5.1, 4.1, 5.1))
      for (na in names(compositional_functions_hex[[frac]])){
        piechart(compositional_functions_hex_bis[[frac]][[na]], na, Hex_list_bis, 
                 name='weight: Jaccard index', name0='Hexanauplia_bis')
        piechart(compositional_functions_hex_p_bis[[frac]][[na]], na, Hex_list_bis,
                 name='weight:  %', name0='Hexanauplia_bis')
        piechart(compositional_functions_hex_a_bis[[frac]][[na]], na, Hex_list_bis,
                 name='weight: Abundance', name0='Hexanauplia_bis')
      }
      dev.off()
    }
  }
  comp_all <- list(compositional_functions, compositional_functions_0,compositional_functions_1,compositional_functions_hex,compositional_functions_hex_bis,
                   compositional_functions_p, compositional_functions_0_p, compositional_functions_1_p,compositional_functions_hex_p,compositional_functions_hex_p_bis,
                   compositional_functions_a, compositional_functions_0_a, compositional_functions_1_a,compositional_functions_hex_a,compositional_functions_hex_a_bis,
                   compositional_functions_ph, compositional_functions_p_ph, compositional_functions_a_ph)
 names(comp_all) <- c('Jaccard genomics', 'Jaccard trophic','Jaccard PFT','Jaccard Hex','Jaccard Hex_bis',
                      'Perc genomics', 'Perc trophic', 'Perc PFT', 'Perc Hex','Perc Hex_bis',
                      'Abundance genomics', 'Abundance trophic', 'Abundance PFT','Abundance Hex','Abundance Hex_bis',
                      'Jaccard Phylum', 
                      'Perc Phylum',
                      'Abundance Phylum') 
 saveRDS(comp_all,'SMAGs_provinces_compositional_functions_genomics_trophic.rds')
} else if (type=='MAGprok'){
  compositional_functions <- extract_funct(jaccard_list, Diaz_list, 'Diaz')
  compositional_functions_p <- extract_funct(percent_list, Diaz_list, 'Diaz')
  compositional_functions_a <- extract_funct(abundance_list, Diaz_list, 'Diaz')
  
  phylums <- unique(tax$Phylum.gtdbtk)
  compositional_functions_ph <- extract_funct(jaccard_list, phylums, 'Phylum.gtdbtk')
  compositional_functions_p_ph <- extract_funct(percent_list, phylums, 'Phylum.gtdbtk')
  compositional_functions_a_ph <- extract_funct(abundance_list, phylums, 'Phylum.gtdbtk')

  
  piechart <- function(slices, p, classif, name, name0){
    sel <- slices!=0
    lbls <- classif[sel]
    pct <- round(slices[sel]/sum(slices[sel])*100)
    lbls <- paste(lbls, pct) # add percents to labels
    lbls <- paste(lbls,"%",sep="") # ad % to labels
    lbls <- paste(lbls, '\n','score= ' ,round(slices[sel], 3), sep='')
    pie(slices[sel],labels = lbls, col=rainbow(length(classif))[sel],
        main=paste("Pie Chart of ",name0," functional groups (",name,"): province ", p, sep=''))
  }
  
  for( frac in fracs){
    pdf(paste('functionnal_composition_MAGprok_provinces_', frac, '.pdf', sep=''), width=7.5, 
        height = 10)
    par(mfrow=c(3,1), mar=c(2.1, 5.1, 4.1, 5.1))
    for (na in names(compositional_functions[[frac]])){
      if (sum(compositional_functions[[frac]][[na]])!=0){
        piechart(compositional_functions[[frac]][[na]], na, Diaz_list, 
                 name='weight: Jaccard index', name0='Diazotrophy')
        
        piechart(compositional_functions_p[[frac]][[na]], na, Diaz_list,
                 name='weight: %', name0='Diazotrophy')
        
        piechart(compositional_functions_a[[frac]][[na]], na,Diaz_list,
                 name='weight: mean abundance', name0='Diazotrophy')
      }
    }
    dev.off()
    pdf(paste('functionnal_composition_phylum_MAGprok_provinces_', frac, '.pdf', sep=''), width=7.5, 
        height = 10)
    par(mfrow=c(3,1), mar=c(2.1, 5.1, 4.1, 5.1))
    for (na in names(compositional_functions[[frac]])){
      if (sum(compositional_functions[[frac]][[na]])!=0){
        piechart(compositional_functions_a_ph[[frac]][[na]], na, phylums, 
                 name='weight: Jaccard index', name0='Phylum')
        
        piechart(compositional_functions_p_ph[[frac]][[na]], na, phylums,
                 name='weight: %', name0='Phylum')
        
        piechart(compositional_functions_a_ph[[frac]][[na]], na, phylums,
                 name='weight: mean abundance', name0='Phylum')
      }
    }
    dev.off()
  }
  comp_all <- list(compositional_functions, 
                   compositional_functions_p, 
                   compositional_functions_a, 
                   compositional_functions_ph, 
                   compositional_functions_p_ph, 
                   compositional_functions_a_ph)
  names(comp_all) <- c('Jaccard Diazotrophy', 
                       'Perc Diazotrophy',
                       'Abundance Diazotrophy',
                       'Jaccard Phylum', 
                       'Perc Phylum',
                       'Abundance Phylum'
                       ) 
  saveRDS(comp_all,'MAGprok_provinces_compositional_functions_diazotrophy.rds')
}

signature_MAGs <- rep(list(NULL),length(fracs))
signature_MAGs_m <- rep(list(NULL),length(fracs))
signature_MAGs_func <- rep(list(NULL),length(fracs))
signature_MAGs_fr <- rep(list(NULL),length(fracs))
names(signature_MAGs)<-fracs
if (type =='MAGprok'){
  for(fra in fracs){
    jc <- jaccard_list[[fra]]
    groups <- colnames(jc)
    signatures <- rep(list(NULL), length(groups))
    signatures_m <- rep(list(NULL), length(groups))
    signatures_f <- rep(list(NULL), length(groups))
    signatures_all <- NULL
    names(signatures) <- groups
    for (j in 1:dim(jc)[1]){
      ids <- 1:length(groups)
      for (i in ids){
        if (jc[j,i]>0.5){
          ok <-T
          ids_e <- ids[ids!=i]
          for (k in ids_e){
            if (jc[j,k]>0.1){
              ok <- F
            }
          }
        } else{
          ok <-F
        }
        if (ok==T){
          if (version=='old'){
            mg <- rownames(jc)[j]
            taxid <- as.character(tax[tax$user_genome==mg,2 ])
            taxid_f <- taxid
            for (kl in 3:8){
              tx <- as.character(tax[tax$user_genome==mg,kl ])
              if ((is.na(tx)) || tx==taxid ){
                taxid_f <- taxid_f
                break
              } else{
                taxid_f <- paste(taxid_f,'_', tx, sep='')
              }
              taxid <-tx
            }
            mg <- paste(mg,taxid_f, sep='_')
            signatures_m[[ groups[i] ]] <- append(signatures_m[[ groups[i] ]], rownames(jc)[j])
            signatures[[ groups[i] ]] <- append(signatures[[ groups[i] ]], mg)
            signatures_all <-append(signatures_all, rownames(jc)[j])
          } else{
            mg <- rownames(jc)[j]
            taxid <- as.character(tax[tax$Genome_Id==mg,30 ])
            taxid_f <- taxid
            for (kl in 31:35){
              tx <- as.character(tax[tax$Genome_Id==mg,kl ])
              if ((is.na(tx)) || tx==taxid ){
                taxid_f <- taxid_f
                break
              } else{
                taxid_f <- paste(taxid_f,'_', tx, sep='')
              }
              taxid <-tx
            }
            mg <- paste(mg,taxid_f, sep='_')
            signatures_m[[ groups[i] ]] <- append(signatures_m[[ groups[i] ]], rownames(jc)[j])
            signatures[[ groups[i] ]] <- append(signatures[[ groups[i] ]], mg)
            signatures_all <-append(signatures_all, rownames(jc)[j])
          }
        }
      }
    }
    signature_MAGs_m[[fra]] <-signatures_m
    signature_MAGs[[fra]] <-signatures
    signature_MAGs_fr[[fra]] <-signatures_all
    
  }
#   for (fra in fracs){
#     jacc_veccs <- jaccard_list[[fra]]
#     for (group in names(signature_MAGs[[fra]])){
#       write(group, paste(group, '_signature_SMAGs.txt', sep=''), append = T)
#       write(signature_MAGs[[fra]][[group]], paste(group, '_signature_SMAGs.txt', sep=''), append = T)
#     }
#   }
  pdf(family="Helvetica",paste(type,'_',version,'_vs_provinces.pdf', sep=''), width = 15, height = 15)
  for (fra in fracs){
    jacc_veccs <- jaccard_list[[fra]]
#     for (group in names(signature_MAGs[[fra]])){
#       write(group, paste(group,'_',version, '_signature_MAGprok.txt', sep=''), append = T)
#       write(signature_MAGs[[fra]][[group]], paste(group,'_',version, '_signature_MAGprok.txt', sep=''), append = T)
#       #write(group, paste(group,'_',version, '_signature_names_MAGprok.txt', sep=''), append = T)
#       write(signature_MAGs_m[[fra]][[group]], paste(group,'_',version, '_signature_names_MAGprok.txt', sep=''), append = T)
#     }
    
    if (type=='SMAGs'){
      br <- seq(0,1, length.out = 100)
    } else if (type=='MAGprok'){
      br <- seq(0,1, length.out = 100)
    }
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    lwid = c(1.5,4)
    lhei = c(1.5,4,1)
    
    if ( (fra == '0.8-5' & type=='MAGprok') | (fra == '0.22-3' & type=='MAGprok')){
      cx <- 0.2
    } else {
      cx <- 0.3
    }
    condi <- !(rownames(jacc_veccs) %in% signature_MAGs_fr[[fra]])
    rownames(jacc_veccs)[condi] <- NA
    jacc_veccs <- jacc_veccs[,order(colnames(jacc_veccs))]
    heatmap.2(jacc_veccs, col= jet.colors(99) ,breaks = br, trace="none", margins=c(10,10), cexRow=cx,
              lmat = lmat, lwid = lwid, lhei = lhei, density.info = 'none', Colv = F)
  }
  dev.off()
} else if (type =='SMAGs'){
  for(fra in fracs){
    jc <- jaccard_list[[fra]]
    groups <- colnames(jc)
    signatures <- rep(list(NULL), length(groups))
    signatures_f <- rep(list(NULL), length(groups))
    signatures_m <- rep(list(NULL), length(groups))
    signatures_all <- NULL
    names(signatures) <- groups
    for (j in 1:dim(jc)[1]){
      ids <- 1:length(groups)
      for (i in ids){
        if (jc[j,i]>0.5){
          ok <-T
          ids_e <- ids[ids!=i]
          for (k in ids_e){
            if (jc[j,k]>0.1){
              ok <- F
            }
          }
        } else{
          ok <-F
        }
        if (ok==T){
          mg <- rownames(jc)[j]
          if (grepl('scaffold', mg)){
            mg <- strsplit(mg, '_scaffold')[[1]][1]
          }
          taxid <- as.character(tax[tax$Genome_Id==mg,14 ])
          taxid_f <- taxid
          for (kl in 15:20){
            tx <- as.character(tax[tax$Genome_Id==mg,kl ])
            if (is.na(tx) || tx==taxid ){
              taxid_f <- taxid_f
              break
            } else{
              taxid_f <- paste(taxid_f,'_', tx, sep='')
            }
            taxid <-tx
          }
          selec <- c(1, 9:12, 72:75, 15:20)
          to_save <- tax[tax$Genome_Id==mg, selec]
          mg <- paste(mg,taxid_f, sep='_')
          mag <- rownames(jc)[j]
          if (grepl('scaffold', mag)){
            mag <- strsplit(mag, split = '_scaffold')[[1]]
          }
          signatures_m[[ groups[i] ]] <- append(signatures_m[[ groups[i] ]], mag)
          signatures_f[[ groups[i] ]] <- rbind(signatures_f[[ groups[i] ]], to_save)
          signatures[[ groups[i] ]] <- append(signatures[[ groups[i] ]], mg)
          signatures_all <- append(signatures_all,rownames(jc)[j] )
        }
      }
    }
    signature_MAGs_m[[fra]] <- signatures_m
    signature_MAGs[[fra]] <- signatures
    signature_MAGs_func[[fra]] <- signatures_f
    signature_MAGs_fr[[fra]] <- signatures_all
  }
#   for (fra in fracs){
#     jacc_veccs <- jaccard_list[[fra]]
#     for (group in names(signature_MAGs[[fra]])){
#       write(group, paste(group, '_signature_SMAGs.txt', sep=''), append = T)
#       write(signature_MAGs[[fra]][[group]], paste(group, '_signature_SMAGs.txt', sep=''), append = T)
#     }
#   }
  pdf(family="Helvetica",paste(type,'_vs_provinces.pdf', sep=''), width = 15, height = 15)
  for (fra in fracs){
    jacc_veccs <- jaccard_list[[fra]]
#     for (group in names(signature_MAGs[[fra]])){
#       write(group, paste(group, '_signature_SMAGs.txt', sep=''), append = T)
#       write(signature_MAGs[[fra]][[group]], paste(group, '_signature_SMAGs.txt', sep=''), append = T)
#       #write(group, paste(group,'_',version, '_signature_names_SMAGs.txt', sep=''), append = T)
#       write(signature_MAGs_m[[fra]][[group]], paste(group,'_',version, '_signature_names_SMAGs.txt', sep=''), append = T)
#     }
    
    if (type=='SMAGs'){
      br <- seq(0,1, length.out = 100)
    } else if (type=='MAGprok'){
      br <- seq(0,1, length.out = 100)
    }
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    lwid = c(1.5,4)
    lhei = c(1.5,4,1)
    
    if ( (fra == '0.8-5' & type=='MAGprok') | (fra == '0.22-3' & type=='MAGprok')){
      cx <- 0.2
    } else {
      cx <- 0.3
    }
    condi <- !(rownames(jacc_veccs) %in% signature_MAGs_fr[[fra]])
    rownames(jacc_veccs)[condi] <- NA
    jacc_veccs <- jacc_veccs[,order(colnames(jacc_veccs))]
    heatmap.2(jacc_veccs, col= jet.colors(99) ,breaks = br, trace="none", margins=c(10,10), cexRow=cx,
              lmat = lmat, lwid = lwid, lhei = lhei, density.info = 'none', Colv = F)
  }
  dev.off()
}

tax$province_sign <- NA
for (fra in fracs){
  for (group in names(signature_MAGs_m[[fra]])){
    for (g in signature_MAGs_m[[fra]][[group]]){
      if (!is.na(tax$province_sign[tax$Genome_Id==g])){
        tax$province_sign[tax$Genome_Id==g]=paste(tax$province_sign[tax$Genome_Id==g],group, sep='_')
      } else{
        tax$province_sign[tax$Genome_Id==g]=group
      }
    }
  }
}

if (type='SMAGs'){
  to_s <- cbind(tax$Genome_Id, tax$province_sign)
} else{
  to_s <- cbind(as.character(levels(tax$Genome_Id))[tax$Genome_Id], tax$province_sign)
}
write.table(to_s, paste(type, '_', version, '_signatures.txt', sep=''), row.names = F)





