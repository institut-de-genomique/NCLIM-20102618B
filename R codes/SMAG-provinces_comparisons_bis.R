library('gplots')
#library('readxl')
library('matlab')

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

data_prov <- readRDS('Provinces_env_parameters_woa_scaled.rds')
data_prov$Fraction[data_prov$Fraction=='43952'] <- '5-20'
data_prov$lab <- paste(data_prov$Fraction, data_prov$Province, sep='_')
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
if (type=='MAGprok'){
  tax <- read.table(paste(type, '_taxo_1.txt',sep=''), sep='\t', header = T)
  tax$N_Fixation_Phylum <- as.character(levels(tax$N_Fixation_Phylum))[tax$N_Fixation_Phylum]
} else{
  tax <- readRDS('SMAGs_Statistics_02.rds')
  tax <- tax[!grepl('METDB',tax$Genome_Id),]
  tax$Photosynthesis_prediction <- tax$Photosynthesis_prediction/4
  tax$Prototrophy_prediction <- tax$Prototrophy_prediction/4  
  tax$`Phagocyte-generalist_prediction` <- tax$`Phagocyte-generalist_prediction`/6
  tax$`Phagocyte-rozellid_prediction` <- tax$`Phagocyte-rozellid_prediction`/6
  tax$`Phagocyte-entamoebid_prediction` <- tax$`Phagocyte-entamoebid_prediction`/6
}

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
letters <- c('A', 'B', 'C', 'D', 'E', 'F')
fractions <- c('0-0.2', '0.22-3', '0.8-5', '5-20', '20-180', '180-2000')
data_prov$prov <- paste(letters[match(data_prov$Fraction, fractions)], data_prov$Province, sep='')

tests_func <- function(fracs, grouping, name){
  final_data <- NULL
  data_ <- rep(list(NULL), length(grouping))
  names(data_)<-grouping
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
    mgs <- perc_smags[sums>4,1]
    ty_s <- function(x){
      return(tax[[name]][ tax[['Genome_Id']]==x ])
    }
    ty_mags <- sapply(mgs, FUN = ty_s )
    dsm0_ab <- dsm_ab[,to_take]
    data <- rep(list(NULL), length(grouping))
    names(data) <- grouping
    for (ty in grouping){
      mgs0 <- mgs[ty_mags==ty]
      data_[[ty]] <- append(data_[[ty]], mgs0)
      dsm0_ab_ty <- dsm0_ab[ty_mags==ty,]
      for (gp in unique(df_frac$prov)){
        occ_vec <- as.numeric(df_frac$prov==gp)
        num <- df_frac$Province[df_frac$prov==gp][1]
        if (!is.null(dim(dsm0_ab_ty))){
          dsm0_ab_ty_p <- dsm0_ab_ty[,occ_vec==1]
          fin_vec <- apply(dsm0_ab_ty_p, 2, sum)
        } else{
          dsm0_ab_ty_p <- dsm0_ab_ty[occ_vec==1]
          fin_vec <- dsm0_ab_ty_p
        }
        data[[ty]] <- rbind(data[[ty]], cbind(fin_vec, rep(gp, length(fin_vec))) )
      }
      data[[ty]] <- as.data.frame(data[[ty]])
      colnames(data[[ty]]) <- c('Ab', 'Prov')
      data[[ty]]$Ab <- as.numeric(levels(data[[ty]]$Ab))[data[[ty]]$Ab]
      v <- pairwise.wilcox.test(data[[ty]]$Ab, data[[ty]]$Prov, paired = F)$p.value
#       print(ty)
#       print(v)
#       print('')
      transi <- NULL
      for (j in colnames(v)){
        transi <- c(transi, paste(rownames(v), j, sep='->'))
      }
      for (j in colnames(v)){
        transi <- c(transi, paste(j, rownames(v), sep='->'))
      }
      v <- as.vector(v)
      final_data <- rbind(final_data, cbind(transi, c(v, v), rep(ty,length(v)*2)))
    }
    
  }
  for (ty in grouping){
    h <- unique(data_[[ty]])
    print(ty)
    print(length(h))
  }
  colnames(final_data)<- c('Trans', 'p-val', 'Clade')
  return(final_data)
}

if (type=='SMAGs'){
  data<-tests_func(fracs = c('180-2000', '20-180', '5-20'), grouping = Hex_list_bis, name = 'Hex_bis')
  data <- data[as.numeric(data[,2])<0.05 & !is.na(data[,2]),]
  write.table(data, 'significant_differences_hexanauplia_wilcox.txt', row.names = F)
} else{
  data<-tests_func(fracs = fractions[2:6], grouping = Diaz_list[c(1,3,4, 5,6, 8)], name = 'N_Fixation_Phylum')
  cond <- as.numeric(data[,2])<0.05 & !is.na(data[,2])
  cond[is.na(cond)]<-F
  data <- data[cond,]
  write.table(data, 'significant_differences_diazotrophy_wilcox.txt', row.names = F)
}

