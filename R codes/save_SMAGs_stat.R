²# library('readxl')
library('car')
type <-'SMAGs'
# tax <- read_xlsx(paste(type, '_Statistics_02.xlsx', sep=''))
tax <- readRDS('SMAGs_Statistics_02.rds')
tax <- tax[!grepl('METDB',tax$Genome_Id),]
tax$Functional_group[tax$Functional_group=='Group_01'] <- 'Group A'
tax$Functional_group[tax$Functional_group=='Group_02'] <- 'Group B'
tax$Functional_group[tax$Functional_group=='Group_03'] <- 'Group C'
tax$Functional_group[tax$Functional_group=='Group_04'] <- 'Group D'
tax$Functional_group[tax$Functional_group=='Group_05'] <- 'Group E'
tax$Functional_group[is.na(tax$Functional_group)] <- 'No group'
tax$Photosynthesis_prediction[is.na(tax$Photosynthesis_prediction)] <- 0
tax$Prototrophy_prediction[is.na(tax$Prototrophy_prediction)] <- 0
tax$`Phagocyte-generalist_prediction`[is.na(tax$`Phagocyte-generalist_prediction`)] <- 0
tax$`Phagocyte-rozellid_prediction`[is.na(tax$`Phagocyte-rozellid_prediction`)] <- 0
tax$`Phagocyte-entamoebid_prediction`[is.na(tax$`Phagocyte-entamoebid_prediction`)] <- 0
tax$PFT <- NA
tax$size_hexanauplia <- NA
tax$pval_size_hexanauplia <- NA
tax$lev_test <- NA
perc_smags <- read.table('SMAGs_percent_mapped_read.txt', header = T)
perc_magprok <- read.table('MAGprok_percent_mapped_read.txt', header=T)

tax$size_phyto <- NA
tax$pval_size_phy <- NA
tax$lev_test_phy <- NA
tax$size_phyto_bis <- NA
tax$pval_size_phy_bis <- NA
tax$lev_test_phy_bis <- NA
i=1
for (Smg in tax$Genome_Id){
  
  if (!is.na(tax$Phytoplankton[tax$Genome_Id==Smg])){
    if (tax$Best_taxonomy_PHYLUM[tax$Genome_Id==Smg]=="Bacillariophyta"){
      tax$PFT[i] <- 'Diatoms'
    } else{
      tax$PFT[i] <- 'Nanophytoplankton'
    }
  }
  
  if (tax$Functional_group[tax$Genome_Id==Smg] %in% c('Small_Animals', 'Group E', 'Outlier') ){
    cond_big <- grepl('SSUU', colnames(perc_smags))
    cond_small <- grepl('QQSS', colnames(perc_smags)) 
    cond_small1 <- grepl('KKQQ', colnames(perc_smags)) | grepl('MMQQ', colnames(perc_smags))
    mean_big_frac <- mean(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_big]))
    mean_small_frac <- mean(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small]))
    mean_small_frac1 <- mean(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small1]))
    mean_small <- mean(c(mean_small_frac,mean_small_frac1 ))
    p_val <- wilcox.test(c(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small])), 
                    as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_big]))
    
    bigs <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_big])
    datab <- cbind(rep('SSUU', length(bigs)), bigs)
    meds <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small])
    datam <- cbind(rep('QSSS', length(meds)), meds)
    smlls <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small1])
    datas <- cbind(rep('KKQQ', length(smlls)), smlls)
    data <- rbind(datab,datam, datas )
    colnames(data) <- c('Group', 'Ab')
    data <- as.data.frame(data)
    data$Ab <- as.numeric(levels(data$Ab))[data$Ab]
    lev <- leveneTest(Ab ~ Group, data = data)
    anv <- oneway.test(Ab~Group, data=data)
    tax$lev_test[i] <- lev$`Pr(>F)`[1]
    tax$pval_size_hexanauplia[i] <- anv$p.value
    #print(p_val$p.value)
    if (mean_small<mean_big_frac){
      tax$PFT[i] <- 'Mesozooplankton'
    } else{
      tax$PFT[i] <- 'Microzooplankton'
    }
    tax$size_hexanauplia[i] <- which.max(c(mean(bigs), mean(meds), mean(smlls)))
  }
  
  if (!is.na(tax$Phytoplankton[tax$Genome_Id==Smg])){
    cond_prok <- grepl('CCKK', colnames(perc_smags))
    cond_small <- grepl('GGZZ', colnames(perc_smags))  | grepl('GGMM', colnames(perc_smags))
    cond_small1 <- grepl('KKQQ', colnames(perc_smags)) | grepl('MMQQ', colnames(perc_smags))
    cond_big1 <- grepl('QQSS', colnames(perc_smags))
    cond_big2 <- grepl('SSUU', colnames(perc_smags))
    mean_prok_frac <- mean(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_prok]))
    mean_small_frac <- mean(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small]))
    mean_small_frac1 <- mean(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small1]))
    mean_small <- mean(c(mean_small_frac,mean_small_frac1 ))
    if (grepl('TOSAG', Smg)){
      Smg <- paste(Smg, '_scaffold', sep='')
    }
    p_val <- wilcox.test(c(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small])), 
                    as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_prok]))
    
    proks <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_prok])
    datab <- cbind(rep('CCKK', length(proks)), proks)
    meds <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small])
    datam <- cbind(rep('GGMM', length(meds)), meds) 
    smlls <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small1])
    datas <- cbind(rep('KKQQ', length(smlls)), smlls)
    bgss <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_big1])
    datab <- cbind(rep('QQSS', length(bgss)), bgss)
    bgss0 <- as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_big2])
    datab0 <- cbind(rep('SSUU', length(bgss0)), bgss0)

    data <- rbind(datab,datam, datas )
    colnames(data) <- c('Group', 'Ab')
    data <- as.data.frame(data)
    data$Ab <- as.numeric(levels(data$Ab))[data$Ab]
    lev <- leveneTest(Ab ~ Group, data = data)
    anv <- oneway.test(Ab~Group, data=data)
    tax$lev_test_phy[i] <- lev$`Pr(>F)`[1]
    tax$pval_size_phy[i] <- anv$p.value

    data <- rbind(datab,datam, datas, datab , datab0)
    colnames(data) <- c('Group', 'Ab')
    data <- as.data.frame(data)
    data$Ab <- as.numeric(levels(data$Ab))[data$Ab]
    lev <- leveneTest(Ab ~ Group, data = data)
    anv <- oneway.test(Ab~Group, data=data)
    tax$lev_test_phy_bis[i] <- lev$`Pr(>F)`[1]
    tax$pval_size_phy_bis[i] <- anv$p.value
   
    #print(p_val$p.value)
    tax$size_phyto_bis[i] <- which.max(c(mean(proks), mean(meds), mean(smlls), mean(bgss), mean(bgss0)))
    tax$size_phyto[i] <- which.max(c(mean(proks), mean(meds), mean(smlls)))
  }
  
  if (is.na(tax$PFT[i])){
    tax$PFT[i] <-'Other'
  }
  i=i+1
}

tax$Hex <- 0
tax$Hex[tax$Best_taxonomy_GENRE=="Marine_Hexanauplia_A"] <-"Marine_Hexanauplia_A"
tax$Hex[tax$Best_taxonomy_GENRE=="Marine_Hexanauplia_B"] <-"Marine_Hexanauplia_B"
tax$pval_size_hexanauplia[!(tax$Best_taxonomy_GENRE %in% c("Marine_Hexanauplia_A", "Marine_Hexanauplia_B"))] <- NA
tax$pval_size_hexanauplia_adj <- p.adjust(tax$pval_size_hexanauplia, 'hochberg')
tax$PFT_hexa <- tax$size_hexanauplia
tax$PFT_hexa[!(tax$Best_taxonomy_GENRE %in% c("Marine_Hexanauplia_A", "Marine_Hexanauplia_B"))] <- NA
tax$PFT_hexa[tax$pval_size_hexanauplia_adj>0.05]<- 'Not classified'
tax$PFT_hexa[tax$PFT_hexa=="2"] <- "Microzooplankton"
tax$PFT_hexa[tax$PFT_hexa=="3"] <- "Microzooplankton"
tax$PFT_hexa[tax$PFT_hexa=="1"] <- "Mesozooplankton"
tax$Hex_bis <-paste(tax$Hex, tax$PFT_hexa)
tax$Hex_bis[tax$Hex_bis=="0 NA"] <- "0"

tax$pval_size_phyto_adj <- p.adjust(tax$pval_size_phy, 'hochberg')
tax$PFT_phy <-tax$size_phyto
tax$PFT_phy[tax$pval_size_phyto_adj>0.05]<- 'Not classified'
tax$PFT_phy[tax$PFT_phy=="2"] <- "Nanophytoplankton big"
tax$PFT_phy[tax$PFT_phy=="3"] <- "Nanophytoplankton small"
tax$PFT_phy[tax$PFT_phy=="1"] <- "Picophytoplankton"


tax$pval_size_phyto_adj_bis <- p.adjust(tax$pval_size_phy_bis, 'hochberg')
tax$PFT_phy_bis <-tax$size_phyto_bis
tax$PFT_phy_bis[tax$pval_size_phyto_adj_bis>0.05]<- 'Not classified'
tax$PFT_phy_bis[tax$PFT_phy_bis=="2"] <- "Nanophytoplankton"
tax$PFT_phy_bis[tax$PFT_phy_bis=="3"] <- "Nanophytoplankton"
tax$PFT_phy_bis[tax$PFT_phy_bis=="1"] <- "Picophytoplankton"
tax$PFT_phy_bis[tax$PFT_phy_bis=="4"] <- "Microphytoplankton"
tax$PFT_phy_bis[tax$PFT_phy_bis=="5"] <- "Microphytoplankton"

tax$PHYLUM_phyto <- tax$Best_taxonomy_PHYLUM
tax$PHYLUM_phyto[is.na(tax$Phytoplankton)]<-NA
tax$PHYLUM_phyto[!is.na(tax$Phytoplankton) & is.na(tax$PHYLUM_phyto)]<- 'Unknown'
tax$PHYLUM_phyto[tax$PHYLUM_phyto=='Bacillariophyta']<-'Diatom'
tax$PHYLUM_phyto[tax$PHYLUM_phyto!='Diatom' & !is.na(tax$PHYLUM_phyto)] <- 'Algae'
tax$Phyto_bis <- paste(tax$PHYLUM_phyto, tax$PFT_phy)
tax$Phyto_bis[tax$Phyto_bis=="NA NA"] <- "0"
tax$Phyto_ter <- paste(tax$PHYLUM_phyto, tax$PFT_phy_bis)
tax$Phyto_ter[tax$Phyto_ter=="NA NA"] <- "0"
saveRDS(tax, 'SMAGs_Statistics_02.rds')


type <-'MAGprok'
taxp <- read.table(paste(type, '_taxo_1.txt',sep=''), sep='\t', header = T)

taxp$size_phyto <- NA
taxp$pval_size_phy <- NA
taxp$lev_test_phy <- NA
taxp$size_phyto_bis <- NA
taxp$pval_size_phy_bis <- NA
taxp$lev_test_phy_bis <- NA


i=1
for (Smg in taxp$Genome_Id){
  
  
  
  if (taxp$Class.gtdbtk[taxp$Genome_Id==Smg]=='Cyanobacteriia'){
    cond_prok <- grepl('CCKK', colnames(perc_magprok))
    cond_small <- grepl('GGMM', colnames(perc_magprok)) | grepl('GGZZ', colnames(perc_magprok))
    cond_small1 <- grepl('KKQQ', colnames(perc_magprok)) | grepl('MMQQ', colnames(perc_magprok))
    cond_big1 <- grepl('QQSS', colnames(perc_magprok))
    cond_big2 <- grepl('SSUU', colnames(perc_magprok))
    mean_prok_frac <- mean(as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_prok]))
    p_val <- wilcox.test(c(as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_small])), 
                    as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_prok]))
    
    proks <- as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_prok])
    datab <- cbind(rep('CCKK', length(proks)), proks)
    meds <- as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_small])
    datam <- cbind(rep('GGMM', length(meds)), meds)
    smlls <- as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_small1])
    datas <- cbind(rep('KKQQ', length(smlls)), smlls)
    bgss <- as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_big1])
    datab <- cbind(rep('QQSS', length(bgss)), bgss)
    bgss0 <- as.numeric(perc_magprok[perc_magprok$AAA_Genome_Id==Smg,cond_big2])
    datab0 <- cbind(rep('SSUU', length(bgss0)), bgss0)    

    data <- rbind(datab,datam, datas )
    colnames(data) <- c('Group', 'Ab')
    data <- as.data.frame(data)
    data$Ab <- as.numeric(levels(data$Ab))[data$Ab]
    lev <- leveneTest(Ab ~ Group, data = data)
    anv <- oneway.test(Ab~Group, data=data)
    taxp$lev_test_phy[i] <- lev$`Pr(>F)`[1]
    taxp$pval_size_phy[i] <- anv$p.value
 
    data <- rbind(datab,datam, datas, datab, datab0 )
    colnames(data) <- c('Group', 'Ab')
    data <- as.data.frame(data)
    data$Ab <- as.numeric(levels(data$Ab))[data$Ab]
    lev <- leveneTest(Ab ~ Group, data = data)
    anv <- oneway.test(Ab~Group, data=data)
    taxp$lev_test_phy_bis[i] <- lev$`Pr(>F)`[1]
    taxp$pval_size_phy_bis[i] <- anv$p.value
    #print(p_val$p.value)
    
    taxp$size_phyto_bis[i] <- which.max(c(mean(proks), mean(meds), mean(smlls), mean(bgss), mean(bgss0)))
    taxp$size_phyto[i] <- which.max(c(mean(proks), mean(meds), mean(smlls)))
  }
  i=i+1
}
taxp$pval_size_phyto_adj <- p.adjust(taxp$pval_size_phy, 'hochberg')
taxp$PFT_phy <-taxp$size_phyto
taxp$PFT_phy[taxp$pval_size_phyto_adj>0.05]<- 'Not classified'
taxp$PFT_phy[taxp$PFT_phy=="2"] <- "Nanophytoplankton big"
taxp$PFT_phy[taxp$PFT_phy=="3"] <- "Nanophytoplankton small"
taxp$PFT_phy[taxp$PFT_phy=="1"] <- "Picophytoplankton"
taxp$Phyto_bis <- paste('Cyanobacteria', taxp$PFT_phy)
taxp$Phyto_bis[taxp$Phyto_bis=="Cyanobacteria NA"] <- "0"

taxp$pval_size_phyto_adj_bis <- p.adjust(taxp$pval_size_phy_bis, 'hochberg')
taxp$PFT_phy_bis <-taxp$size_phyto_bis
taxp$PFT_phy_bis[taxp$pval_size_phyto_adj_bis>0.05]<- 'Not classified'
taxp$PFT_phy_bis[taxp$PFT_phy_bis=="2"] <- "Nanophytoplankton"
taxp$PFT_phy_bis[taxp$PFT_phy_bis=="3"] <- "Nanophytoplankton"
taxp$PFT_phy_bis[taxp$PFT_phy_bis=="1"] <- "Picophytoplankton"
taxp$PFT_phy_bis[taxp$PFT_phy_bis=="4"] <- "Microphytoplankton"
taxp$PFT_phy_bis[taxp$PFT_phy_bis=="5"] <- "Microphytoplankton"
taxp$Phyto_ter <- paste('Cyanobacteria', taxp$PFT_phy_bis)
taxp$Phyto_ter[taxp$Phyto_ter=="Cyanobacteria NA"] <- "0"
saveRDS(taxp,'MAGprok_Statistics.rds')

all_MAGS <- data.frame(Genome_Id=c(as.character(taxp$Genome_Id), as.character(tax$Genome_Id)), Phyto_bis=c(taxp$Phyto_bis, tax$Phyto_bis),
                       Phyto_ter=c(taxp$Phyto_ter, tax$Phyto_ter),
                       KINGDOM=c(as.character(taxp$Domain.anvio), tax$Best_taxonomy_KINGDON),
                       PHYLUM=c(as.character(taxp$Phylum.gtdbtk), tax$Best_taxonomy_PHYLUM), 
                       CLASS=c(as.character(taxp$Class.gtdbtk), tax$Best_taxonomy_CLASS),
                       ORDER=c(as.character(taxp$Order.gtdbtk), tax$Best_taxonomy_ORDER), 
                       FAMILY=c(as.character(taxp$Family.gtdbtk), tax$Best_taxonomy_FAMILY),
                       GENRE=c(as.character(taxp$Genus.gtdbtk), tax$Best_taxonomy_GENRE))

colnames(perc_magprok)[1] <-'SMAG'
perc_smags <- perc_smags[,match(colnames(perc_magprok), colnames(perc_smags))]
perc_smags_all <- rbind(perc_magprok, perc_smags)
saveRDS(all_MAGS, 'all_MAGs_statistics.rds')
saveRDS(perc_smags_all, 'all_MAGs_percent_mapped_read.rds')
write.table(perc_smags_all, 'all_MAGs_percent_mapped_read.txt')

data_SMAGs <- read.table(paste('Occurence_SMAGs.txt',sep=''), header = T)
data_MAGprok <- read.table(paste('Occurence_MAGprok.txt',sep=''), header = T)
colnames(data_MAGprok)[1] <-'SMAG'
data_SMAGs <- data_SMAGs[,match(colnames(data_MAGprok), colnames(data_SMAGs))]
data_smags_all <- rbind(data_MAGprok, data_SMAGs[2:dim(data_SMAGs)[1],])
data_smags_all <- rbind(data_SMAGs[1,], data_smags_all)
write.table(data_smags_all, 'Occurence_all_MAGs.txt')

