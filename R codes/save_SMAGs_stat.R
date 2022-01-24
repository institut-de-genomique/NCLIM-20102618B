library('readxl')
library('car')
type <-'SMAGs'
tax <- read_xlsx(paste(type, '_Statistics_02.xlsx', sep=''))
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
    p_val <- t.test(c(as.numeric(perc_smags[perc_smags$SMAG==Smg,cond_small])), 
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
saveRDS(tax, 'SMAGs_Statistics_02.rds')


type <-'MAGprok'
tax <- read.table(paste(type, '_taxo_1.txt',sep=''), sep='\t', header = T)
