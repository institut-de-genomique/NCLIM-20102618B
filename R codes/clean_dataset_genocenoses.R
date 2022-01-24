library("readxl")
library('randomForest')

df <- read_excel("Genocenoses_env_parameters_all_woa.xlsx")
fractions <- c('180-2000', '20-180', '43952','0.8-5','0.22-3', '0-0.2')
df$id <- paste(df$Fraction, df$Genocenose, sep='_')
lagoons <- c('180-2000_4', '180-2000_6', '20-180_4')
df <- df[df$Fraction!='all',]
to_remove <- NULL
for (i in unique(df$id)){
  d <- df[df$id==i,]
  if (dim(d)[1]<4 | grepl('NA', i) | i %in% lagoons){
    to_remove <- append(to_remove, i)
  }
}
saveRDS(to_remove, 'excluded_niches.rds')
df <- df[!(df$id %in% to_remove),]
df0 <- df
df1 <- df[!duplicated(df$Station),]
df <- as.data.frame(df)
df0 <- as.data.frame(df0)
df1 <- as.data.frame(df1)
means_sds <- rep(list(NA), 9)
for (u in 6:14){
  df1[,u] <- randomForest::na.roughfix(df1[,u])
  me <- mean(df1[,u], na.rm = T)
  sde <- sd(df1[,u], na.rm = T)
  means_sds[[u]] <- c(me, sde)
  df0[,u] <- (df0[,u]-me)/sde
}
saveRDS(df, 'Genocenoses_env_parameters_woa.rds')
saveRDS(df0, 'Genocenoses_env_parameters_woa_scaled.rds')
saveRDS(means_sds, 'means_sds_woa_tara.rds')
