library("readxl")
library('maptools')
library('plotrix')

setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5")
cambria <- readRDS('cambria.rds')
source('axis_map0.R')
source('hide_arctic.R')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
cross_val <- read_xlsx('models_cross_validation.xlsx')
cross_val <- cross_val[1:38,]
cross_val$Fraction[cross_val$Fraction=='43952']='5-20'

cross_val$auc_bt[is.na(cross_val$auc_bt)]<-0
labs1 <- paste(cross_val$Fraction, cross_val$Gen)
labs2 <- paste(cross_val$Fraction, cross_val$Gen, sep='_')
v <- matrix(c(cross_val$auc_gam, cross_val$auc_bt, 
            cross_val$auc_nn, cross_val$auc_rf), ncol=38, byrow = T,
            dimnames = list(c("gam", "gbm", "nn", 'rf'), labs1))
fraction <-unique(cross_val$Fraction)
letters<-rev(c('A', 'B', 'C', 'D', 'E', 'F'))
valid_model <- function(u){
  return(sum(u>0.65))
}
n_valid <- apply(v, 2, valid_model)
valids <- as.numeric(n_valid>2)+1
saveRDS(labs2[valids==2],'valids.rds')
tes <- c('', '*')
h <- barplot(v, beside = T)
dev.off()
xlim <- max(h)
order <- c(1:38)
v <- v[, order]
valids <- valids[order]
labs <- NULL
for (u in colnames(v)){
  frc <- strsplit(u, ' ')[[1]][1]
  let <- letters[match(frc, fraction)]
  labs <- append(labs, paste(let, strsplit(u, ' ')[[1]][2], sep=''))
}
colnames(v)<- labs
pdf(family="Helvetica",'barplot_auc_cross_val.pdf', width = 15, height = 10)
par(mar=c(10,10,10,10))
xu <- barplot(v, beside = T, col=c('red','green', 'blue', 'violet' ), 
        names.arg = toupper(colnames(v)), legend.text = TRUE,
        args.legend=list(x=xlim+15,y=1,bty = "n"), cex.names = 1.5,
        las=2, ylim=c(0,1.1), cex.axis = 1.5)
abline(h=0.65, col='black')
text(x = (xu[3,]+xu[2,])/2, y=1.03,labels =  tes[valids], cex = 2)
dev.off()
fractions0 <- list(1:3, 4:6, 7:9, 10:16, 17:22, 23:27 )
saveRDS(fractions0, 'fractions0.rds')

df_init <- read_excel("Genocenoses_env_parameters_all_woa.xlsx")
to_remove <- readRDS('excluded_niches.rds')
df <- readRDS('Genocenoses_env_parameters_woa_scaled.rds')
df$Fraction[df$Fraction=='43952']<-'5-20'
df$labs <- paste(df$Fraction, df$Genocenose)
df$labs0 <- labs[match(df$labs, labs1)]

df_init$id <- paste(df_init$Fraction, df_init$Genocenose, sep='_')
df_init$Fraction[df_init$Fraction=='43952']<-'5-20'
df_init$labs <- paste(df_init$Fraction, df_init$Genocenose)
df_init$labs0 <- labs[match(df_init$labs, labs1)]
df_init <- df_init[df_init$Fraction!='all',]


counts <- NULL
for (l in labs[valids==1]){
  count <- sum(df$labs0==l)
  counts <-append(counts, count)
}

counts1 <- NULL
for (l in labs[valids==2]){
  count <- sum(df$labs0==l)
  counts1 <-append(counts1, count)
}

mat_to_write <- matrix(c(mean(counts), median(counts),sd(counts), 
                         mean(counts1),median(counts1), sd(counts1)), ncol=2)

pdf(family="Helvetica",'niches_validated_vs_non_validated.pdf')
hist(counts1, col=scales::alpha('green', 0.5), breaks=20, ylim=c(0, 9))
hist(counts, col=scales::alpha('blue', 0.5), add=T)
dev.off()

write.table(mat_to_write, file = 'niches_validated_vs_non_validated.txt', row.names = F, col.names = F)

df0 <- df_init[,c(1,2,4,5,20)]
stations <- df0[!duplicated(df0$Station),c(1,3,4)]
#df0 <- df0[df0$labs0 %in% labs[as.logical(valids-1)],]
fractions <-c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
valid_colors <- NULL
for (st in unique(stations$Station)){
  valid <- NULL
  df1 <- df0[df0$Station==st,]
  for (let in letters){
    if (sum(grepl(let, df1$labs0))==1){
      
      if (df1$labs0[grep(let, df1$labs0)] %in% labs[as.logical(valids-1)]){
        valid <- append(valid, 'green')
      } else{
        valid <- append(valid, 'red')
      }
    } else{
      if (sum(grepl(fractions[letters==let], df1$Fraction))==1){
        valid <- append(valid, 'black')
      } else{
        valid <- append(valid, 'white')
      }
    }
  }
  valid_colors <- append(valid_colors, valid)
}
valid_colors <- matrix(valid_colors, ncol=6, byrow = T)
stations <- cbind(stations, valid_colors)

pdf(family="Helvetica",'valid_vs_non_valid_niches.pdf',width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
for (st in unique(stations$Station)){
  co <- stations$Station==st
  floating.pie(stations$Long[co], stations$Lat[co], x=rep(1/6, 6), col=valid_colors[co,], radius=3)
}
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",'valid_vs_non_valid_niches_SUR.pdf',width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
for (st in unique(stations$Station)){
  if (grepl('SUR', st)){
    co <- stations$Station==st
    floating.pie(stations$Long[co], stations$Lat[co], x=rep(1/6, 6), col=valid_colors[co,], radius=3.5)
  }
}
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",'valid_vs_non_valid_niches_DCM.pdf',width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
for (st in unique(stations$Station)){
  if (grepl('DCM', st)){
    co <- stations$Station==st
    floating.pie(stations$Long[co], stations$Lat[co], x=rep(1/6, 6), col=valid_colors[co,], radius=3.5)
  }
}
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",'legend_valid_vs_non_valid_niches.pdf')
plot(0,0,axes=F, col=scales::alpha('white',0), xlab='', ylab='', xlim=c(0,5), ylim=c(0,5))
floating.pie(2.5,2.5, x=rep(1/6, 6), col='white')
dev.off()

pdf(family="Helvetica",'legend_valid_vs_non_valid_niches1.pdf')
plot(0,0,axes=F, col=scales::alpha('white',0), xlab='', ylab='', xlim=c(0,5), ylim=c(0,5))
floating.pie(2.5,2.5, x=rep(1/6, 6), col=c('green', 'white', 'black', 'white', 'red', 'white'))
dev.off()
