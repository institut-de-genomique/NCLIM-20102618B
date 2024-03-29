library('gplots')
library('RColorBrewer')
library('reshape2')
library('data.table')

source('vioplot.R')
cambria <- readRDS('cambria.rds')
RL_w <- read.table('data_relinf_dalex.txt')
letters <- rev(c('A', 'B', 'C', 'D', 'E', 'F'))
Fractions = c('180-2000', '20-180', '43952', '0.8-5', '0.22-3', '0-0.2')
colnames(RL_w)[c(1,2)]<- c('Frac', 'Metacommunity')
RL_w$clust <- paste(letters[match(RL_w$Frac,Fractions)], RL_w$Metacommunity, sep='')
RL_w_mat <- reshape2::acast(RL_w, RL_w$clust~RL_w$vars, value.var = 'Inf.')

br <- c(seq(0, 75, length.out = 101))
set = rev(colorRampPalette(c('darkred', 'red','orange', 'green','blue' , 'lightblue'))(100))
order <- c(25:27, 22:24, 19:21, 12, 15:18, 13,14, 6:11, 1:5)
order1 <- c(7,4,5,2, 3,1,6)
pdf(family="Helvetica",'heatmap_relinf_niches.pdf', width = 10, height = 10)
heatmap.2(RL_w_mat[order,order1], trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1, col=set ,margins = c(13.5,13.5),
          breaks=br, symkey = F, cexRow=1)
dev.off()

names <- rownames(RL_w_mat[order,])
relinf_cc <-read.table('drivers_solo_niches.txt')
relinf_cc <- relinf_cc[,-1]
rownames(relinf_cc)<-names

relinf_cc_table <- melt(setDT(relinf_cc, keep.rownames = TRUE), 'rn')
relinf_cc_table <- as.matrix(relinf_cc_table)
colnames(relinf_cc_table)<-c('clust', 'Param', 'Inf.')
relinf_cc_table<- as.data.frame(relinf_cc_table)
relinf_cc_table$Inf. <- as.numeric(levels(relinf_cc_table$Inf.))[relinf_cc_table$Inf.]
pairwise.wilcox.test(relinf_cc_table$Inf., g = relinf_cc_table$Param)

for (l in letters){
  ta <- relinf_cc_table[grep(l, relinf_cc_table$clust),]
  print(l)
  print(pairwise.wilcox.test(relinf_cc_table$Inf., g = relinf_cc_table$Param))
}


relinf_cc<- as.matrix(relinf_cc[,2:8])
rownames(relinf_cc)<-names
pdf(family="Helvetica",'heatmap_relinf_niches_cc.pdf', width = 10, height = 10)
heatmap.2(as.matrix(relinf_cc*100), trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1, col=set ,margins = c(13.5,13.5),
          breaks=br, symkey = F, cexRow=1)
dev.off()

data_rl <- as.data.frame(relinf_cc)
x1 <- 100*data_rl$T
x2 <- 100*data_rl$Sal
x3 <- 100*data_rl$Si
x4 <- 100*data_rl$NO3
x5 <- 100*data_rl$Phos
x6 <- 100*data_rl$Fe
# x7 <- data_rl$Inf.[data_rl$vars=='Chla']
x8 <- 100*data_rl$SI_NO3
N <- 7 
colores = brewer.pal(N, "Dark2")
pdf(family="Helvetica",'relative_influence_violin_cc.pdf')
vioplot(x1, x2, x3,x4,x5,x6,x8, names=c('T', 'Sal', 'Si', 'NO3', 'PO4', 'Fe','SI_NO3'), col =colores, ylim=c(0,100))
dev.off()

jackson <- NULL
for (i in 1:7){
  me <- mean(relinf_cc[,i])*100
  med <-  median(relinf_cc[,i])*100
  jackson <- append(jackson, c(colnames(relinf_cc)[i],me, med))
}

jackson <- matrix(jackson, nrow=7, byrow = T)
write.table(jackson, 'relative_influences_drivers_cc.txt')
