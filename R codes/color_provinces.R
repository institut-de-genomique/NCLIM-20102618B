library('RColorBrewer')

cambria <- readRDS('cambria.rds')
CA <- read.table('Projections/covered_areas_all.txt')
colnames(CA) <- c("Metacommunity" , "Fraction", 
                  "Model","Area woa","Area 2006 (km2)"  ,
                  "Area 2090 (km2)",   "Delta area (2090-2006)")
letters <- c('A', 'B', 'C', 'F', 'E', 'D')
fracs0 <- unique(CA$Fraction)
CA$Genomic_province <- paste(letters[match(CA$Fraction, fracs0)], CA$Metacommunity,
                             sep='')
provinces <- CA$Genomic_province[order(unique(CA$Genomic_province))]


qual_col_pals1 = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector1 = as.vector(unlist(mapply(brewer.pal, qual_col_pals1$maxcolors, rownames(qual_col_pals1))))
# set.seed(3)
# col_vec <- col_vector1[sample(x = 1:99, size = 27)]
pdf(family="Helvetica",'colors.pdf')
plot(0,0, col='white', xlim=c(0,30), ylim=c(0,30))
legend('topright', legend=1:length(col_vector1),fill = col_vector1, ncol = 5)
dev.off()

colors_provinces <- rep(list(NULL), 6)
names(colors_provinces) <- c('180-2000', '20-180', '43952', '0.8-5', '0.22-3', '0-0.2')

#colors_provinces$`0.8-5` <- col_vector1[c(1:3, 60, 5:7)]
colors_provinces$`0.8-5` <- col_vector1[c(1, NA, 2, 3,NA, NA, NA, 60, 5:7)]
#colors_provinces$`0.22-3` <- col_vector1[c(9,11,14,17, 22,23)]
colors_provinces$`0.22-3` <- col_vector1[c(NA, 9,11,NA,14,17, 22,23)]
#colors_provinces$`180-2000` <- col_vector1[c(49,50,53)]
colors_provinces$`180-2000` <- col_vector1[c(49,NA,NA,NA,50,NA,NA,20)]
#colors_provinces$`20-180` <- col_vector1[c(29,59,63)]
colors_provinces$`20-180` <- col_vector1[c(29,NA,NA,NA,59,63)]
# colors_provinces$`43952` <- col_vector1[c(61, 57,58)]
colors_provinces$`43952` <- col_vector1[c(NA,NA,61, 57,NA,58)]
# colors_provinces$`0-0.2` <- col_vector1[c(68,69, 72:74)]
colors_provinces$`0-0.2` <- col_vector1[c(NA,NA,68,69,NA, 72:74)]
saveRDS(colors_provinces, 'colors_provinces.rds')

dat <- NULL
c=1
letters <- c('F', 'E', 'D', 'C', 'B', 'A')
for (fr in names(colors_provinces)){
  prvs <- paste(letters[c], c(1:length(colors_provinces[[fr]])), sep='')
  cols <- colors_provinces[[fr]]
  dt <- cbind(prvs, cols)
  dat <- rbind(dat, dt)
  c=c+1
}
dat <- dat[!is.na(dat[,2]),]
colnames(dat)<- c('Province', 'Color')
write.table(dat, 'color_provinces.txt', row.names=F)



