library('circlize')

transitions <- readRDS('transitions_bis.rds')
trans <- readRDS('considered_transitions.rds')
trans <- as.character(sapply(trans, FUN = function(x){strsplit(x = x, split = ' ')[[1]][1]}))
data_trans <- rep(list(NULL), 6)
weights_cos <- readRDS('model-mean_weights_cos.rds')
frcs <- c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
names(data_trans)<- frcs
for (u in 1:6){
  trs <- transitions[[u]]
  for (tr in unique(trs)[unique(trs)!="FALSE"]){
    p1 <- strsplit(tr, '->')[[1]][1]
    p2 <- strsplit(tr, '->')[[1]][2]
    score <- sum(111*111*weights_cos*c(trs==tr)*1)
    data_trans[[frcs[[u]]]] <- rbind(data_trans[[u]], c(p1,p2, score))
  }
}
col_prov <- read.table('color_provinces.txt', header = T)
comps <- readRDS('SMAGs_provinces_compositional_functions_genomics_trophic.rds')
comps_dia <- readRDS('MAGprok_provinces_compositional_functions_diazotrophy.rds')
col_list <- list(c('darkorange', 'darkorange4', 'darkgoldenrod1', 'yellow'),c('darkviolet', 'mediumpurple1', 'magenta') ,
                 c('darkgreen', 'mediumseagreen', 'green', 'greenyellow'),
                 c('dodgerblue2', 'mediumpurple4', 'darkblue', 'darkturquoise', 'cyan', 'lightslateblue', 'slateblue4', 'skyblue', 
                   'blueviolet', 'cadetblue'),
                 c('red', 'darkred', 'coral', 'deeppink', 'indianred1', 'orangered', 'rosybrown', 'violetred', 'bisque'),
                 c('saddlebrown', 'peru', 'sandybrown', 'rosybrown', 'moccasin'))
names(col_list)<- frcs
letters <- list('F', 'E', 'D', 'C', 'B', 'A')
names(letters)<- frcs


tax <- read.table('MAGprok_taxo_1.txt', sep='\t', header = T)
tax$N_Fixation_Phylum <- as.character(levels(tax$N_Fixation_Phylum))[tax$N_Fixation_Phylum]
Diaz_list <- unique(tax$N_Fixation_Phylum)
Diaz_list <- Diaz_list[Diaz_list!=""]

tax_s <- readRDS('SMAGs_Statistics_02.rds')
Hex_list_bis <- unique(tax_s$Hex_bis)
Hex_list_bis <- Hex_list_bis[!is.na(Hex_list_bis) & Hex_list_bis!='0']
Hex_list_bis <- Hex_list_bis[c(4, 3, 2, 1, 5)]

sig_diaz <- read.table('significant_differences_diazotrophy_wilcox.txt', header = T)
sig_hex <-  read.table('significant_differences_hexanauplia_wilcox.txt', header = T)
sig_hex <- sig_hex[!is.na(sig_hex$Trans),]
circ_trans <- function(type, frac, name0, compos, fac, climates){
  data <- compos[[type]][[frac]]
  data0 <- matrix(unlist(data), byrow = T, nrow=length(data))
  sectors = names(data)
  if (frac!='0.8-5'){
    sectors=sectors[order(sectors)]
  } else{
    sectors <- sectors[c(7,5,2, 4, 1, 6, 3)]
  }
  sel <- apply(data0, 2, sum)
  sel <- which(sel!=0)
  data <-data[sectors]
  
  x_lists0=rep(list(c(0,4)), length(data))
  x_lists=data
  trsi <- data_trans[[frac]]
  # s1 = factor(sectors)
  if (name0=='hexanauplia'){
    colos=c('darkred','red','lightpink' ,'lightblue','blue')
  } else{
    colos=rainbow(dim(data0)[2])[sel]
  }
  colos_prov=col_prov$Color[match(sectors, col_prov$Province)]
  mx =signif(max(log(unlist(data)*fac+1)), 1)
  mxo =max(unlist(data))
#   seco <- seq(0, round(mx,1), length.out = 4)
#   labs <- (exp(seco)-1)/100
  if (name0=='diaz'){
    sig <- sig_diaz
    grouping=Diaz_list[sel]
    if (frac %in% c('20-180', '180-2000')){
      labs= c(0.001, 0.01, 0.1)
      seco= log(c(0.001, 0.01, 0.1)*fac+1)
#       if (signif(mxo,1)!=0.1 ){
#         labs <- c(labs, mxo)
#         seco <- c(seco, mx)
#       }
      
    } else{
      labs= c(0.00001,0.0001,0.001, 0.01, 0.1)
      seco= log(c(0.00001,0.0001,0.001, 0.01, 0.1)*fac+1)
#       if (signif(mxo,1)!=0.1 & mxo>7*labs[which.min(abs(mxo-labs))]){
#         labs <- c(labs, mxo)
#         seco <- c(seco, mx)
#       }
    }
    
  } else{
    sig <- sig_hex
    grouping=Hex_list_bis
    labs= c( 0.01, 0.1, 1, 10)
    seco= log(c( 0.01, 0.1, 1, 10)*fac+1)
#     if (mxo<7 & mxo !=10){
#       labs <- c(labs, mxo)
#       seco <- c(seco, mx)
#     }
  }
  labs <- labs[order(seco)]
  seco <- seco[order(seco)]
  trsi <- data_trans[[frac]]
  trsi <- as.data.frame(trsi)
  colnames(trsi)<- c('from', 'to', 'value')
  trsi$value <- as.numeric(levels(trsi$value))[trsi$value]
  let <- letters[[frac]]
  trs_possibles <-trans[grepl(let, trans)]
  trsi$trans <- paste(trsi$from, trsi$to, sep='->')
  pdf(paste('circos_transitions_',name0,'_',frac, '.pdf', sep=''))
  circos.initialize(sectors, xlim = c(0, 8))
  circos.track(ylim = c(0,2 ),track.height=0.13, panel.fun = function(x, y) {
    i =CELL_META$sector.numeric.index
    ylim = CELL_META$ylim
    circos.text(2, (ylim[1]+ylim[2])/2, labels = climates[i], cex = 1.3)
  }, bg.border=NA)
  circos.track(ylim = c(0,mx+1*0.3*mx ),track.height=0.32, panel.fun = function(x, y) {
    i =CELL_META$sector.numeric.index
    value = x_lists[[i]][sel]
    #print(value)
    tp = log(value*fac+1)
    xs=(2:(length(value)+1))/2
    circos.barplot(value = tp, pos =  xs, col = colos, bar_width = 0.5)
    circos.yaxis(side='left',at = seco, labels = signif(labs, 1), labels.cex = 1)
    if (frac %in% c('0.22-3', '0.8-5')){
      g = -3
    } else{
      g= -1.2
    }
    circos.text(x=g, y=(mx+1*0.3*mx)/2,facing='inside', labels='%')
    
    sector <- sectors[i]
    if (sector %in% trsi$to){
      c=rep(1, length(grouping))
      j=rep(1, length(grouping))
      xxs<-NULL
      yys <- NULL
      cols <- NULL
      for (tr in trsi$trans[trsi$to==sector]){
        if (tr %in% sig$Trans){
          print(tr)
          print(sector)
          gs <- sig$Clade[sig$Trans==tr]
          for (a in gs){
            ind <- which(grouping==a)
#             print(a)
#             print(ind)
#             print(grouping)
            tot <- c(1,-1)
            if (frac=='0.8-5' & name0=='diaz'){
              ec=0.2
            } else{
              ec=0.1
            }
            xxs <- append(xxs, xs[ind]+tot[c[ind]%%2+1]*ec)
            yys <- append(yys, tp[ind]+0.3*j[ind]*tp[ind])
            if (tr %in% trs_possibles){
              cl = col_list[[frac]][match(tr,trs_possibles)]
            } else{
              cl <- 'gray'
            }
            # print(cl)
            cols <- append(cols, cl)
            if (c[ind]%%2==0){
              j[ind]=j[ind]+1
            }
            c[ind]=c[ind]+1
          }
          
        }
      }
#       print(xxs)
#       print(yys)
      circos.points(x=xxs,y=yys ,col=rep('black',length(cols)),bg=cols, cex=1, pch=25)
      #circos.points(x=xxs,y=yys ,col=cols, cex=1.5, pch='*')
      # circos.points(x=xxs,y=yys ,col=rep('black',length(cols)),bg=cols, cex=1, pch=23)
    }
    # print('')
    #print(warnings())
#     if (length(xxs)>2){
#       print(xxs)
#       xxs[seq(1,length(xxs), 2)] <- xxs[seq(1,length(xxs), 2)]-0.25
#       xxs[seq(2,length(xxs), 2)] <- xxs[seq(2,length(xxs), 2)]+0.25
#     }
    
    mxa <- 4
    for (l in seco){
      if (mx+2*0.2*mx >= l){
        for (i in seq(1,60,2)){
          circos.lines(x = seq(i*mxa/60, (i+1)*mxa/60, length.out = 10), y = rep(l, 10), lwd = 0.5 , type = 'l')
        }
      }
    }
  }, bg.border=NA)
  circos.track(ylim = c(0, 2), track.height=0.13, panel.fun = function(x,y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    i =CELL_META$sector.numeric.index
    breaks = x_lists0[[i]]
    n_breaks = length(breaks)
    circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
                breaks[-1], rep(ylim[2], n_breaks - 1),
                col = as.character(colos_prov[i]), border = NA)
    circos.text(2, (ylim[1]+ylim[2])/2, labels = sectors[i], cex = 1.3)
  }, bg.border=NA)
#   cols=NULL
#   for (i in 1:dim(trsi)[1]){
#     p1 <- trsi[i,1]
#     p2 <- trsi[i,2]
#     t <- paste(p1, p2, sep='->')
#     if (t %in% trs_possibles){
#       cols <- append(cols, scales::alpha(col_list[[frac]][match(t,trs_possibles)], 0.5))
#     } else{
#       cols <- append(cols, 'gray')
#     }
#   }
  
  #chordDiagram(trsi, col=cols)
  #trsi <- trsi[order(trsi$value, decreasing = F),] 
  co_o <- rep(0, length(sectors))
  co_a <- rep(0, length(sectors))
  for (i in 1:dim(trsi)[1]){
    p1 <- trsi[i,1]
    p2 <- trsi[i,2]
    t <- paste(p1, p2, sep='->')
    widths_arr = (as.numeric(trsi[i,3])*1.5-0.3)/max(as.numeric(trsi[,3]))+0.3
    l_wds=(as.numeric(trsi[i,3])*20-1.5)/max(as.numeric(trsi[,3]))+1.5
    if (t %in% trs_possibles){
      
      circos.link(p1, 0+co_o[which(sectors==p1)]*0.5, p2, 4-co_a[which(sectors==p2)]*0.5,
                  col = col_list[[frac]][match(t,trs_possibles)], arr.width = widths_arr,
                  directional=1, lwd   = l_wds)
                  #lwd = as.numeric(trsi[i,3])*2/max(as.numeric(trsi[,3])), )lwd=as.numeric(trsi[i,3])*10/max(as.numeric(trsi[,3]))
      co_o[which(sectors==p1)]= co_o[which(sectors==p1)]+1
      co_a[which(sectors==p2)]=co_a[which(sectors==p2)]+1
    } else{
      ld <-  as.numeric(trsi[i,3])*10/max(as.numeric(trsi[,3]))
      # if (ld>0.5){
        circos.link(p1, 0+co_o[which(sectors==p1)]*0.33, p2, 4-co_a[which(sectors==p2)]*0.33,
                    col = 'gray', arr.width = widths_arr,
                    directional=1, lwd   =l_wds)
        #lwd = as.numeric(trsi[i,3])*2/max(as.numeric(trsi[,3])), )lwd=as.numeric(trsi[i,3])*10/max(as.numeric(trsi[,3]))
        co_o[which(sectors==p1)]= co_o[which(sectors==p1)]+1
        co_a[which(sectors==p2)]=co_a[which(sectors==p2)]+1
      #}
    }
  }
  dev.off()
}

types <- c(rep('Abundance Hex_bis', 3), rep('Abundance Diazotrophy', 5))
name0s <- c(rep('hexanauplia', 3), rep('diaz', 5))
fracs <- c('180-2000', '20-180', '5-20', '180-2000', '20-180', '5-20', '0.8-5', '0.22-3')
list_climate <- annotations0 <- list('180-2000'=c('Polar', 'Temperate','Tropico-equatorial'), 
                                     '20-180'=c('Polar', 'Tropico-equatorial', 'Temperate'),
                                     '5-20'=c('Temperate', 'Equatorial',  'Tropical'),
                                     '0.8-5' =c('Polar', 'Subtropical', 'Subtropical', 'Temperate', 'Equatorial', 'Subtropical', 'Tropical'),
                                     '0.22-3'=c('Subtropical', 'Equatorial', 'Tropical','Temperate', 'Temperate', 'Polar'),
                                     '0-0.2'=c('Tropical', 'Subtropical', 'Equatorial', 'Temperate', 'Temperate'))

for (i in 1:length(fracs)){
  type=types[i]
  frac=fracs[i]
  name0=name0s[i]
  clims=list_climate[[frac]]
  if (name0=='hexanauplia'){
    circ_trans(type,frac,name0, comps, 100, clims)
  } else if (name0=='diaz'){
    if (frac %in% c('20-180', '180-2000')){
      circ_trans(type,frac,name0, comps_dia, 100, clims)
    } else{
      circ_trans(type,frac,name0, comps_dia, 100000, clims)
    }
  }
}


