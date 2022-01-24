rgb_map <- function(d, dimension, outputfile){
  cambria <- readRDS('cambria.rds')
  pca <- cmdscale(d, k=3,eig = TRUE, add = TRUE)
  x <- pca$points[,1]
  y <- pca$points[,2]
  z <- pca$points[,3]
  
  e1 <- pca$eig[1]
  e2 <-pca$eig[2]
  e3 <- pca$eig[3]
  
  e1_var <- e1/sum(pca$eig)
  e2_var <- e2/sum(pca$eig)
  e3_var <- e3/sum(pca$eig)
  print(sum(e1_var, e2_var, e3_var))
  print(outputfile)
  # 
  # x_n <- (x - min(x))/max(x - min(x))*255
  # y_n <- (y - min(y))/max(y - min(y))*255
  # z_n <- (z - min(z))/max(z - min(z))*255
  x_n <- (x/max(abs(x))+1)*255/2
  if (dimension==3 | dimension ==2){
    y_n <- ((e2/e1)*y/max(abs(y))+1)*255/2
  } else{
    y_n <- rep(0, length(x))
  }
  if (dimension==3){
    z_n <- ((e3/e1)*z/max(abs(z))+1)*255/2
  } else{
    z_n <- rep(0, length(x))
  }
  # 
  # plot(x_n, e2*y_n/e1, col=rgb(x_n, e2*y_n/e1, e3*z_n/e1, maxColorValue = 255))
  # 
  # colors <-rgb(x_n, e2*y_n/e1, e3*z_n/e1, maxColorValue = 255)
  # if (is.nan(sum(z_n))){
  #   z_n <- rep(0, length(x))
  # }
  # if (is.nan(sum(y_n))){
  #   y_n <- rep(0, length(x))
  # }
  colors <-rgb(x_n, y_n, z_n, maxColorValue = 255)
  source('hide_arctic.R')
  pdf(family="Helvetica",file=paste(outputfile),width=10,height=4.065)
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  lats <- NULL
  longs <- NULL
  locs <- NULL
  done_stat <- NULL
  col_dcm <- NULL
  col_srf <- NULL
  v_loc <- NULL
  for(j in 1:dim(d)[1]){
#     if (proj=='tara'){
#       c = c(df2$Lat.[j], df2$Long.[j])
#       points(x=c[2], y=c[1], col=NA,bg= colors[j] ,pch=21,cex=20)
#       text(x=c[2], y=c[1], labels = df2$Station[j]  ,cex=7)
#     } else if (proj =='bioadv'){
#       st=str_sub(start=1, end=nchar(colnames(d)[j])-3, colnames(d)[j])
#       c = c(df2$Lat.[df2$Station==st], df2$Long.[df2$Station==st])
#       points(x=c[2], y=c[1], col=NA,bg= scales::alpha(colors[j], 1) ,pch=21,cex=20)
#       # text(x=c[2], y=c[1], labels = st  ,cex=7)
#     }
    st <- colnames(d)[j]
    h <- which(Lomb_sts$Station==st)
    c= c(Lomb_sts$Latitude[h], Lomb_sts$Longitude[h])
    points(x=c[2], y=c[1], col=NA,bg= colors[j] ,pch=21,cex=20)
  }
  
  dev.off()
  #rgb_map <- list(e1_var, e2_var, e3_var)
}
upper.half.circle <- function(x,y,r,nsteps=100,...){
  rs <- seq(0,pi,len=nsteps)
  xc <- x+r*cos(rs)
  yc <- y+r*sin(rs)
  polygon(xc,yc, border = NA,...)
}

lower.half.circle <- function(x,y,r,nsteps=100,...){
  rs <- seq(0,pi,len=nsteps)
  xc <- x-r*cos(rs)
  yc <- y-r*sin(rs)
  polygon(xc,yc,border = NA,...)
}
rgb_map0 <- function(d, dimension, outputfile){
  cambria <- readRDS('cambria.rds')
  pca <- cmdscale(d, k=3,eig = TRUE, add = TRUE)
  x <- pca$points[,1]
  y <- pca$points[,2]
  z <- pca$points[,3]
  
  e1 <- pca$eig[1]
  e2 <-pca$eig[2]
  e3 <- pca$eig[3]
  
  e1_var <- e1/sum(pca$eig)
  e2_var <- e2/sum(pca$eig)
  e3_var <- e3/sum(pca$eig)
  print(sum(e1_var, e2_var, e3_var))
  print(outputfile)
  # 
  # x_n <- (x - min(x))/max(x - min(x))*255
  # y_n <- (y - min(y))/max(y - min(y))*255
  # z_n <- (z - min(z))/max(z - min(z))*255
  x_n <- (x/max(abs(x))+1)*255/2
  if (dimension==3 | dimension ==2){
    y_n <- ((e2/e1)*y/max(abs(y))+1)*255/2
  } else{
    y_n <- rep(0, length(x))
  }
  if (dimension==3){
    z_n <- ((e3/e1)*z/max(abs(z))+1)*255/2
  } else{
    z_n <- rep(0, length(x))
  }
  # 
  # plot(x_n, e2*y_n/e1, col=rgb(x_n, e2*y_n/e1, e3*z_n/e1, maxColorValue = 255))
  # 
  # colors <-rgb(x_n, e2*y_n/e1, e3*z_n/e1, maxColorValue = 255)
  # if (is.nan(sum(z_n))){
  #   z_n <- rep(0, length(x))
  # }
  # if (is.nan(sum(y_n))){
  #   y_n <- rep(0, length(x))
  # }
  colors <-rgb(x_n, y_n, z_n, maxColorValue = 255)
  
  pdf(family="Helvetica",file=paste(outputfile),width=10,height=4.065)
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  lats <- NULL
  longs <- NULL
  locs <- NULL
  done_stat <- NULL
  col_dcm <- NULL
  col_srf <- NULL
  v_loc <- NULL
  for(j in 1:dim(d)[1]){
    #     if (proj=='tara'){
    #       c = c(df2$Lat.[j], df2$Long.[j])
    #       points(x=c[2], y=c[1], col=NA,bg= colors[j] ,pch=21,cex=20)
    #       text(x=c[2], y=c[1], labels = df2$Station[j]  ,cex=7)
    #     } else if (proj =='bioadv'){
    #       st=str_sub(start=1, end=nchar(colnames(d)[j])-3, colnames(d)[j])
    #       c = c(df2$Lat.[df2$Station==st], df2$Long.[df2$Station==st])
    #       points(x=c[2], y=c[1], col=NA,bg= scales::alpha(colors[j], 1) ,pch=21,cex=20)
    #       # text(x=c[2], y=c[1], labels = st  ,cex=7)
    #     }
    st <- colnames(d)[j]
    h <- which(df$Station==st)
    co= c(df$Lat[h], df$Long[h])
    if (grepl(pattern = 'DCM', st)){
      lower.half.circle(co[2],co[1],col= colors[j], r=3)
    } else if (grepl(pattern = 'SUR',st)){
      upper.half.circle(co[2],co[1],col= colors[j], r=3)
    }
    #points(x=c[2], y=c[1], col=NA,bg= colors[j] ,pch=21,cex=20)
  }
  
  dev.off()
  #rgb_map <- list(e1_var, e2_var, e3_var)
}
