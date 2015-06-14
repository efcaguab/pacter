## Supplementary material S2
## R script for displaying results of SSM analysis of observation network data
## 30.04.2013
rm(list=ls())
require(animation)

## Read grid information
gridfile <- 'grid.cfg'
xmin <- scan(gridfile,skip=1,nmax=1)
ymin <- scan(gridfile,skip=3,nmax=1)
xmax <- scan(gridfile,skip=5,nmax=1)
ymax <- scan(gridfile,skip=7,nmax=1)
ngx <- scan(gridfile,skip=9,nmax=1)
ngy <- scan(gridfile,skip=11,nmax=1)
gx <- seq(xmin,xmax,length=ngx)
gy <- seq(ymin,ymax,length=ngy)

## Bounds
boundlon <- c(-162.1255,-162.118)
boundlat <- c(5.8735,5.8765)
latkm <- 111*1000 ## One degree of lat in metres
lonkm <- latkm*cos(mean(boundlat)*pi/180)
lon <- gx/lonkm+boundlon[1]
lat <- gy/latkm+boundlat[1]

## Read smoo
scriptname <- 'onssm'
normalize <- function(distr){
  size <- dim(distr)
  rd <- as.vector(distr)
  rd <- rd/sum(rd)
  rds <- sort(rd,index.return=TRUE)
  rdsc <- cumsum(rds$x)
  rd[rds$ix] <- rdsc
  rd
}
smoodat <- scan(paste(scriptname,'.rep',sep=''))
nt <- length(smoodat)/(ngy*ngx)
smoo <- array(0,dim=c(nt,ngx,ngy))
smoounorm <- array(0,dim=c(nt,ngx,ngy))
fac <- ngy*ngx
for(i in 1:nt){
  ind <- ((i-1)*fac+1):(i*fac)
  smoounorm[i,,] <- matrix(smoodat[ind],ngx,ngy)
  smoounorm[i,,] <- smoounorm[i,,]/sum(smoounorm[i,,])
  smoo[i,,] <- normalize(smoounorm[i,,])
}

## Track calculated from mean
repmat <- function(X,m,n){
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}
md <- apply(smoounorm,c(1,2),sum)
mx <- repmat(t(as.matrix(gx)),nt,1)
meanx <- apply(md*mx,1,sum)
meanlon <- meanx/lonkm+boundlon[1]
md <- apply(smoounorm,c(1,3),sum)
mx <- repmat(t(as.matrix(gy)),nt,1)
meany <- apply(md*mx,1,sum)
meanlat <- meany/latkm+boundlat[1]

## Plot 2D animation of movements ##
filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
# RV - 10-03-2012
# note replacement of .Internal(filledcontour(as.double(x),...)
# with .filled.contour() as of R-2.15.0
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

xrec <- as.numeric(read.table('data.dat',skip=12,nrows=1))
yrec <- as.numeric(read.table('data.dat',skip=14,nrows=1))
lonstn <- xrec/lonkm+boundlon[1]
latstn <- yrec/latkm+boundlat[1]

lvls <- c(0,0.01,0.05,0.5,0.95,1)
nlvl <- length(lvls)-1
col <- gray(seq(1,0.2,length=nlvl))
w <- 0.0004
h <- 0.00015
ll <- c(-162.12,5.8766)
xscl <- c(0,100/lonkm)-162.1253
yscl <- 5.8755+c(0,0)
graphics.off()
dev.new(height=5,width=8)
par(las=1)
par(xpd=NA)
for(i in 1:nt){
  filled.contour3(lon,lat,smoo[i,,],col=col,levels=lvls,xlab='Longitude',ylab='Latitude',cex.names=0.2)
  points(lonstn,latstn,pch=22,bg='black')
  par(xpd=TRUE)
  for(j in 1:(nlvl-1)){
    rect(ll[1]+(j-1)*w,ll[2],ll[1]+j*w,ll[2]+h,col=col[j+1])
    text(ll[1]+(j-1)*w+0.5*w,ll[2]+1.7*h,as.character(1-lvls[j+1]))
  }
  points(meanlon[i],meanlat[i],pch='.',cex=2)
  inds <- max(c(1,i-10)):i
  lines(meanlon[inds],meanlat[inds])
  lines(xscl,yscl)
  points(xscl,yscl,pch=3)
  text(mean(xscl),yscl[1]+0.0002,'100 m')
  legend('topleft',legend=c('Stations','40 min tail'),pch=c(22,-1),col=c(1,1),lty=c(-1,1),pt.bg=1)
  text(ll[1],ll[2]+3.2*h,'Confidence region',pos=4,cex=1.3,offset=0)
  text(-162.1261,ll[2]+3*h,'Reef fish, Palmyra atoll',pos=4,cex=1.4,offset=0)
  Sys.sleep(0.1)
}
