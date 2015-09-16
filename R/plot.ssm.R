#' Title
#'
#' @param ssm_results 
#'
#'
plotssm <- function(ssm_results, animation = F){
  
  ngx <- ssm_results$state_space_grid$ngx
  ngy <- ssm_results$state_space_grid$ngy
  xmin <- ssm_results$state_space_grid$xmin
  xmax <- ssm_results$state_space_grid$xmax
  ymin <- ssm_results$state_space_grid$ymin
  ymax <- ssm_results$state_space_grid$ymax
  gx <- seq(xmin,xmax,length=ngx)
  gy <- seq(ymin,ymax,length=ngy)
  lonkm <- ssm_results$state_space_grid$lon2m
  latkm <- ssm_results$state_space_grid$lat2m
  boundlon <- ssm_results$state_space_grid$lon_span
  boundlat <- ssm_results$state_space_grid$lat_span
  lon <- gx/lonkm+boundlon[1]
  lat <- gy/latkm+boundlat[1]
  
  smoo <- ssm_results$positions$smoo_prob
  meanlon <- ssm_results$positions$meanlon
  meanlat <- ssm_results$positions$meanlat
  nt <- ssm_results$positions$nt
  
  xrec <- ssm_results$state_space_data$x_station
  yrec <- ssm_results$state_space_data$y_station
  lonstn <- xrec/lonkm+boundlon[1]
  latstn <- yrec/latkm+boundlat[1]
  
  lvls <- c(0,0.01,0.05,0.5,0.95,1)
  nlvl <- length(lvls)-1
  col <- gray(seq(1,0.2,length=nlvl))
  w <- 0.0004
  h <- 0.00015
  
  if(animation){
    graphics.off()
    dev.new(height=5,width=8)
    par(las=1)
    par(xpd=NA)
    for(i in 1:nt){
      filled.contour3(lon,lat,smoo[i,,],col=col,levels=lvls,xlab='Longitude',ylab='Latitude',cex.names=0.2)
      # points(lonstn,latstn,pch=22,bg='black')
      par(xpd=TRUE)
      #     for(j in 1:(nlvl-1)){
      #       rect(ll[1]+(j-1)*w,ll[2],ll[1]+j*w,ll[2]+h,col=col[j+1])
      #       text(ll[1]+(j-1)*w+0.5*w,ll[2]+1.7*h,as.character(1-lvls[j+1]))
      #     }
      points(meanlon[i],meanlat[i],pch='.',cex=2)
      inds <- max(c(1,i-100)):i
      lines(meanlon[inds],meanlat[inds])
      # lines(xscl,yscl)
      #     points(xscl,yscl,pch=3)
      # text(mean(xscl),yscl[1]+0.0002,'100 m')
      legend('topleft',legend=c('Stations','40 min tail'),pch=c(22,-1),col=c(1,1),lty=c(-1,1),pt.bg=1)
      # text(ll[1],ll[2]+3.2*h,'Confidence region',pos=4,cex=1.3,offset=0)
      # text(-162.1261,ll[2]+3*h,'Reef fish, Palmyra atoll',pos=4,cex=1.4,offset=0)
      Sys.sleep(0.1)
    }
    
  }
  
  reservoir <- qmap (location = c (lon = 1.820509, lat = 45.564856),
                     zoom = 14, color = "bw", legend = "topleft")
  
  plot <- reservoir +
    #     geom_point (data = data.frame(longitude = meanlon,
    #                                   latitude = meanlat), 
    #                 aes (x = longitude,
    #                      y = latitude, colour = hours (date_time))) +
    geom_segment (data = data.frame(longitude = meanlon,
                                    latitude = meanlat), 
                  aes (x = longitude,
                       y = latitude,
                       xend = lead (longitude),
                       yend = lead (latitude))) +
    geom_point (data = data.frame(longitude = meanlon,
                                    latitude = meanlat), 
                  aes (x = longitude,
                       y = latitude,
                       xend = lead (longitude),
                       yend = lead (latitude)))
  print(plot)
  
}

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
