#' Title
#'
#' @param ssm_results 
#' @param state_space_grid 
#'
#' @return
#' @export
#'
#' @examples
get_ssm_smooth <- function(ssm_results, state_space_grid){
  
  smoodat <- tail(ssm_results$coefficients, -6)
  ngx <- state_space_grid$ngx
  ngy <- state_space_grid$ngy
  xmin <- state_space_grid$xmin
  xmax <- state_space_grid$xmax
  ymin <- state_space_grid$ymin
  ymax <- state_space_grid$ymax
  gx <- seq(xmin,xmax,length=ngx)
  gy <- seq(ymin,ymax,length=ngy)
  lonkm <- state_space_grid$lon2m
  latkm <- state_space_grid$lat2m
  boundlon <- state_space_grid$lon_span
  boundlat <- state_space_grid$lat_span
  lon <- gx/lonkm+boundlon[1]
  lat <- gy/latkm+boundlat[1]
  
  nt <- length(smoodat)/(ngy * ngx)
  smoo <- array(0,dim=c(nt,ngx,ngy))
  smoounorm <- array(0,dim=c(nt,ngx,ngy))
  fac <- ngy*ngx
  
  for(i in 1:nt){
    ind <- ((i-1)*fac+1):(i*fac)
    smoounorm[i,,] <- matrix(smoodat[ind],ngx,ngy)
    smoounorm[i,,] <- smoounorm[i,,]/sum(smoounorm[i,,])
    smoo[i,,] <- normalize(smoounorm[i,,])
  }
  
  md <- apply(smoounorm,c(1,2),sum)
  mx <- repmat(t(as.matrix(gx)),nt,1)
  meanx <- apply(md*mx,1,sum)
  meanlon <- meanx/lonkm+boundlon[1]
  md <- apply(smoounorm,c(1,3),sum)
  mx <- repmat(t(as.matrix(gy)),nt,1)
  meany <- apply(md*mx,1,sum)
  meanlat <- meany/latkm+boundlat[1]
  
  return(list(meanlon = meanlon, 
              meanlat = meanlat,
              smoo_prob = smoo,
              nt = nt))
}

normalize <- function(distr){
  size <- dim(distr)
  rd <- as.vector(distr)
  rd <- rd/sum(rd)
  rds <- sort(rd,index.return=TRUE)
  rdsc <- cumsum(rds$x)
  rd[rds$ix] <- rdsc
  rd
}

## Track calculated from mean
repmat <- function(X,m,n){
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}