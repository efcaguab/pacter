
ssm_grid <- function(acoustic_array, 
                     lon_grid = list(from = NA, to = NA, extra = 0.2),
                     lat_grid = list(from = NA, to = NA, extra = 0.2), 
                     cell_width = NULL, 
                     n_cells = 2500){
  
  # fill lists if user specified only some of the options
  default_grid <- list(from = NA, to = NA, extra = 0.2)
  lon_grid <- ttutils:::merge.list(lon_grid, default_grid)
  lat_grid <- ttutils:::merge.list(lat_grid, default_grid)
  
  # Sanity checks
  stopifnot(is.numeric(unlist(lon_grid)),
            is.numeric(unlist(lat_grid)),
            is.null(cell_width) | is.numeric(cell_width),
            is.numeric(n_cells))
  
  if(n_cells > 10000) warning("High number of cells, computation might bee big for RAM memory")
  
  # Variables to convert latitude and longitude to meters
  lat2m <- 111*1000
  lon2m <- lat2m * cos(mean(acoustic_array$lat) * pi/180) 
  
  # Initial span based on the array location
  lon_span <- range(acoustic_array$lon)
  lat_span <- range(acoustic_array$lat)
  
  # Modify the span based on user imput or defaults
  lon_span <- fix_limit(lon_grid, lon_span)
  lat_span <- fix_limit(lat_grid, lat_span)
  
  xmin <- 0
  ymin <- 0
  xmax <- (lon_span[2] - lon_span[1]) * lon2m
  ymax <- (lat_span[2] - lat_span[1]) * lat2m
  
  if(is.null(cell_width)){
    cell_width <- ceiling((xmax + ymax) / round(sqrt(n_cells) * 2))
  }
  
  ngx <- floor(xmax / round(cell_width)) 
  ngy <- floor(ymax / round(cell_width))
  
  structure(list(xmin = xmin, 
                 ymin = ymin, 
                 xmax = xmax, 
                 ymax = ymax,
                 ngx = ngx, 
                 ngy = ngy, 
                 lon_span = lon_span,
                 lat_span = lat_span,
                 lat2m = lat2m,
                 lon2m = lon2m), 
            class = "ssm_grid")
}



# Modify limits based on user imput

fix_limit <- function (limits, span) {
  if(is.na(limits$from)){
    span[1] <- `-`(span[1], (span[2] - span[1]) * limits$extra)
  } else {
    span[1] <- limits$from
  }
  if(is.na(limits$to)){
    span[2] <- `+`(span[2], (span[2] - span[1]) * limits$extra)
  } else {
    span[2] <- limits$to
  }
  span
}

write_ssm_grid <- function(ssm_grid, 
                           file = file.path (tempdir(), "grid.cfg")){
  
  file.create (file)
  # check that the input object is an ssm_grid
  if(!inherits(ssm_grid, "ssm_grid")) stop("")
  
  # add contents file
  write("# xmin", file)
  write(ssm_grid$xmin, file, append = TRUE)
  write("# ymin", file, append = TRUE)
  write(ssm_grid$ymin, file, append = TRUE)
  write("# xmax", file, append = TRUE)
  write(ssm_grid$xmax, file, append = TRUE)
  write("# ymax", file, append = TRUE)
  write(ssm_grid$ymax, file, append = TRUE)
  write("# ngx", file, append = TRUE)
  write(ssm_grid$ngx, file, append = TRUE)
  write("# ngy", file, append = TRUE)
  write(ssm_grid$ngy, file, append = TRUE)
  
}
