std_positions <- function (station_positions) {
  # Transform any of the possible station positions into a data frame with lat and lon 
  
  # If it's a tbl_df convert to a normal data.frame
  if(any(class(station_positions) == "data.frame")) 
    class(station_positions) <- "data.frame"
  
  # If it's a SpacialPoints object extract the coordinates
  if(any(class(station_positions) == "SpatialPoints"))
    station_positions <- coordinates(station_positions)
  
  # Check that it has at least two columns
  if(ncol(station_positions) < 2) stop("'station_positions' has less than two columns")
  
  # Check that there are columns called latitude/longitude
  if(all(c("lat", "lon") %in% colnames(station_positions))){
    out_pos <- data_frame(lon = station_positions[,"lon"],
                          lat = station_positions[, "lat"])
  } else if(all(c("latitude", "longitude") %in% colnames(station_positions))){
    out_pos <- data_frame(lon = station_positions[,"longitude"],
                          lat = station_positions[, "latitude"])
  } else {
    if(ncol(station_positions) > 2) 
      warning("'station_positions' has more than two columns, only the first two are used")
    out_pos <- data_frame(lon = station_positions[, 1],
                          lat = station_positions[, 2])
  }
  
  # Coerce to numeric
  if(!(is.numeric(out_pos$lat) & is.numeric(out_pos$lon))) {
    warning("Coordinates are not numbers, coercing to 'numeric'")
    out_pos %<>% mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  }
}

