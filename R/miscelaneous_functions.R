get_positions <- function (station_names, station_positions) {
  # Get positions in a standard way given possible input formats, and check that
  # everything is correct
  
  # Standarise station positions
  station_positions <- std_positions(station_positions) 
  # Check that there is only one location per receiver 
  if(length(station_names) != nrow(station_positions))
    stop ("Number of stations in 'station_names' differ to the number of coordinates in 'station_positions'")
  
  positions <- data_frame (rec = station_names) %>%
    bind_cols(station_positions)
  
  # If there are NA's
  if(positions %>% 
       lapply(function(x) any(is.na(x))) %>%
       unlist() %>% any()) stop("One or more of the coordinates or station names is missing or is not a valid value")
  
  return(positions)
}


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
  
  return(out_pos)
}

