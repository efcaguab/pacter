#' Create objects of class `acoustic_array`
#' 
#' @description Creates an object of the class `acoustic_array`. It ensures
#'   compatibility with other functions in the package
#'   
#' @param station_names vector containing the station names (or receiver ids) in
#'   the array
#' @param station_positions numeric matrix, data.frame or SpatialPoints object 
#'   with the coordinates (Longitude/Latitude) of the stations. Each row is a 
#'   station and there must be only one coordinate per station
#'   
#' @return An object of the class `acoustic array`. It is a special case of a 
#'   data.frame with columns 'station', 'lon' and 'lat'
#'   
#' @examples
#' plot_crayons()
#' 
#' @export
acoustic_array <- function(station_names, station_positions) {
  require(magrittr)
  
  # Get positions in a standard way given possible input formats, and check that
  # everything is correct
  
  # Standarise station positions
  station_positions <- std_positions(station_positions) 
  # Check that there is only one location per receiver 
  if(length(station_names) != nrow(station_positions))
    stop ("Number of stations in 'station_names' differ to the number of coordinates in 'station_positions'")
  
  positions <- dplyr::data_frame (station = station_names) %>%
    dplyr::bind_cols(station_positions)
  
  # If there are NA's
  if(positions %>% 
       lapply(function(x) any(is.na(x))) %>%
       unlist() %>% any()) stop("One or more of the coordinates or station names is missing or is not a valid value")
  
  structure(positions, class = c("data.frame", "acoustic_array")) 
}


std_positions <- function (station_positions) {
  # Transform any of the possible station positions into a data frame with lat
  # and lon
  
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
    out_pos <- dplyr::data_frame(lon = station_positions[,"lon"],
                                 lat = station_positions[, "lat"])
  } else if(all(c("latitude", "longitude") %in% colnames(station_positions))){
    out_pos <- dplyr::data_frame(lon = station_positions[,"longitude"],
                                 lat = station_positions[, "latitude"])
  } else {
    if(ncol(station_positions) > 2) 
      warning("'station_positions' has more than two columns, only the first two are used")
    out_pos <- dplyr::data_frame(lon = station_positions[, 1],
                                 lat = station_positions[, 2])
  }
  
  # Coerce to numeric
  if(!(is.numeric(out_pos$lat) & is.numeric(out_pos$lon))) {
    warning("Coordinates are not numbers, coercing to 'numeric'")
    out_pos %<>% mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  }
  
  out_pos
}