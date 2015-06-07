#' Create objects of class `detections`
#' 
#' @description This function creates an object of the class `detections` to be 
#'   passed as argument of other functions that require information about the 
#'   location of the stations. It ensures that the formating and structure of 
#'   the detection data is consistent and makes it compatible with all other 
#'   functions in this package
#'   
#' @param det_stations a vector with the names of the receiver id/stations name 
#'   where the tag was detected. Alternatively, if the argument `data` is 
#'   provided it can be the name (quoted or unquoted) of the corresponding 
#'   column in data
#' @param det_times vector containing the times of the detections. If the vector
#'   is not an Date-Time (POSIXct) object, it would be coerced to POSIXct using 
#'   \code{\link[readr]{parse_datetime}}. Alternatively, if the argument `data` 
#'   is provided it can be the name (quoted or unquoted) of the corresponding 
#'   column in data
#' @param det_tag an (optional) vector containing the transmitter or animal ids.
#'   Alternatively, if the argument `data` is provided it can be the name 
#'   (quoted or unquoted) of the corresponding column in data
#' @param det_sensor an (optional) vector containing sensor information. 
#'   Alternatively, if the argument `data` is provided it can be the name 
#'   (quoted or unquoted) of the corresponding column in data
#' @param an (optional) data frame containing the detection information. If 
#'   provided, det_stations, det_times, det_tag and det_sensor must correspond 
#'   to (quoted or unquoted) column names
#'   
#' @return A `detections` object. A special case of a `data.frame` with standard
#'   columns: `station`, `time`, `tag`, `sensor`
#'   
#' @examples
#' plot_crayons()
#' 
#' @export
detections <- function(det_times, det_stations,
                       det_tag = NULL, det_sensor = NULL, data) {
  require(magrittr)
  
  if(!missing(data)){
    if(all(is.character(substitute(det_times)), is.character(substitute(det_stations)))){
      det_times <- data[, det_times]
      det_stations <- data[, det_stations]
      if(!is.null(det_tag)) det_tag <- data[, det_tag]
      if(!is.null(det_sensor)) det_sensor <- data[, det_sensor]
    }else if(all(is.name(substitute(det_times)), is.name(substitute(det_stations)))){
      det_times <- eval(substitute(det_times), data)
      det_stations <- eval(substitute(det_stations), data)
      if(!is.null(det_tag)) 
        det_tag <- eval(substitute(det_stations), data)
      if(!is.null(det_sensor)) 
        det_sensor <- eval(substitute(det_sensor), data)
    } else {
      stop ("All etection data must be specified either as characters or names")
    }
  }
  
  if(!any(class(det_times) == "POSIXct")) {
    warning("'det_times' is not a POSIXct object. Coervertion will at UTC timezone will be attempted")
    det_times <- readr::parse_datetime(det_times)
    if(any(is.na(det_times))) warning("Problems parsing date-times detected")
  }
  
  # put into a data frame
  if(is.null(det_sensor)){
    if(is.null(det_tag)){
      detections <- dplyr::data_frame(station = det_stations, time = det_times)
      warning ("Tag id's were no provided")
    } else {
      detections <- dplyr::data_frame(station = det_stations, 
                                      time = det_times, 
                                      tag = det_tag)
    }
  } else {
    detections <- dplyr::data_frame(station = det_stations, 
                                    time = det_times, 
                                    tag = det_tag,
                                    sensor = det_sensor)
  }
  
  # Check if there is NA's
  if(detections %>% lapply(function(x) any(is.na(x))) %>% unlist() %>% any())
    warning("Missing or invalid values detected in 'det_stations' or 'det_times'. NA's will be removed" )
  
  # remove NAs
  detections %<>% 
    dplyr::filter(!is.na(station), !is.na(time))
  
  # assign the detections class
  structure(detections, class = c("data.frame", "detections"))
}
