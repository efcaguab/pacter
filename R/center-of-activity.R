#' Calculate center of activity
#' 
#' @description
#' Calculates the center of activity positions
#' 
#' @param station_names vector containing the station names (or receiver ids) 
#'   per array
#' @param station_positions numeric matrix, data.frame or SpatialPoints object 
#'   with the coordinates (Longitude/Latitude) of the stations. Each row is a 
#'   station and there must be only one coordinate per station
#' @param det_stations vector with the names of the receiver id/stations name 
#'   where the tag was detected. All receiver ids/station names in this vector 
#'   must have been defined in the station_names argument
#' @param det_times vector containing the times of the detections. If the vector
#'   is not an Date-Time (POSIXct) object, it would be coerced to POSIXct using 
#'   \code{\link[readr]{parse_datetime}}
#' @param breaks a vector of cut points or number giving the number of 
#'   intervals which the detections are to be grouped OR an interval 
#'   specification. Intervals can be one of "sec", "min", "hour", "day", 
#'   "DSTday", "week", "month", "quarter" or "year", optionally preceded by an 
#'   integer and a space, or followed by "s"
#' @param method "observation-weighted", "model-weighted" or average. If 
#'   "model-weighted" is selected, a logistic function must be specified in the 
#'   model argument. See details section for more information
#' @param mean_type use "arithmetic" or "harmonic" estimator for the mean 
#'   position
#' @param model logistic relationship of the detection probability as a function
#'   of distance from a receiver. Must be specified as a named list with the 
#'   fields P0 (detection probability at a distance 0), D50 (distance 
#'   corresponding to a detection probability of 50\%), and D95 (distance 
#'   corresponding to a detection probability of 5\%)
#'   
#' @return 
#' None  
#' @details 
#' Three methods were evaluated for estimating a tag’s position at time
#' t when detected by two or more receivers during a time interval (deltat). The 
#' first method (termed ‘Average’) was to compute the mean position (latitude, 
#' longitude) of the detecting receivers. The second was [7] method of using the
#' mean of receiver locations weighted by the number of detections 
#' (‘Observation-Weighted’, Table 6). The third modified the 
#' ‘Observation-Weighted’ method to incorporate regression model estimates of 
#' tag distance from a receiver (‘Model-Weighted’, Table 6). The weighting term
#' was computed as the ratio of estimated distance between the tag and detecting
#' receiver i (d i ) relative to the farthest estimated distance of the tag from
#' any of the detecting receivers (d max ). The estimated distance (d i ) is
#' determined by inputting the observed detection rate into the logistic
#' function. Following [7], both arithmetic and harmonic mean estimators for
#' each method were evaluated.
#' 
#' @references
#' Farmer, N. a, Ault, J. S., Smith, S. G., & Franklin, E. C. (2013). Methods for assessment of short-term coral reef fish movements within an acoustic array. Movement Ecology, 1(1), 7. \url{http://doi.org/10.1186/2051-3933-1-7}
#' 
#' Simpfendorfer, C., & Heupel, M. (2002). Estimation of short-term centers of activity from an array of omnidirectional hydrophones and its use in studying animal movements. Canadian Journal of Fisheries and Aquatic Sciences, 32, 23–32. \url{http://doi.org/10.1139/F01-191}
#' 
#' @examples
#' plot_crayons()
#' 
#' @export
# 

coa <- function (station_names, 
                 station_positions,
                 det_stations,
                 det_times,
                 breaks = "5 min", 
                 method = c("observation-weighted", "model-weighted", "average"),
                 mean_type = c("arithmetic", "harmonic"),
                 model = NULL) {
  
  # Get positions
  positions <- get_positions(station_names, station_positions)
  detections <- get_detections(det_stations, det_times)
  
  # Check that every receiver in detections has a position
  if (!all(detections$rec %in% positions$rec)) 
    stop("Not all station-names/receiver-ids in 'det_stations' are present in 'station_names'")

  detections %<>% 
    mutate(breaks = cut (time, breaks)) %>% # Cut by the given breaks
    group_by(breaks, rec)
  
  if (method == "average"){
    detections %<>% summarise(weight = 1) %>% group_by()
  } else if (method == "observation-weighted") {
    detections %<>% summarise(weight = n()) %>% group_by()
  } else if (method == "model-weighted") {
    detections %<>% 
      summarise(det = n()) %>%
      group_by() %>%
      mutate(det_prob = det/max(det),
             distance = dist_dp(det_prob, model)) %>%
      group_by(breaks) %>%
      mutate(max_dist = max(distance)) %>%
      group_by() %>%
      mutate(weight = max_dist / distance)
  }
  
  detections %<>%
    inner_join(positions) %>%
    group_by(breaks) %>%
    summarise(lon = vers_mean(lon, weight, mean_type),
              lat = vers_mean(lat, weight, mean_type))
    
}
