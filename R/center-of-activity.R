#' Calculate center of activity
#' 
#' @description Calculates the center of activity in an observation network.
#'   
#' @param detections an object of the class 'detections' created with the function \code{\link{detections}}
#' @param acoustic_array an object of the class 'acoustic_array' created with the function \code{\link{acoustic_array}}

#' @param breaks a vector of cut points or number giving the number of intervals
#'   which the detections are to be grouped OR an interval specification. 
#'   Intervals can be one of "sec", "min", "hour", "day", "DSTday", "week", 
#'   "month", "quarter" or "year", optionally preceded by an integer and a 
#'   space, or followed by "s"
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
#' @param det_lon,det_lat numeric vectors containing the longitude of each
#'   detection specified in 'det_stations' and 'det_times'
#' @param ... other arguments passed to \code{coa}
#'   
#' @return 
#' A Data frame
#' @details 
#' 
#' Three methods for center of activity are providedwere evaluated for estimating a tag’s position at time
#'   t when detected by two or more receivers during a time interval (deltat).
#'   The first method (termed ‘Average’) was to compute the mean position
#'   (latitude, longitude) of the detecting receivers. The second was [7] method
#'   of using the mean of receiver locations weighted by the number of
#'   detections (‘Observation-Weighted’, Table 6). The third modified the 
#'   ‘Observation-Weighted’ method to incorporate regression model estimates of 
#'   tag distance from a receiver (‘Model-Weighted’, Table 6). The weighting
#'   term was computed as the ratio of estimated distance between the tag and
#'   detecting receiver i (d i ) relative to the farthest estimated distance of
#'   the tag from any of the detecting receivers (d max ). The estimated
#'   distance (d i ) is determined by inputting the observed detection rate into
#'   the logistic function. Following [7], both arithmetic and harmonic mean
#'   estimators for each method were evaluated.
#'   
#'   All receiver ids/station names in this vector 
#'   must have been defined in the station_names argument
#'   
#' @references Farmer, N. a, Ault, J. S., Smith, S. G., & Franklin, E. C. 
#'   (2013). Methods for assessment of short-term coral reef fish movements
#'   within an acoustic array. Movement Ecology, 1(1), 7. 
#'   \url{http://doi.org/10.1186/2051-3933-1-7}
#'   
#'   Simpfendorfer, C., & Heupel, M. (2002). Estimation of short-term centers of
#'   activity from an array of omnidirectional hydrophones and its use in
#'   studying animal movements. Canadian Journal of Fisheries and Aquatic
#'   Sciences, 32, 23–32. \url{http://doi.org/10.1139/F01-191}
#'   
#' @examples
#' plot_crayons()
#' 
#' @export
coa <- function (detections,
                 acoustic_array,
                 breaks = "5 min", 
                 method = c("observation-weighted", "model-weighted", "average"),
                 mean_type = c("arithmetic", "harmonic"),
                 model = NULL) {
  
  # Get positions
  stopifnot(inherits(detections, "detections"),
            inherits(acoustic_array, "acoustic_array"))
  
  
  # Check that every receiver in detections has a position
  if (!all(detections$station %in% acoustic_array$station)) 
    stop("Not all station-names/receiver-ids in 'det_stations' are present in 'station_names'")

  # load magrittr
  require(magrittr)
  
  detections %<>% 
    dplyr::mutate(breaks = cut (time, breaks)) %>% # Cut by the given breaks
    dplyr::group_by(breaks, station)
  
  if (method[1] == "average"){
    detections %<>% dplyr::summarise(weight = 1) %>% dplyr::group_by()
  } else if (method[1] == "observation-weighted") {
    detections %<>% dplyr::summarise(weight = n()) %>% dplyr::group_by()
  } else if (method[1] == "model-weighted") {
    detections %<>% 
      dplyr::summarise(det = n()) %>%
      dplyr::group_by() %>%
      dplyr::mutate(det_prob = det/max(det),
                    distance = dist_dp(det_prob, model)) %>%
      dplyr::group_by(breaks) %>%
      dplyr::mutate(max_dist = max(distance)) %>%
      dplyr::group_by() %>%
      dplyr::mutate(weight = max_dist / distance)
  }
  
  detections %<>%
    dplyr::inner_join(acoustic_array, by = "station") %>%
    dplyr::group_by(breaks) %>%
    dplyr::summarise(lon = vers_mean(lon, weight, mean_type),
                     lat = vers_mean(lat, weight, mean_type))
  
  structure(detections, class = c("data.frame", "coa", "position"))
}

#' @describeIn coa A wrapper for \code{coa()} in which there is one longitude/latitude
#'   pair for each detection
#' @export
coa_d <- function(det_stations, det_times, det_lon, det_lat, ...) {
  
  if(!zero_range(c(length(det_stations),
                  length(det_times),
                  length(det_lon),
                  length(det_lat)))) stop("'det_stations', 'det_times', 'det_lon', and 'det_lat' must have the same length")
  
  # Join in a data frame
  det <- data.frame (name = det_stations, 
                     time = det_times, 
                     lon = det_lon, 
                     lat = det_lat)
  
  # Get the location of the stations
  stations <- det %>% select (name, lon, lat) %>% distinct ()
  
  # Calculate the coa
  coa(det$name, det$time, stations$name, stations %>% select(lon, lat), ...)
}