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
#'   parse_datetime()
#' @param interval a vector of cut points or number giving the number of 
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
#' @examples
#' plot_crayons()
#' 
#' @export

# # center-of-activity
# load("./data/processed-data/detections_corr.RData")
# 
# stations <- detections_corr %>%
#   select(name, latitude, longitude) %>%
#   distinct()
# 
# rec <- stations$name
# rec_pos <- stations %>% select(latitude, longitude)
# 
# det <- detections_corr$name
# det_times <- detections_corr$date_time
# 
# interval = "30 min"

coa <- function (station_names, 
                                station_positions,
                                det_stations,
                                det_times,
                                interval = "5 min", 
                                method = c("observation-weighted", "model-weighted", "average"),
                                mean_type = c("arithmetic", "harmonic"),
                                model = NULL) {
   
  # Check that all detections have a receiver position ----------
  # Check that times are correct ----------
  if(any(class(det_times)) != "POSIXct") det_times <- parse_datetime(det_times)
  # Check that there is only one location per receiver ----------
  
  
  class(station_pos) <- "data.frame"
  
  detections <- data_frame(rec = det, time = det_times)
  positions <- data_frame(rec = rec, lat = rec_pos[, 1], lon = rec_pos[, 2])
  
  detections %<>% 
    inner_join(positions) %>%
    mutate(breaks = cut (time, interval))
  
  detections %<>% 
    group_by(breaks) %>%
    summarise(lat = mean(lat), 
              lon = mean(lon))
}