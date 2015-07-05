#' Calculate center of activity
#' 
#' @description Calculates the center of activity in an observation network.
#'   
#' @param detections an object of the class 'detections' created with the
#'   function \code{\link{detections}}
#' @param acoustic_array an object of the class 'acoustic_array' created with
#'   the function \code{\link{acoustic_array}}
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
#' @param range_test an object of the class 'range_test' created with the
#'   function \code{\link{range_test}} containing the range test information
#' @param det_lon,det_lat numeric vectors containing the longitude of each 
#'   detection specified in 'det_stations' and 'det_times'
#' @param ... other arguments passed to \code{coa}
#'   
#' @return A Data frame
#' @details
#' 
#' Three methods for center of activity are providedwere evaluated for
#' estimating a tag’s position at time t when detected by two or more receivers
#' during a time interval (deltat). The first method (termed ‘Average’) was to
#' compute the mean position (latitude, longitude) of the detecting receivers.
#' The second was [7] method of using the mean of receiver locations weighted by
#' the number of detections (‘Observation-Weighted’, Table 6). The third
#' modified the ‘Observation-Weighted’ method to incorporate regression model
#' estimates of tag distance from a receiver (‘Model-Weighted’, Table 6). The
#' weighting term was computed as the ratio of estimated distance between the
#' tag and detecting receiver i (d i ) relative to the farthest estimated
#' distance of the tag from any of the detecting receivers (d max ). The
#' estimated distance (d i ) is determined by inputting the observed detection
#' rate into the logistic function. Following [7], both arithmetic and harmonic
#' mean estimators for each method were evaluated.
#' 
#' All receiver ids/station names in this vector must have been defined in the
#' station_names argument
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
#' 
sentinel <- function(times, efficiency, reference_distance, data){
  
  if(!missing(data)){
    if(all(is.character(substitute(times)), is.character(substitute(efficiency)))){
      times <- data[, times]
      efficiency <- data[, efficiency]
    }else if(all(is.name(substitute(times)), is.name(substitute(efficiency)))){
      times <- eval(substitute(times), data)
      efficiency <- eval(substitute(efficiency), data)
    } else {
      stop ("All detection efficiency data must be specified either as characters or names")
    }
  }
  
  if(!any(class(times) == "POSIXct")) {
    warning("'times' is not a POSIXct object. Coervertion will at UTC timezone will be attempted")
    times <- readr::parse_datetime(times)
    if(any(is.na(times))) warning("Problems parsing date-times detected")
  }
  
  # put into a data frame
  sent <- dplyr::data_frame(efficiency = efficiency, time = times)

  # Check if there is NA's
  if(sent %>% lapply(function(x) any(is.na(x))) %>% unlist() %>% any())
    warning("Missing or invalid values detected in 'efficiency' or 'times'. NA's will be removed" )
  
  # remove NAs
  sent %<>% 
    dplyr::filter(!is.na(efficiency), !is.na(time))
  
  sent <- list(ts = sent, 
               reference_distance = reference_distance)
  
  class(sent) <- append(class(sent), "sentinel")
  sent
}