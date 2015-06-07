# get_positions <- function (station_names, station_positions) {
#   # Get positions in a standard way given possible input formats, and check that
#   # everything is correct
#   
#   # Standarise station positions
#   station_positions <- std_positions(station_positions) 
#   # Check that there is only one location per receiver 
#   if(length(station_names) != nrow(station_positions))
#     stop ("Number of stations in 'station_names' differ to the number of coordinates in 'station_positions'")
#   
#   positions <- data_frame (rec = station_names) %>%
#     bind_cols(station_positions)
#   
#   # If there are NA's
#   if(positions %>% 
#        lapply(function(x) any(is.na(x))) %>%
#        unlist() %>% any()) stop("One or more of the coordinates or station names is missing or is not a valid value")
#   
#   return(positions)
# }
# 
# 



# get_detections <- function (det_stations, det_times) {
#   # Get detections in a standard way. Check that dates and formats are correct
#   
#   # Check that times are correct
#   if(!any(class(det_times) == "POSIXct")) {
#     warning("'det_times' is not a POSIXct object. Convertion with UTC timezone will be attempted")
#     det_times <- readr::parse_datetime(det_times)
#     if(any(is.na(det_times))) warning("Problems parsing date-times detected")
#   }
#   
#   detections <- data_frame(rec = det_stations, time = det_times)
#   
#   # Check if there is NA's
#   if(detections %>% lapply(function(x) any(is.na(x))) %>% unlist() %>% any())
#     warning("Missing or invalid values detected in 'det_stations' or 'det_times'. NA's will be removed" )
#   
#   detections %<>% 
#     filter(!is.na(rec), !is.na(time))
# }

std_model <- function(model) {
  
  # If it is a list make it a numeric vector
  if (any(class(model) == "list")) {
    if(all(c("P0", "D50", "D95") %in% names(model))) {
      model <- c(model[["P0"]], model[["D50"]], model[["D95"]])
    } else {
      warning("Elements in 'model' list are not named default order (P0, D50, D95) was assumed")
      model <- c(model[[1]], model[[2]], model[[3]])
    }
  }
  
  # If is a numeric vector make it a data frame
  if (any(class(model) == "numeric")){
    if (length(model) != 3) stop ("Incorrect model specification. It should contain three elements")
    data = data.frame(x = c(0, model[2], model[3]),
                      y = c(model[1], 0.5, 0.05))
    
    if(data$y[1] > 1 | data$y[1] < 0) stop("Detection probability at 0m is outside the interval [0,1]")
    if(data$x[3] < data$x[2]) stop("Distance at which detection probability is 0.05 is larger than the distance at which detection probability is 0.5")
    
    suppressWarnings(model <- glm (y ~ x, data = data, family = "binomial"))
  }
  
  # If it's a model check that it's all good
  if (any(class(model) == "glm")) {
    if(length(coefficients(model)) != 2) stop("The model must have only one intercept and one dependent variable")
    if(family(model)$family != "binomial") warning("family of 'model' is not binomial")
    if(family(model)$link != "logit") warning("family of 'model' is not logit")
  } else {
    stop("The type of the model provided is not valid")
  }
  
  return(model)
}

# Mathematic -------------------

vers_mean <- function(x, w, mean_type = "arithmetic") {
  # Calculates weighted arithmetic or harmonic means

  if(mean_type[1] == "harmonic") {
    y <- 1 / weighted.mean(1 / x, w)
  } else {
    y <- weighted.mean (x, w)
  } 
  y
}

dist_dp <- function(det_prob, model) {
  # Calculates the distance that correspond to a specific detection probability
  # in a logistic function
  
  if (is.null(model)) stop("The chosen method is 'model weighted' but 'model' was not provided")
  
  model <- std_model(model)  # Standarize model input
  co <- coefficients(model)  # Get the coeficients
  names(co) <- NULL
  
  dist <- (logit(det_prob) - co[1]) / co[2]  # Invert the formula
  dist[dist < 1] <- 1  # Only return positive distances...
  return(dist)
}

# Determine if all values in a vector are the same
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
