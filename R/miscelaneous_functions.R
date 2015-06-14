
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

dist_dp <- function(det_prob, range_test) {
  # Calculates the distance that correspond to a specific detection probability
  # in a logistic function
  
  rt <- range_test
  if (is.null(rt)) stop("The chosen method is 'model weighted' but 'model' was not provided")
  stopifnot(inherits(rt, "range_test"))
  
  co <- coefficients(rt)  # Get the coeficients
  names(co) <- NULL
  
  dist <- (qlogis(det_prob) - co[1]) / co[2]  # Invert the formula
  dist[dist < 1] <- 1  # Only return positive distances...
  return(dist)
}

# Determine if all values in a vector are the same
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
