#' Create objects of class `range_test`
#' 
#' @description jhjh 
#' 
#' @param x logistic relationship of the detection probability as a function of
#'   distance from a receiver. It can be specified as 1) a named list with the 
#'   fields P0 (detection probability at a distance 0), D50 (distance 
#'   corresponding to a detection probability of 50\%), and D95 (distance 
#'   corresponding to a detection probability of 5\%). OR 2) a numeric vector
#'   with three elements corresponding to P0, D50 and D95. OR 3) a
#'   \code{\link{glm}} model of with a binomial family
#'   
#' @return j
#' @export
#' 
range_test <- function(x){
  UseMethod("range_test")
}
#' @export
range_test.list <- function(range_test){
  
  rt <- range_test
  if(all(c("P0", "D50", "D95") %in% names(rt))) {
    rt <- c(rt[["P0"]], rt[["D50"]], rt[["D95"]])
  } else {
    warning("Elements in 'model' list are not named default order (P0, D50, D95) was assumed")
    rt <- c(rt[[1]], rt[[2]], rt[[3]])
  }
  
  range_test(rt)
}
#' @export
range_test.numeric <- function(range_test){
  
  rt <- range_test
  
  if (length(rt) != 3) stop ("Incorrect model specification. It should contain three elements")
  data = data.frame(x = c(0, rt[2], rt[3]),
                    y = c(rt[1], 0.5, 0.05))
  
  if(data$y[1] > 1 | data$y[1] < 0) stop("Detection probability at 0m is outside the interval [0,1]")
  if(data$x[3] < data$x[2]) stop("Distance at which detection probability is 0.05 is larger than the distance at which detection probability is 0.5")
  
  suppressWarnings(rto <- glm (y ~ x, data = data, family = "binomial"))
  
  rto$param <- list(P0 = rt[1], D50 = rt[2], D95 = rt[3])
  class(rto) <- append(class(rto), "range_test")
  rto
}
#' @export
range_test.glm <- function(range_test){
  
  rt <- range_test
  
  if(length(coefficients(rt)) != 2) stop("The model must have only one intercept and one dependent variable")
  if(family(rt)$family != "binomial") warning("family of 'model' is not binomial")
  if(family(rt)$link != "logit") warning("family of 'model' is not logit")
  
  co <- coefficients(rt)
  rt$param <- list(P0 = plogis(co[1] + 0),
                   D50 = (qlogis(0.5) - co[1]) / co[2],
                   D95 = (qlogis(0.05) - co[1]) / co[2])
  class(rt) <- append(class(rt), "range_test")
  rt
}



