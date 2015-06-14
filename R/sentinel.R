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