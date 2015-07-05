
ssm_data <- function (detections, acoustic_array, range_test, 
                      sentinel, delta_t, ssm_grid = ssm_grid(acoustic_array),
                      start_time = NULL, end_time = NULL){
  
  n_station <- dplyr::n_distinct(acoustic_array$station)
  x_station <- (acoustic_array$lon - ssm_grid$lon_span[1]) * ssm_grid$lon2m
  y_station <- (acoustic_array$lat - ssm_grid$lat_span[1]) * ssm_grid$lat2m
  
  # if no start and end time are provided use the first and last detection
  if (is.null(start_time)){
    start_time <- min(detections$time)
  } else if (start_time > min(detections$time)) {
    stop("start_time is larger than the first detection")
  }
  if (is.null(end_time)){
    end_time <- max(detections$time)
  } else if (end_time < max(detections$time)) {
    stop("end_time is smaller than the last detection")
  }
  
  time_at_liberty <- as.numeric(difftime(end_time, start_time, units = "secs"))
  # adjust the end time so that it's a multiple of delta_t
  time_at_liberty <- ceiling(time_at_liberty / delta_t) * delta_t
  
  det_arr <- acoustic_array %>%
    dplyr::mutate(n_row = 1:nrow(acoustic_array)) %>%
    dplyr::select(station, n_row) %>%
    dplyr::inner_join(detections, by = "station") %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(secs = as.numeric(difftime(time, start_time, units = "secs")),
                  approx_secs = round(secs / delta_t) * delta_t) %>%
    dplyr::select(n_row, approx_secs) %>%
    dplyr::distinct() 
  
  n_det <- nrow(det_arr)
  
  rec_n <- det_arr$n_row
  rec_t <- det_arr$approx_secs
  
  ref_tag <- start_time  %>% 
    add(seq(0, time_at_liberty, by = delta_t)) %>%
    approx(sentinel$ts$time, sentinel$ts$efficiency, .) %>%
    magrittr::extract("y")
  
  data <- list(P0 = range_test$param$P0,
               D50 = range_test$param$D50,
               D95 = range_test$param$D95,
               dist_ref = sentinel$reference_distance,
               n_station = n_station, 
               x_station = x_station,
               y_station = y_station,
               delta_t = delta_t,
               time_at_liberty = time_at_liberty,
               n_det = n_det,
               rec_n = rec_n,
               rec_t = rec_t,
               ref_tag = ref_tag$y)
  
  class(data) <- append(class(data), "ssm_data")
  data
}

write_ssm_data <- function(ssm_data,
                           file = file.path (tempdir(), "data.dat")){
  
  file.create (file)
  # check that the input object is an ssm_grid
  if(!inherits(ssm_data, "ssm_data")) stop("")
  
  # data file
  write ("# State Space Model of passive telemetry movement", file = file)
  write ("# pmax", file, append = TRUE)
  write (ssm_data$P0, file, append = TRUE)
  write ("# D50", file, append = TRUE)
  write (ssm_data$D50, file, append = TRUE)
  write ("# D95", file, append = TRUE)
  write (ssm_data$D95, file, append = TRUE)
  write ("# dref", file, append = TRUE)
  write (ssm_data$dist_ref, file, append = TRUE)
  write ("# number of receivers", file, append = TRUE)
  write (ssm_data$n_station, file, append = TRUE)
  write ("# x", file, append = TRUE)
  write (ssm_data$x_station, file, append = TRUE, 
         ncolumns = length(ssm_data$x_station))
  write ("# y", file, append = TRUE)
  write (ssm_data$y_station, file, append = TRUE, 
         ncolumns = length(ssm_data$y_station))
  write ("# Transmitter data", file, append = TRUE)
  write ("# dt", file, append = TRUE)
  write (ssm_data$delta_t, file, append = TRUE)
  write ("# T (total time at liberty)", file, append = TRUE)
  write (ssm_data$time_at_liberty, file, append = TRUE)
  write ("# Number of observations", file, append = TRUE)
  write (ssm_data$n_det, file, append = TRUE)
  write ("# Reception receiver number", file, append = TRUE)
  write (ssm_data$rec_n, file, append = TRUE, ncolumns = length(ssm_data$rec_n))
  write ("# Reception times", file, append = TRUE)
  write (ssm_data$rec_t, file, append = TRUE, ncolumns = length(ssm_data$rec_t))
  write ("# Covariate from reference tag", file, append = TRUE)
  write (ssm_data$ref_tag, file, append = TRUE, ncolumns = length (ssm_data$ref_tag))
}