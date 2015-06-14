# load("~/github/VPS-SSM/data/processed-data/detections_corr.RData")
# library(dplyr)
# library(magrittr)
# 
# dete <- detections(date_time, name, data = detections_corr)
# a_array <- detections_corr %>%
#   select(name, latitude, longitude) %>%
#   distinct() %$%
#   acoustic_array(name, cbind(longitude, latitude))
# rt <- range_test(c(1, 50, 200))
# state_space_grid <- ssm_grid(a_array)
# delta_t <- 240
# sent <- distinct(detections, time) %$%
#   sentinel(time, 0.7, 60)
# state_space_data <- ssm_data(dete, a_array, rt, sent, 240)

# # center-of-activity
# load("~/github/VPS-SSM/data/processed-data/detections_corr.RData")


# # 
# station_names <- stations$name
# station_positions <- stations %>% dplyr::select(latitude, longitude)
# station_positions_2 <- as.matrix(station_positions_1)
# station_positions_3 <- station_positions_2
# colnames(station_positions_3) <- NULL
# det_stations <- detections_corr$name
# det_times <- detections_corr$date_time
# # 
# # interval = "30 min"
# # 
# # nls(y ~ 1 / (1 + exp(-b * (x -c))),
# #     data = data.frame(y = c(1, 0.5, 0.05), x = c(0,200, 500)))
# # 
# # model_1 <- glm (y ~ x,
# #                 data = data.frame(y = c(1, 0.5, 0.05), x = c(0,200, 500)), family = "binomial")
# # 
# # model_2 <- list(P0 = 1, D50 = 200, D95 = 500)
# # model_4 <- list(1, D509 = 200, D95 = 500)
# # 
# # model_3 <- c(1,200,500)
