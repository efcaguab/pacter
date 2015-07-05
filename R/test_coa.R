# load("~/github/VPS-SSM/data/processed-data/detections_corr.RData")
# library(dplyr)
# library(magrittr)
# 
# detections_corr %<>%
#   filter(as.Date(date_time) == unique(as.Date(date_time))[2])
# 
# dete <- detections(date_time, name, data = detections_corr)
# a_array <- detections_corr %>%
#   select(name, latitude, longitude) %>%
#   distinct() %$%
#   acoustic_array(name, cbind(longitude, latitude))
# rt <- range_test(c(1, 50, 200))
# # state_space_grid <- pacter:::ssm_grid(a_array)
# sent <- distinct(dete, time) %$%
#   sentinel(time, 0.7, 60)
# 
# R2admb::setup_admb()
# a <- ssm(dete, a_array, rt, sent, delta_t = 30, verbose = TRUE, n_cells = 500)
# 
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
