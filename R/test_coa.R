# library(magrittr)
# library(dplyr)
# library(geosphere)
# load("~/github/VPS-SSM/data/processed-data/gps_positions_corr.RData")
# blo <- readRDS("~/github/VPS-SSM/data/processed-data/analysis_blocks.rds")
# snt <- readRDS("~/github/VPS-SSM/data/processed-data/sentinel_covariate.rds")
# det <- readRDS("~/github/VPS-SSM/data/processed-data/clumped_detections.rds")
# gps <- gps_positions_corr
# # SSM with negative information
# 
# #Let's have a look at how the model behaves when including the negative information
# 
# source("~/github/pacter/R/plot.ssm.R")
# 
# x <- blo[2,]
# 
# det %<>%
# filter(date_time > x$sta, date_time < x$end)
# 
# dete <- detections(det$date_time_avg, det$name)
# a_array <- det %>%
# select(name, latitude, longitude) %>%
# distinct() %$%
# acoustic_array(name, cbind(longitude, latitude))
# rt <- list(param = list(P0 = x$pmax, D50 = x$D50, D95 = x$D95))
# class(rt) <- "range_test"
# sent <- sentinel(snt$date_time_avg, snt$covariate, x$ref_dist)
# delta_t <- 25
# n_cells <- 1000
# R2admb::setup_admb()
# ssm_results <- ssm(dete, a_array, rt, 
#                    sent, 
#                    delta_t = delta_t, 
#                    verbose = T, 
#                    n_cells = n_cells,
#                    params = list(logbx = -9, 
#                                  logby = -9, 
#                                  logitcor = 0, 
#                                  logsigma = 3), 
#                    neggative_info = T)
# 
# pos <- data.frame(date_time = seq(0, ssm_results$state_space_data$time_at_liberty,
#                                   by = ssm_results$state_space_data$delta_t) +
#                     ssm_results$state_space_data$start_time, 
#                   longitude = ssm_results$positions$meanlon,
#                   latitude  = ssm_results$positions$meanlat)
# 
# exp_pos <- data_frame(date_time = seq(x$sta, x$end, by = "10 sec"), 
#                       ssm_lon = approx(pos$date_time, 
#                                        pos$longitude, 
#                                        date_time, 
#                                        rule = 2)$y,
#                       ssm_lat = approx(pos$date_time,
#                                        pos$latitude, 
#                                        date_time,
#                                        rule = 2)$y,
#                       gps_lon = approx(gps$date_time,
#                                        gps$longitude,
#                                        date_time,
#                                        rule = 2)$y,
#                       gps_lat = approx(gps$date_time,
#                                        gps$latitude,
#                                        date_time,
#                                        rule = 2)$y,
#                       dist = distMeeus(as.matrix(cbind(ssm_lon, ssm_lat)),
#                                        as.matrix(cbind(gps_lon, gps_lat))))
# 
# output <- data.frame(index = index,
#                      ME = sum(exp_pos$dist)/nrow(exp_pos),
#                      MSE = sum(exp_pos$dist^2)/nrow(exp_pos))
# 
# 
# 
# plotssm(ssm_results)
