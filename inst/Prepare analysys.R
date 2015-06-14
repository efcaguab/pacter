# Load fish data
load (file = "./Processed Data/Detections.RData")
load (file = "./Processed Data/Sentinel.RData")
levels (DP.reef$Reef) <- c("Malathu", "AlJadir")


library (doParallel)
library (reshape2)
library (ggplot2)
library (plyr)
library (dplyr)
library (doMC)
library (data.table)
library (lubridate)
registerDoMC (cores = 31)

det <- tbl_df(det)
# We will divide the analysis by reef. So if a fish moved between reefs we will only analyse stuff from one reef at a time. We'll create new levels that contain FishID and Reef at the same time
stations <- read.csv ("./Raw Data/Metadata/Stations_20140701.csv")
det <- merge(det, stations[, c(1, 4)])
det$FishReef <- with (det, paste(FishID, Reef))

# We first decide which fish and which species to do first. We'll chose the most
# data rich fish. For each fish and each day we count in how many stations it
# was detected. Then we add it all for each species. Species with larger numbers
# will have move more or have more individuals or have been monitored for a
# longer time
det$Day <- as.Date (det$Date)
nRec <- group_by (det, FishReef, Day) %>% summarise (Number = n_distinct(Station))

# Find out if the fish stayed at a single reef or moved around
nRec.total <- group_by (nRec, FishReef) %>% summarise (Total = sum (Number))
nRec.total <- mutate (nRec.total, FishID = unlist(lapply(strsplit(FishReef, " "), function (x) x[[1]])))
tagged.fish <- tbl_df(tagged.fish)
nRec.total <- merge(nRec.total, select (tagged.fish, FishID, Species), by = "FishID")
nRec.species <- group_by (nRec.total, Species) %>% summarise (TotalSp = sum(Total))
nRec.total <- merge (nRec.total, nRec.species) %>% arrange (desc (TotalSp), desc (Total))

det <- tbl_df (det)

# We will only analyse data using the state space model for Malathu and AlJadir and for fish detected by more than one receiver in a given day
nRec$FishID <-  unlist(lapply(strsplit(nRec$FishReef, " "), function (x) x[1]))
nRec$Reef <-  unlist(lapply(strsplit(nRec$FishReef, " "), function (x) x[2]))
ToAnalyze <- filter (nRec, Reef == "Malathu" | Reef == "AlJadir")

# Create directory structure ----------------------------------------------

# base.dir <- "../../../../../admb/fish_models/"
# # Relevel the detections fish factor
# fish.folders <- NA
# fish.date.folders <- vector ("list", nrow(fish.day))
# # Detections per fish per day
# un <- unique(ToAnalyze$FishReef)
# d_ply(ToAnalyze, .(FishReef), function (x, base.dir){
#   fish.folder <- paste (base.dir, as.character(x$FishReef[1]), "/", sep = "")
#   dir.create (fish.folder)
#   d_ply (x, .(Day), function (y, fish.folder){
#     fish.date.folder <- paste (fish.folder, "/", y$FishID[1], "_", y$Day[1], "/", sep = "")
#     dir.create (fish.date.folder)
#   }, fish.folder = fish.folder)
# }, base.dir = base.dir, .progress = "text")

# # Exploration -------------------------------------------------------------
# 
# # Check which fish has been detected in one place only
# receiver.fish <- ddply (det, "FishID", function(x) Stations <- length (unique (x$Station)))
# names (receiver.fish)[2] <- "Stations"
# receiver.fish$Stations <- as.factor (receiver.fish$Stations)
# receiver.fish <- dlply (receiver.fish, "Stations", function(x) x)
# 
# # Check the detection patterns by species
# # Include species in the detection DF
# det <- merge (det, tagged.fish[, c(1,6)], by = "FishID")
# species.fish <- dlply (tagged.fish, "Species", function (x) x$FishID)
# 
# registerDoParallel(cores = 12)
# foreach (i = 1:length(species.fish)) %dopar% {
#   if (nrow(det[det$FishID %in% species.fish[[i]], ]) > 0){
#     png (paste ("./Plots/PerSpecies/DotPlot ", as.character(names(species.fish[i])), ".png", sep = ""), width = 1600, height = 1050)
#       print (ggplot(det[det$FishID %in% species.fish[[i]], ], aes (x = Date, y = Station)) + geom_point(size = 1) + facet_wrap (~ FishID, scales = "free_y") + ggtitle(as.character(names(species.fish[i]))))
#     dev.off()
#     }
# }

# Sort fish that will be processed first ----------------------------------

# INITIALIZE ARRAY INFORMATION 
# Load range test information. Shape corresponds to average topographic and environmental conditions for daytime (the best range)
# load (file = "./Raw Data/Range Test/RangeTestResults.RData")
# RT <- data.frame (DIST = 1:length(DM), DP = DM)
# RT.Ped <- nls (DP ~ pmax / (1 + exp (log (19) * (DIST - D50) / (D95 - D50))),
#                data = RT, start = list (pmax = 0.86, D50 = 200, D95 = 400), control = nls.control(maxiter = 100), trace = TRUE)

## PARAMETERS THAT ARE THE SAME FOR ALL MODELS
# Found that pmax = 1, D50 = 193.2201, D95 = 360.79
pmax <- 1
D50 <- 193.2201
D95 <- 360.79
base.dir <- "~/admb/fish_models/"

# Transmitter data
# delta t
dt <- 240
# Tal - Time at liberty
Tal <- 24*3600 # One day in seconds

# Create interval for station locations
station.events <- read.csv ("./Raw Data/Metadata/Mooring assignments 22.06.14.csv")

# Reorganize data frame
station.events <- data.frame (Receiver = paste ("VR2W", station.events$Rec_SN, sep = "-"),
                              Station = station.events$Mooring,
                              Date.in = as.POSIXct (paste (station.events$Date_in, station.events$Time_in), format = "%d.%m.%y %H:%M", tz = "Asia/Riyadh"),
                              Date.out = as.POSIXct (paste (station.events$Date_out, station.events$Time_out), format = "%d.%m.%y %H:%M", tz = "Asia/Riyadh")) %>%
  mutate (period = new_interval (Date.in %>% floor_date (unit = "day"), Date.out %>% ceiling_date (unit = "day"))) # Create interval


# Cycle through each fish
foreach (i=1:nrow (nRec.total)) %dopar% {
  # Selecr days to be analysed for that fish
  days.fish <- filter (ToAnalyze, FishReef == nRec.total$FishReef[i])
  
  # Create fish directory 
  fish.folder <- paste (base.dir, as.character(days.fish$FishID[1]), "/", sep = "")
  dir.create (fish.folder)
  
  # Cycle trough each day
  foreach (j=1:nrow (days.fish)) %do% {
#  foreach (j=1:6) %dopar% {
    # Select detections that correspond to that fish and that day
    det.fish.day <- filter (det, FishReef == nRec.total$FishReef[i], as.Date (Date) == days.fish$Day[j]) %>%
      arrange (Date)
    # Determine reef
    reef <- unique (det.fish.day$Reef)
    # If the fish was detected in two or no reef
    if (length (reef) != 1) {
      message ("Wrong Reef for ", nRec.total$FishID[i], " at ", days.fish$Day[j])
      return (1)
      # Otherise chech which reef it is and asign corresponding bounding box
    } else if (reef == 'Malathu'){
      boundlon <- c(39.900193, 39.915900)
      boundlat <- c(19.743328, 19.756033)
    } else if (reef == 'AlJadir') {
      boundlon <- c( 39.945928, 39.960731)
      boundlat <- c(19.783289, 19.793611)
    } else {
      message ("Wrong Reef Name")
      return (1)
    }
    
    # Fill other parameters
    xmin <- 0
    ymin <- 0
    latkm <- 111*1000
    lonkm <- latkm*cos(mean(boundlat)*pi/180) 
    xmax <- (boundlon[2]-boundlon[1])*lonkm
    ymax <- (boundlat[2]-boundlat[1])*latkm
    cellwidth <- 50
    ngx <- floor (xmax/cellwidth) 
    ngy <- floor (ymax/cellwidth)
    gx <- seq(xmin,xmax,length=ngx)
    gy <- seq(ymin,ymax,length=ngy)
    
    # Find out which receivers were in 
    station.events$present <- (first (det.fish.day$Date) %>% floor_date (unit = "day")) %within% station.events$period
    stations.present.day <- tbl_df (station.events) %>% select (Station, present) %>% group_by (Station) %>% summarise (present = any (present))
    stations.day <- left_join (stations, stations.present.day) %>%
      filter (Reef == reef, present == TRUE) %>%
      mutate (Station = factor (Station))
    
    # Number of receivers
    nrec <- nrow (stations.day)
    
    # coordinates of the stations asssuming that all of them are available 
    x <- (stations.day$Lon - boundlon[1])*lonkm
    y <- (stations.day$Lat - boundlat[1])*latkm
    
    # Number of observations
    nobs <- nrow (det.fish.day)
    # Reception receiver number
    recn <- match (det.fish.day$Station, stations.day$Station)
    # Reception times
    times <- (((det.fish.day$Date %>% as.numeric ()) - 
                 (det.fish.day$Date %>% first() %>% floor_date("day") %>% as.numeric ())) / 240) %>%
      round() * 240
    
    # Reference tag covariates 
    sent.day <- filter (DP.reef, as.character (Reef) == as.character (reef), abs (difftime (floor_date (Date, "day"), first (det.fish.day$Date) %>% floor_date("day"), units = "secs")) <= 24*3600)
    sentinel.available <- (det.fish.day$Date %>% first () %>% floor_date(unit="day")) %in%
      (DP.reef$Date %>% floor_date (unit = "day"))
    # Create covariate DF
    covar <- data.frame (Date = seq (first (det.fish.day$Date) %>% floor_date ("day"), first (det.fish.day$Date) %>% floor_date ("day") + 3600*24, by = "240 sec"))
    
    # If there is not sentinel tag data
    if (!sentinel.available) {
      covar <- mutate (covar, seg = as.numeric (Date - min(Date)))
      covar$covariate[with (covar, (seg > 17700 & seg < 28500) | (seg > 60900 & seg < 71700))] <- 0.85
      covar$covariate[with (covar, seg >= 28500 & seg <= 60900)] <- 1
      covar$covariate[with (covar, is.na (covariate))] <- 0.6
    } else { # If there is, use it
      covar$covariate <- approx (sent.day$Date, sent.day$Index, covar$Date)$y
    }
    
    # Covariate from reference tags
    covref <- covar$covariate
    
    # Create Date directory
    fish.date.folder <- paste (fish.folder, "/", det.fish.day$FishID[1], "_", det.fish.day$Day[1], "/", sep = "")
    dir.create (fish.date.folder)
    # Copy tpl file and create the .dat and .cfg file
    file.copy ("./CodePedersen/onssm.tpl", fish.date.folder)
    data.file <- file.path (fish.date.folder, "data.dat")
    grid.file <- file.path (fish.date.folder, "grid.cfg")
    file.create (data.file)
    file.create (grid.file)
    
    # data file
    write (paste ("# Fish movement Red Sea - ", first(det.fish.day$FishReef), first (det.fish.day$Day)), file = data.file)
    write ("# pmax", data.file, append = TRUE)
    write (pmax, data.file, append = TRUE)
    write ("# D50", data.file, append = TRUE)
    write (D50, data.file, append = TRUE)
    write ("# D95", data.file, append = TRUE)
    write (D95, data.file, append = TRUE)
    write ("# dref", data.file, append = TRUE)
    write (D50, data.file, append = TRUE)
    write ("# number of receivers", data.file, append = TRUE)
    write (nrec, data.file, append = TRUE)
    write ("# x", data.file, append = TRUE)
    write (x, data.file, append = TRUE, ncolumns = length (x))
    write ("# y", data.file, append = TRUE)
    write (y, data.file, append = TRUE, ncolumns = length (y))
    write ("# Transmitter data", data.file, append = TRUE)
    write ("# dt", data.file, append = TRUE)
    write (dt, data.file, append = TRUE)
    write ("# T (total time at liberty)", data.file, append = TRUE)
    write (Tal, data.file, append = TRUE)
    write ("# Number of observations", data.file, append = TRUE)
    write (nobs, data.file, append = TRUE)
    write ("# Reception receiver number", data.file, append = TRUE)
    write (recn, data.file, append = TRUE, ncolumns = length (recn))
    write ("# Reception times", data.file, append = TRUE)
    write (times, data.file, append = TRUE, ncolumns = length (times))
    write ("# Covariate from reference tag", data.file, append = TRUE)
    write (covref, data.file, append = TRUE, ncolumns = length (covref))
    
    # grid file
    write ("# xmin", grid.file)
    write (xmin, grid.file, append = TRUE)
    write ("# ymin", grid.file, append = TRUE)
    write (ymin, grid.file, append = TRUE)
    write ("# xmax", grid.file, append = TRUE)
    write (xmax, grid.file, append = TRUE)
    write ("# ymax", grid.file, append = TRUE)
    write (ymax, grid.file, append = TRUE)
    write ("# ngx", grid.file, append = TRUE)
    write (ngx, grid.file, append = TRUE)
    write ("# ngy", grid.file, append = TRUE)
    write (ngy, grid.file, append = TRUE)
  }

}

save (nRec.total, ToAnalyze, file = "./Processed Data/MovementDataFrames.RData")





# 
# # Take out initial behaviour ----------------------------------------------
# 
# # To avoit taking into consideration anomalus behaviour due to the tagging
# # procedure we ignore one week after the surgery
# det <- ddply(det, "FishID", function(x){
#   x[as.Date (x$Date) > (as.Date (min (x$Date)) + 7), ] 
# })
# 
# # Remove "non-moving" fish ------------------------------------------------
# 
# # Beacuse it is not possible to quantify movement from fish that did not moved
# # at all we will remove fish that did not moved at all i.e. that stayed at a
# # single receiver for the whole duration of the experiment. We will also test
# # (by relating the detecton paterns to sentinel tag patterns) if fish was dead
# # or not
# 
# A <- ddply (det, "Transmitter", function(x) length (unique (x$Receiver)), .parallel = TRUE)
# A <- merge (A, tagged.fish[, c(1, 4, 6)])
# A[order(A$V1),]
# 
# ggplot(A) + geom_histogram (aes(x = V1, fill = Reef), binwidth = 1) + facet_wrap (~ Species)