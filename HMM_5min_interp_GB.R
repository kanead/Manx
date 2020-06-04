#' HMM code for Manx Shearwater tracking study
#' Great Blasket at 5 mins
rm(list = ls())

#' load packages
library(tidyverse)
library(momentuHMM)
library(rgdal)

#' load in data file with appended chlorophyll data from Movebank
mydata <-
  read_csv("data/movebank/Manx Shearwater 5 min interval-2041876219758145264.csv",
           col_names = TRUE)
head(mydata)

#' select relevant columns
mydata <-
  mydata %>% select(
    `location-long`,
    `location-lat`,
    timestamp,
    comments,
    `MODIS Ocean Aqua OceanColor 4km Monthly Chlorophyll A (OCI)`
  )
head(mydata)

#' rename the columns
mydata <- rename(mydata, long = `location-long`)
mydata <- rename(mydata, lat = `location-lat`)
mydata <- rename(mydata, ID = comments)
mydata <-
  rename(mydata, chloro = `MODIS Ocean Aqua OceanColor 4km Monthly Chlorophyll A (OCI)`)
head(mydata)
summary(mydata$chloro)

#' split the track up by colony
#' first get the colony name on its own from the IDs
mydata$colony <- sub("\\_.*", "", mydata$ID)
head(mydata)

#' extract the GB birds
mydata_GB <- filter(mydata, colony == "GB")

#' remove chlorophyll values with NA
#' check how this affects the length of the data frame
length(mydata_GB$chloro) #' 19220
mydata_GB <- mydata_GB %>%
  dplyr::mutate(chloro = ifelse(is.na(chloro), 0, chloro))

mydata_GB <- group_by(mydata_GB, ID) %>%
  dplyr::mutate(first2 = min(which(chloro == 0 |
                                     row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  dplyr::select(-first2)
mydata_GB <- mydata_GB %>% filter(chloro > 0)
#' new dataframe should be shorter and have no NAs or 0 values
length(mydata_GB$chloro) #' 18823
summary(mydata_GB$chloro)

#' extract the date stamp to use later
time <- mydata_GB$timestamp

#' do the same with lat long
lon_lat <- select(mydata_GB, long, lat)
#' write.csv(lon_lat, "lon_lat_GB.csv", row.names = FALSE)

#' project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(mydata_GB[, 1:2],
                         proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <-
  spTransform(llcoord, CRS("+proj=utm +zone=29 ellps=WGS84")) # 29 = IRE or 30 = UK
#' add UTM locations to data frame
mydata_GB$x <- attr(utmcoord, "coords")[, 1]
mydata_GB$y <- attr(utmcoord, "coords")[, 2]

#' Prep data for HMM
mydata_GB <- mydata_GB %>% dplyr::select(ID, x, y, chloro)
mydata_GB <- data.frame(mydata_GB)
mydata_GB <-
  momentuHMM::prepData(
    data = mydata_GB,
    type = "UTM",
    coordNames = c("x", "y"),
    covNames = "chloro"
  )

#' plot the data using UTM
plot(mydata_GB$x, mydata_GB$y)

#' take a look at the histogram of step length to inform initial values
hist(mydata_GB$step)

# Indices of steps of length zero, used for zero mass value
whichzero <- which(mydata_GB$step == 0)
# Proportion of steps of length zero in the data set
length(mydata_GB) / nrow(mydata_GB)

#' take a look at the turning angle
hist(mydata_GB$angle, breaks = seq(-pi, pi, length = 15))

#' fit HMM - 3 state model
#' without covariate

stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <-
  list(
    step = c(
      2000,
      500,
      100,
      2000,
      500,
      100,
      0.000318759,
      0.000318759,
      0.000318759
    ),
    angle = c(40, 10, 10)
  )
# fit model
m1 <-
  momentuHMM::fitHMM(
    data = mydata_GB,
    nbStates = 3,
    dist = dist,
    Par0 = Par0_m1,
    estAngleMean = list(angle = FALSE),
    stateNames = stateNames,
    retryFits = 1
  )
m1
plot(m1, plotCI = TRUE)

#' fit HMM - 3 state model
#' with monthly chlorophyll covariate

formula <- ~ chloro
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m2 <-
  list(
    step = c(
      2000,
      500,
      100,
      2000,
      500,
      100,
      0.000318759,
      0.000318759,
      0.000318759
    ),
    angle = c(40, 10, 10),
    formula = formula
  )
m2 <-
  momentuHMM::fitHMM(
    data = mydata_GB,
    nbStates = 3,
    dist = dist,
    Par0 = Par0_m2,
    estAngleMean = list(angle = FALSE),
    stateNames = stateNames,
    formula = formula,
    retryFits = 1
  )
m2
plot(m2, plotCI = TRUE)

#' compare models
AIC(m2, m1)
AIC(m1) - AIC(m2)

#' Model      AIC
#' 1    m2 293498 # with chloro
#' 2    m1 293508 # without chloro

#' save models
saveRDS(m1, "hmm_models/GB_no_chloro.rds")
saveRDS(m2, "hmm_models/GB_with_chloro.rds")

#' extract states and plot them
states <- viterbi(m2)
m2$data$state <- states
ggplot(data = m2$data, aes(x, y, color = as.factor(state))) + geom_point()

#' export the model details including states as a csv
m2$data <- rename(m2$data, ID_burst = ID)
m2$data$ID <- sub("_[^_]+$", "", m2$data$ID_burst)
levels(as.factor(m2$data$ID))
head(m2$data)
m2$data$datetime <- time
gb_data_states <-
  m2$data %>% select(x, y, datetime, chloro, state, ID, ID_burst)
head(gb_data_states)
tail(gb_data_states)
write.csv(gb_data_states, file = "results/GB_5min_states.csv", row.names = F)

#' stick on long and lat
gb_data_states <- read_csv("results/GB_5min_states.csv")
gb_data_states$long <- lon_lat$long
gb_data_states$lat <- lon_lat$lat

#' check the plots to see if the coords match
filter(gb_data_states, ID_burst == "GB_EZ04713_138_0/0") %>% ggplot(., aes(x, y)) + geom_point()
filter(gb_data_states, ID_burst == "GB_EZ04713_138_0/0") %>% ggplot(., aes(long, lat)) + geom_point()

#' export it again
head(gb_data_states)
write.csv(gb_data_states, file = "results/GB_5min_states.csv", row.names = F)

#' load in model if needed
hmm <- readRDS("hmm_models/GB_with_chloro.rds")
hmm
plot(hmm,plotCI = TRUE, cex.lab = 1.2)

#' Plot stationary state probabilities
plotStationary(hmm, plotCI = TRUE, legend.pos = "topleft")

#' check the number of trips we're left with
hmm_data <- read_csv("results/GB_5min_states.csv")
hmm_data$year <- format(hmm_data$datetime, "%Y")
hmm_data$newID <- paste(hmm_data$ID, hmm_data$year, sep = "_")
length(levels(as.factor(hmm_data$newID))) # 19
levels <- levels(as.factor(hmm_data$newID))
table(hmm_data$newID)
write.csv(levels, "GB_levels_HMM.csv", row.names = F)
#' check sample size excluding resting
filter(hmm_data,state<3)
