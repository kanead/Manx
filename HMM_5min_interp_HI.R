#' HMM code for Manx Shearwater tracking study
#' High Island at 5mins
rm(list = ls())

#' load packages
library(tidyverse)
library(momentuHMM)
library(rgdal)

#' load in data file with appended chlorophyll data from Movebank
mydata <-
  read_csv(
    "data/movebank/Manx Shearwater 5 min interval-2041876219758145264.csv",
    col_names = TRUE
  )
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

#' extract the HI birds
mydata_HI <- filter(mydata, colony == "HI")

#' remove chlorophyll values with NA
#' check how this affects the length of the data frame
length(mydata_HI$chloro) #' 85660
mydata_HI <- mydata_HI %>%
  dplyr::mutate(chloro = ifelse(is.na(chloro), 0, chloro))

mydata_HI <- group_by(mydata_HI, ID) %>%
  dplyr::mutate(first2 = min(which(chloro == 0 |
                                     row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  dplyr::select(-first2)
mydata_HI <- mydata_HI %>% filter(chloro > 0)
#' new dataframe should be shorter and have no NAs or 0 values
length(mydata_HI$chloro) #' 80731
summary(mydata_HI$chloro)

#' the NA removal makes some IDs less than 3 observations long
#' Need to remove these for the HMM
#' Keep it consitent with the data prep by removing those with n < 10
#' how many IDs with fewer than 10 relocations?
which(table(mydata_HI$ID) < 10)

length(mydata_HI$chloro) 
mydata_HI <- mydata_HI %>%
  group_by(ID) %>%
  filter(n() > 10)  # require each level with more than 10 obs
length(mydata_HI$chloro) 

#' remove weird deployment to check if the effect of
#' chlorophyll is still there; it is (see below)
# levels(as.factor(mydata_HI$ID))
# mydata_HI <- mydata_HI %>% 
#  filter(!ID %in% "HI_FB48028_086_0/0") 

#' extract the date stamp to use later
time <- mydata_HI$timestamp

#' do the same with lat long
lon_lat <- select(mydata_HI, long, lat)
#' write.csv(lon_lat, "lon_lat_HI.csv", row.names = FALSE)

#' project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(mydata_HI[, 1:2],
                         proj4string = CRS("+proj=longlat +datum=WGS84"))
utmcoord <-
  spTransform(llcoord, CRS("+proj=utm +zone=29 ellps=WGS84")) # 29 = IRE or 30 = UK
#' add UTM locations to data frame
mydata_HI$x <- attr(utmcoord, "coords")[, 1]
mydata_HI$y <- attr(utmcoord, "coords")[, 2]

#' Prep data for HMM
mydata_HI <- mydata_HI %>% dplyr::select(ID, x, y, chloro)
mydata_HI <- data.frame(mydata_HI)
mydata_HI <-
  momentuHMM::prepData(
    data = mydata_HI,
    type = "UTM",
    coordNames = c("x", "y"),
    covNames = "chloro"
  )

#' plot the data using UTM
plot(mydata_HI$x, mydata_HI$y)

#' take a look at the histogram of step length to inform initial values
hist(mydata_HI$step)

# Indices of steps of length zero, used for zero mass value
whichzero <- which(mydata_HI$step == 0)
# Proportion of steps of length zero in the data set
length(mydata_HI) / nrow(mydata_HI)

#' take a look at the turning angle
hist(mydata_HI$angle, breaks = seq(-pi, pi, length = 15))

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
      7.433747e-05,
      7.433747e-05,
      7.433747e-05
    ),
    angle = c(40, 10, 10)
  )
# fit model
m1 <-
  momentuHMM::fitHMM(
    data = mydata_HI,
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
      7.433747e-05,
      7.433747e-05,
      7.433747e-05
    ),
    angle = c(40, 10, 10),
    formula = formula
  )
m2 <-
  momentuHMM::fitHMM(
    data = mydata_HI,
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
#' 1  1242879  m2  # with chloro
#' 2  1242985  m1  # without chloro
#' 
#' if HI_FB48028_086_0/0 is omitted: 
#' Model     AIC
#' 1    m2 1236484 # with chloro
#' 2    m1 1236555 # without chloro

#' save models
saveRDS(m1, "hmm_models/HI_no_chloro.rds")
saveRDS(m2, "hmm_models/HI_with_chloro.rds")

#' extract states and plot them
states <- viterbi(m2)
m2$data$state <- states
ggplot(data = m2$data, aes(x,y, color = as.factor(state))) + geom_point()

#' export the model details including states as a csv
m2$data <- rename(m2$data, ID_burst = ID)
m2$data$ID <- sub("_[^_]+$", "", m2$data$ID_burst) 
levels(as.factor(m2$data$ID))
head(m2$data)
m2$data$datetime <- time
hi_data_states <- m2$data %>% select(x,y,datetime,chloro,state,ID,ID_burst)
head(hi_data_states)
tail(hi_data_states)
write.csv(hi_data_states,file="results/HI_5min_states.csv",row.names = F)

#' stick on long and lat
hi_data_states <- read_csv("results/HI_5min_states.csv")
hi_data_states$long <- lon_lat$long
hi_data_states$lat <- lon_lat$lat

#' check the plots to see if the coords match
filter(hi_data_states, ID_burst == "HI_FB48036_103_27/0") %>% ggplot(., aes(x, y)) + geom_point()
filter(hi_data_states, ID_burst == "HI_FB48036_103_27/0") %>% ggplot(., aes(long, lat)) + geom_point()

#' export it again
head(hi_data_states)
write.csv(hi_data_states, file = "results/HI_5min_states.csv", row.names = F)

#' load in model if needed
hmm <- readRDS("hmm_models/HI_with_chloro.rds")
hmm
plot(hmm,plotCI = TRUE)

#' Plot stationary state probabilities
plotStationary(hmm, plotCI = TRUE, legend.pos = "center")

#' check the number of trips we're left with
hmm_data <- read_csv("results/HI_5min_states.csv")
hmm_data$year <- format(hmm_data$datetime, "%Y")
hmm_data$newID <- paste(hmm_data$ID, hmm_data$year, sep = "_")
length(levels(as.factor(hmm_data$newID))) # 65
table(hmm_data$newID)
#' FB48028_086 was tracked in 2014 and 2015
levels <- levels(as.factor(hmm_data$newID))
write.csv(levels, "HI_levels_HMM.csv", row.names = F)
#' check sample size excluding resting
filter(hmm_data,state<3)


