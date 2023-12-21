### Applied Geodata Science ###
### Phenology Model ###


# libraries ---------------------------

#install.packages("daymetr")
library(tidyverse) # careful, masks terra::extract
library(terra) # process raster data
library(sf) # prcess shape file data
library(MODISTools) # API, access MODIS data
library(cowplot) # combining plots
library(patchwork) # also combining plots
library(signal) # for smoothing
library(leaflet)
# phenocam network API
# (download time series of vegetation greenness and derived phenology metrics)
library(phenocamr)
library(GenSA)
library(daymetr)


######################
# get data -------------------------

# download greenness time series,
# calculate phenology (phenophases),
# amend with DAYMET data
phenocamr::download_phenocam(
  site = "harvard$",
  veg_type = "DB",
  roi_id = "1000",
  daymet = TRUE,
  phenophase = TRUE,
  trim = 2022,
  out_dir = tempdir())

harvard_phenocam_data <- readr::read_csv(
  file.path(tempdir(), "harvard_DB_1000_3day.csv"),
  comment = "#")

# reading in harvard phenology only retaining
# spring (rising) phenology for the GCC 90th
# percentile time series (the default)
harvard_phenology <- readr::read_csv(
    file.path(tempdir(), "harvard_DB_1000_3day_transition_dates.csv"),
    comment = "#") |>
  dplyr::filter(direction == "rising",
                gcc_value == "gcc_90")


# plot phenology data -----------

ggplot(harvard_phenocam_data) +
  geom_line(aes(as.Date(date), smooth_gcc_90),
            colour = "grey25") +
  geom_point(data = harvard_phenology,
             aes(as.Date(transition_25), threshold_25)) +
  labs(x = "", y = "GCC") +
  theme_bw() +
  theme(legend.position = "none")
######################
# GDD <- get growing degree days ----------------------

# GDD: cumulative sum of temperatures above a specified threshold (usually T0=5)
# More advanced models exist. Such models may include not only temperature but also
# radiation, precipitation and or temporal lags or so-called chilling requirements,
# frost during the preceding winter months. For an overview of these more advanced models
# I refer to Basler (2016) and Hufkens et al. (2018).

# return mean daily temperature and formal dates (for plotting)
harvard_temp <- harvard_phenocam_data |>
  group_by(year) |>
  dplyr::mutate(tmean = (tmax..deg.c. + tmin..deg.c.)/2) |>
  dplyr::mutate(date = as.Date(date),
                gdd = cumsum(ifelse(tmean >= 5, tmean - 5, 0))) |>
  dplyr::select(date,
                year,
                tmean,
                gdd) |>
  ungroup()

# convert the harvard phenology data and only retain required data
harvard_phenology <- harvard_phenology |>
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) |>
  select(
    year,
    doy,
    transition_25,
    threshold_25
  )

# IMPORTANT: (see plot)
# In 2010, Spring leaf development started at day of year (DOY) 114.
# GDD for that day was 130.44°C.

# GOAL:
# find a GDD that returns accurate leaf-out-day predictions for multiple years.

# plot GDD ------------------

# daily temperature (red = above threshhold, blue below)
# cumulative GDD

# grab only the 2010 value of spring phenology
harvard_phenology_2010 <- harvard_phenology |>
  dplyr::filter(year == 2010)

harvard_gdd_value <- harvard_temp |>
  dplyr::filter(date == harvard_phenology_2010$transition_25)

p <- ggplot(harvard_temp) +
  geom_line(aes(date, tmean)) +
  geom_point(aes(date, tmean,
                 colour = tmean > 5,
                 group = 1)) +
  geom_vline(data = harvard_phenology_2010,
             aes(xintercept = as.Date(transition_25))) +
  scale_colour_discrete(type = c("blue","red")) +
  labs(x = "", y = "Temperature (deg. C)") +
  xlim(c(as.Date("2010-01-01"),
         as.Date("2010-06-30"))) +
  theme_bw() +
  theme(legend.position = "none")

p2 <- ggplot(harvard_temp) +
  geom_line(aes(date, gdd)) +
  geom_point(aes(date, gdd,
                 colour = tmean > 5,
                 group = 1)) +
  scale_colour_discrete(type = c("blue","red")) +
  geom_vline(data = harvard_phenology_2010,
             aes(xintercept = as.Date(transition_25))) +
  geom_hline(data = harvard_gdd_value,
             aes(yintercept = gdd),
             lty = 2) +
  labs(x = "", y = "GDD (deg. C)") +
  xlim(c(as.Date("2010-01-01"),
         as.Date("2010-06-30"))) +
  ylim(c(0, 1000)) +
  theme_bw()  +
  theme(legend.position = "none")

# compositing with library patchwork
p + p2  +
  plot_layout(ncol = 1) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")")



######################
# GDD model optimization -------------

# function is fed with
# 1) temp timeseries
# 2) vector with 2 parameters:
#   a. temperature threshhold (eg T0=5)
#   b. gdd at which leaf out is predicted


gdd_model <- function(temp, par) {
  # split out parameters from a simple vector of parameter values
  temp_threshold <- par[1]
  gdd_crit <- par[2]

  # accumulate growing degree days for temperature data
  gdd <- cumsum(ifelse(temp > temp_threshold, temp - temp_threshold, 0))

  # figure out when the number of growing degree days exceeds min value
  # required for leaf development, only return the first value
  doy <- unlist(which(gdd >= gdd_crit)[1])

  return(doy)
}


# TEST model:
# with T0 = 5 and GDD = 130.44 (see above),
# DOY 114 should be predicted for 2010

# confirm that the model function returns expected results (i.e. DOY 114)
# (we filter out the year 2010, but removing the filter would run the
# model for all years!)
prediction <- harvard_temp |>
  dplyr::filter(year == 2010) |>
  group_by(year) |>
  summarize(pred = gdd_model(temp = tmean,
                             par = c(5, 130.44)))

print(prediction)


######################
# GDD model calibration (parameter tuning) -------------

# CAREFUL!! this is a slow implementation.don't copy for large datasets!

# evaluation function rmse_gdd
source("./R/rmse_gdd.R")

# starting model parameters
par = c(0, 130)

# limits to the parameter space
lower <- c(-10,0)
upper <- c(45,500)

# data needs to be provided in a consistent
# single data file, a nested data structure
# will therefore accept non standard data formats
data <- list(drivers = harvard_temp,
             validation = harvard_phenology)

# optimize the model parameters (usually with caret::train!) --> 2.81 °C treshhold, 228 GDD
optim_par = GenSA::GenSA(
  par = par,
  fn = rmse_gdd, # own function goes here
  lower = lower,
  upper = upper,
  control = list(max.call = 4000),
  data = data
  )$par


# run the model for all years to get the phenology predictions
predictions <- harvard_temp |>
  group_by(year) |>
  summarize(
    prediction = gdd_model(
      temp = tmean,
      par = optim_par
    )
  )


# plot prediction for all years with optimized parameters --------------

# join predicted with observed data
validation <- left_join(predictions, harvard_phenology)

ggplot(validation) +
  geom_smooth(
    aes(
      doy,
      prediction
    ),
    colour = "grey25",
    method = "lm"
  ) +
  geom_point(
    aes(
      doy,
      prediction
    )
  ) +
  geom_abline(
    intercept=0,
    slope=1,
    linetype="dotted"
  ) +
  labs(
    x = "Observed leaf-out date (DOY)",
    y = "Predicted leaf-out date (DOY)"
  ) +
  theme_bw()  +
  theme(
    legend.position = "none"
  )

######################
# Spatial scaling ---------------

# Download daily temp data for bigger area
daymetr::download_daymet_tiles(
  tiles = 11935,
  start = 2012,
  end = 2012,
  param = c("tmin","tmax"),
  path = paste0(here::here(), "/data_raw/"),
  silent = TRUE
)

# calculate the daily mean values
# ERROR
r <- daymetr::daymet_grid_tmean(
  path = paste0(here::here(), "/data_raw/"),
  product = 11935,
  year = 2012,
  internal = TRUE
)

# reproject to lat lon
r <- terra::project(r, "+init=epsg:4326")

# can't read r --> TINOS FILE
r1 <- readRDS("./data/first_r_in_pheno.RDS") # needs to be projected (line above)
r2 <- readRDS("./data/r.RDS") # already projected

# subset to first 180 days
# (smaller memory and we're only interested in start of growing season)
ma_nh_temp <- terra::subset(r2,1:180)

predicted_phenology <- terra::app(
  ma_nh_temp,
  fun = gdd_model,
  par = optim_par
)


# plot results -----------------------

library(leaflet)

# set te colour scale manually
pal <- colorNumeric(
  "magma",
  values(predicted_phenology),
  na.color = "transparent"
)

# build the leaflet map
# using ESRI tile servers
# and the loaded demo raster
leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") |>
  addProviderTiles(providers$Esri.WorldTopoMap, group = "World Topo") |>
  addRasterImage(
    predicted_phenology,
    colors = pal,
    opacity = 0.8,
    group = "Phenology model results"
  ) |>
  addLayersControl(
    baseGroups = c("World Imagery","World Topo"),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE),
    overlayGroups = c("Phenology model results")
  ) |>
  addLegend(
    pal = pal,
    values = values(predicted_phenology),
    title = "DOY")
