## -------------------------------------------------------------
##
## Script name: Anomaly Data Formatting
##
## Script purpose: Formatting anomaly data from GEE
##
## Author: Spencer R Keyser
##
## Date Created: 2025-08-06
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
## -------------------------------------------------------------
##
## Notes:
##
##
## -------------------------------------------------------------

## Defaults
options(scipen = 10, digits = 10)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(here)

## -------------------------------------------------------------

## Load in the anomalies data
anom <- readr::read_csv("C:/Users/srk252/Documents/GIS_Data/Sierra_Climate_Anomalies/Temperature_Anomalies_by_ARU.csv")
glimpse(anom)

anom <- anom |> 
  select(-`system:index`, -.geo)

## Look at the overall trend in anomaly
anom |> 
  group_by(anomaly_year) |> 
  summarise(anom_mn = mean(anomaly_value)) |> 
  ggplot(aes(x = anomaly_year, y = anom_mn)) +
  geom_point() + 
  geom_line() + 
  theme_bw() + 
  xlab("Year") + 
  ylab("Mean T. Anomaly from 1980-2010 (Deg. C)")

anom |> 
  ggplot(aes(x = anomaly_year, y = anomaly_value, color = Y)) +
  geom_point(position = "jitter") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_viridis_c(option = "C", direction = -1) +
  scale_x_continuous(breaks = c(2011, 2013, 2015, 2017, 2019, 2021, 2023)) +
  xlab("Year") + 
  ylab("Max T. Anomaly from 1980-2010 (Deg. C)") + 
  labs(color = "Latitude") + 
  theme_bw()

test <- anom |> 
  group_by(X, Y) |> 
  summarise(anom_mn_site = mean(anomaly_value))

plot(test$Y, test$anom_mn_site)

## Trends in climate

