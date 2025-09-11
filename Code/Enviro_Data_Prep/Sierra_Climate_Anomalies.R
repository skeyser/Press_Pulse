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
library(terra)
library(sf)
library(exactextractr)

## -------------------------------------------------------------

## Load in the anomalies data
## Waiting on the next batch of extracted data
## Need to think about buffer sizes for the species extractions
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Expanding this to workflow for extraction
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Order the .tifs for the anomalies by date
anom.path <- sort(list.files(here("Data/Spatial_Data/Climate_Anomaly_Rasters/"), pattern = "Tmax_Anomaly*", full.names = T))

## read in the tifs
anom.r <- lapply(anom.path, rast)
anom.r <- rast(anom.r)

## Give informative names
names(anom.r) <- anom.r@pntr@.xData$get_sourcenames()

## Load in the ARU locations
locs <- st_read(here("Data/Spatial_Data/ARU_Locs_2021_2024.shp"))
glimpse(locs)

## generate the buffer - starting with 120 but can set variable buffers if needed
locs.buff <- locs |> st_buffer(dist = units::set_units(120, "m"))

## terra extraction
anom.ex <- exact_extract(anom.r, locs.buff, fun = "mean")

## Bind with the data
locs <- cbind(locs, anom.ex)

## Locs in long format
locs_long <- locs |> 
  tidyr::pivot_longer(cols = contains("Tmax_Anomaly"), names_to = "Year", values_to = "Anomaly") |> 
  mutate(Year = as.numeric(stringr::str_extract(Year, "\\d+")))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: carpentry and plotting with full extracted data
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Check out the extractions
glimpse(locs)

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

## Original figure year x anom x lat
anom_plot_points <- locs_long |> 
  ggplot(aes(x = Year, y = Anomaly, color = Lat)) +
  geom_point(position = "jitter") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_viridis_c(option = "C", direction = -1) +
  scale_x_continuous(breaks = c(2011, 2013, 2015, 2017, 2019, 2021, 2023)) +
  xlab("Year") + 
  ylab("Max T. Anomaly from 1980-2010 (Deg. C)") + 
  labs(color = "Latitude") + 
  theme_bw()

## Some alpha coding
anom_plot_pointsmooth <- locs_long |> 
  ggplot(aes(x = Year, y = Anomaly, color = Lat)) +
  geom_point(alpha = 0.3, position = position_jitter(width = 0.2)) + 
  geom_smooth(aes(group = 1), color = "black", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_viridis_c(option = "D", direction = -1) +
  scale_x_continuous(breaks = seq(2011, 2023, by = 2)) +
  xlab("Year") + 
  ylab("Max T. Anomaly from 1980-2010 (Deg. C)") + 
  labs(color = "Latitude") + 
  theme_bw()

## Smoothing
anom_plot_pointsmooth <- locs_long |> 
  mutate(Lat_category = case_when(
           Lat >= 40 ~ "High",
           Lat >= 38 & Lat < 40 ~ "Mid",
           Lat < 38 ~ "Low"
         ),
         Lat_category = factor(Lat_category, levels = c("Low", "Mid", "High"))) |> 
  ggplot(aes(x = Year, y = Anomaly, color = Lat_category)) +
  geom_point(alpha = 0.3, position = position_jitter(width = 0.2)) + 
  geom_smooth(aes(group = Lat_category), method = "loess", se = F) +  # Changed this line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("#FF9E44", "#1E88E5", "#004D40")) +
  scale_x_continuous(breaks = seq(2011, 2023, by = 2)) +
  xlab("Year") + 
  ylab("Max T. Anomaly from 1980-2010 (Deg. C)") + 
  labs(color = "Latitude Category") + 
  theme_bw()
