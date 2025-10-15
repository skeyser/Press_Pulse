## -------------------------------------------------------------
##
## Script name: Geographic Variable extraction
##
## Script purpose: Extract geographic variables from 21-24 ARUS
##
## Author: Spencer R Keyser
##
## Date Created: 2025-10-09
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

## -------------------------------------------------------------

## Data loading
topo.paths <- list.files("C:/Users/srk252/Documents/GIS_Data/Topo/", full.names = T)
topo.r <- lapply(topo.paths, rast)
topo.r <- do.call(c, topo.r)

plot(topo.r)

## Load in the ARU locations
locs <- st_read(here("Data/Spatial_Data/ARU_Locs_2021_2024.shp"))
glimpse(locs)
head(unique(sort(locs$Cll_Unt)))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Buffering and extraction
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## generate the buffer - starting with 120 but can set variable buffers if needed
locs.buff <- locs |> 
  st_buffer(dist = units::set_units(120, "m")) |> 
  st_transform(crs = crs(topo.r))

## terra extraction
topo.ex <- exact_extract(topo.r, locs.buff, fun = "mean")

## Bind with the data
locs <- cbind(locs, topo.ex)

## Save the extraction
clncol <- function(x) gsub("mean.", "", x)
locs <- locs |> 
  st_drop_geometry() |> 
  rename_all(clncol)
  

## Write the file
write.csv(locs, file = here("Data/Topo_Data_ARU_21_24.csv"))
