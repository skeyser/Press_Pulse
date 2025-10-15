## -------------------------------------------------------------
##
## Script name: ARU Yearly Locations
##
## Script purpose: Generate a shapefile for all the ARU locations
## over all years for extracting variables
##
## Author: Spencer R Keyser
##
## Date Created: 2025-09-04
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
library(purrr)
library(ggplot2)
library(here)
library(CAbioacoustics)
library(CAbioacousticsextras)
library(sf)
## -------------------------------------------------------------

## Bring in the external functions
source(here("./Code/Enviro_Data_Prep/CBI_Extract_Funs.R"))

## Query all of the ARU locations
aru_locs <- cabio_loc_query(years = c(2021,
                                      2022,
                                      2023,
                                      2024)) |> 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2])

## Write the shapefile for ingesting into GEE
st_write(aru_locs, here("./Data/Spatial_Data/ARU_Locs_2021_2024.shp"), append = FALSE)

