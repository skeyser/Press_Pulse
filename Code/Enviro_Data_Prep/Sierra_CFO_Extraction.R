## -------------------------------------------------------------
##
## Script name: cFO Extraction
##
## Script purpose: Extracting data on forest characteristics from
## ARU locations
##
## Author: Spencer R Keyser
##
## Date Created: 2025-09-24
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
library(sf)
library(terra)
library(exactextractr)
## -------------------------------------------------------------

## List the files for CFO
cfo.list <- list.files("C:/Users/srk252/Documents/GIS_Data/CFO/", full.names = T, pattern = "cfo")

cfo.r <- lapply(cfo.list, rast)
cfo.r <- do.call(c, cfo.r)
names(cfo.r) <- gsub("CFO-California-|-Summer2020-00010m", "", names(cfo.r))

plot(cfo.r)

## Load CFO data## Load CFrasterizeWinO data
# cancov <- rast("C:/Users/srk252/Documents/GIS_Data/CFO/CFO-California-CanopyCover-Summer2020-00010m.tif")
# plot(cancov)

## Testing the validity of supplying the same arus for all 4 years
locs <- st_read(here("Data/Spatial_Data/ARU_Locs_2021_2024.shp"))
locs <- locs |> st_transform(crs = crs(cfo.r))

## Generate the buffer
locs <- locs |> st_buffer(dist = units::set_units(120, "m"))

## Fast extract
cfo <- exact_extract(cfo.r, locs, fun = "mean")

locs <- cbind(locs, cfo)

colnames(locs) <- gsub("mean.", "", colnames(locs))

## Write to file
locs <- locs |> st_drop_geometry()
data.table::fwrite(locs, file = here("Data/CFO_ARU_21_24.csv"))

## Correlations
cor(locs |> select(CanopyBaseHeight:SurfaceFuels))
