## -------------------------------------------------------------
##
## Script name: ARU Species Detection Pre-processing
##
## Script purpose: Format and process species detection histories
## from ARUs and link together with spatial and environmental metadata.
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-09
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
options(scipen = 6, digits = 4)

## -------------------------------------------------------------

## Package Loading
# Package Loading
library(here)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(CAbioacoustics)
library(sf)
library(mapview)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Thresholds file exploration
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes:
## Variable types
## 1. Species {character} = Species Common Name
## 2. keep {numeric} = Binary indicator variable
## 3. BirdNET_code {numeric} = Species' code names from BirdNET
## 4. intercept.c {numeric}
## 5. beta.c {numeric}
## 6. intercept.r {numeric}
## 7. beta.r {numeric}
## 8. cutoff85.r {numeric}
## 9. cutoff90.r {numeric}
## 10. cutoff95.r {numeric}
## 11. cutoff975.r {numeric}
## 12. cutoff99.r {numeric}
## 13. cutoff90.r_conf {numeric}
## 14. cutoff95.r_conf {numeric}
## 15. cutoff99.r_conf {numeric}
##
## *************************************************************

## -------------------------------------------------------------
##
## Begin Section: Processing species detections
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes: Source ac_det_filter() from Sierra_functions.R
## 
## To-do:
## 1. Streamline the multi-year species data processing DONE
## 2. Unify with spatial coordinates from CAbioacoustics DONE
## 3. Write yearly files for community data (summarized; naive biodiv. models) WIP
## 4. Write yearly files for community data (date-detection; occ models) DONE
## 5. Add flexibility for the setting wd WIP
##    a. Make sure relative and complete paths work
##    b. Allow user to provide both quotes and non-quotes to path
## 6. Add flexibility to the date string so that users can provide one date for all years
## 7. Add flexibility to the species function so the user can select a priori species to remove
##
## *************************************************************


## Source the ac_det_filter function
source(here("./Code/Acoustic_Data_Prep/Sierra_functions.R"))

## threshold file
thresh <- readr::read_csv("./Data/Thresholds_2021_20230309_BestThreshold_975min.csv")

## Run the function
aru_det_file_gen(det_dir = "C:/Users/srk252/Documents/Rprojs/Press_Pulse/Data/Detections_By_Species/",
                 det_years = c("2021", "2022", "2023", "2024"),
                 seas_format = F,
                 seas_outdir = "C:/Users/srk252/Documents/Rprojs/Press_Pulse/Data/Seasonal_Summaries/",
                 occ_format = T,
                 occ_outdir = "C:/Users/srk252/Documents/Rprojs/Press_Pulse/Data/Occ_Data/Thresh_By_Species/",
                 eff_file = T,
                 coord_link = T,
                 d_thresh = thresh,
                 thresh_scale = "Conf",
                 thresh_transform = TRUE,
                 thresh_trans_dir = "conf2logit",
                 thresh_cut = "99",
                 species_thresh_cut = "BestThresh",
                 time_format = "ymd",
                 no_dets = 2,
                 binary = T,
                 date_range = c("2021-06-01", "2021-06-30",
                                "2022-06-01", "2022-06-30",
                                "2023-06-01", "2023-06-30",
                                "2024-06-01", "2024-06-30"),
                 eff_site_name = "Cell_U",
                 eff_filter = 10,
                 verbose = F)

