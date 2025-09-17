## -------------------------------------------------------------
##
## Script name: CBI Fire Extraction
##
## Script purpose: Extract and summarize fire variables for
## ARU locations across the Sierra Nevada
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-17
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
library(dplyr)
library(stringr)
library(ggplot2)
library(here)
library(terra)
library(sf)
library(purrr)
library(CAbioacoustics)
## -------------------------------------------------------------

## Load functions for CBI data processing
source(here("./Code/Enviro_Data_Prep/CBI_Extract_Funs.R"))

## -------------------------------------------------------------
##
## Begin Section: Extraction function body
##
## -------------------------------------------------------------

## *************************************************************
##
## Section Notes:
## Environmental extraction function
## Acceptable fire products
## fire_interval, most_recent_fire, fire_frequency, number of fires,
## and landscape metrics
##
## To-Do
## Move the main function to separate script so this script
## only involve running
## *************************************************************

aru_fire_prep <- function(fire_prod = NULL, # character vector of desired fire output
                          locs_from_cabio = T, #override flag for using CAbioacoustic locations and metadata
                          custom_locs = NULL, # Data.frame with coordinates for custom locations                          
                          survey_years = c(2021,
                                           2022
                                           #2023,
                                           #2024
                                           ), # Survey year is only applicable for locs_from_cabio = TRUE
                          id_col = "deployment_name",
                          year_col = NULL,
                          x_col = "Long", # chr for x coordinate col name
                          y_col = "Lat", # chr for y coordinate col name
                          .crs = 4326, # default WGS84 for coordinates
                          des_out, # desired output
                          spat_ex, # chr for type (point, buff, hex)
                          buff_size = 120, # vector of buffer sizes
                          intervals = c("1-5", "6-10", "11-35"),
                          landscape_metrics = T,
                          allow_raster_time_gaps = FALSE
){
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CAbioacoustics Spatial
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  aru_locs <- prep_locs(locs_from_cabio = locs_from_cabio,
                        custom_locs = custom_locs,
                        survey_years = survey_years,
                        x_col = x_col,
                        y_col = y_col,
                        .crs = .crs)
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CBI Categorical Extraction and prep
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cbi_stack <- cbi_clean(cbi_path = "C:/Users/srk252/Documents/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack.tif")
  
  ## Project the points
  aru_locs <- st_transform(aru_locs,
                           crs = crs(cbi_stack))
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Intial input validation
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  validate_inputs(fire_prod = fire_prod, 
                  locs_from_cabio = locs_from_cabio, 
                  custom_locs = custom_locs, 
                  survey_years = survey_years,
                  ras_stack = cbi_stack,
                  intervals = intervals,
                  allow_gaps = allow_raster_time_gaps)
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Fire year interval stacking
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cbi_int <- process_intervals(cbi_stack = cbi_stack,
                               intervals = intervals,
                               locations = aru_locs)
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Fire Severity
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## *************************************************************
  ##
  ## Section Notes: Should be modified to accommodate more than 1 year
  ## with a shifting baseline of the year for fire. It's easy enough to 
  ## do this in a small loop within this function but the interval 
  ## function returns the actual years which makes a merge more difficult.
  ##
  ## Still a WIP...but the yearly extractions are fine
  ## Added in a identifier for the survey year for book keeping
  ## To-Do
  ## 1. Edit functionality for the single point extraction
  ## 2. Move to separate function
  ## 3. Modify the full script for better handling of fire products
  ## *************************************************************
  if("fire_severity" %in% fire_prod){
    message("fire_severity specified...extracting fire severity data.")
    fire_sev <- extract_fire_sev(cbi_sev = cbi_int, aru_locs = aru_locs,
                               buff_size = buff_size, id_col = id_col)
  } else {
    fire_sev <- NULL
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Calculate fire variables if specified
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  custom_fire_mets <- create_fire_metrics(cbi_stack = cbi_stack, fire_prod = fire_prod, intervals = intervals)
  
  ## Function needs to be returned but will move when code below is fixed.
  if(length(custom_fire_mets) == 0){
    fire_buff_merge <- NULL
  }
  
  ## *************************************************************
  ##
  ## Section Notes: WIP
  ## To-Do
  ## 1. Need to build a better function for toggling buffer vs non-buffer
  ## 2. Need a multi-year shift for ARU surveys collected at different years
  ## 3. Multi-buffer size extraction give the slow estimation speed
  ## 4. Optional region clips? Maybe the user can add this before with the CBI
  ## stack itself?
  ##
  ## *************************************************************
  # ## Report the status of s2 geometry for convenience
  # if(!is.null(buff_size)){
  #   if(sf_use_s2() == T){
  #     message("s2 geometry enabled. Buff_size interpreted as meters.")
  #   } else {
  #     message("s2 disabled, if locations are geodetic (lat/lon) units interpretted as degrees.")
  #   }
  # }
  # 
  #
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Buffer Extraction from Custom Fire Variables
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # if(any(fire_prod %in% c("time_since_fire", "fire_freq", "fire_ret_int"))){
  #   variable_name <- fire_prod[!fire_prod == "fire_severity"]
  #   fire_buff_out <- vector(mode = "list", length = length(variable_name))
  #   
  #   for(var in 1:length(variable_name)){
  #     
  #     var.tmp <- variable_name[var]
  #     
  #     if(var.tmp == "time_since_fire"){
  #       fire_buff_extract <-
  #         buff_size %>%
  #         # create column names
  #         str_c(var.tmp, ., sep = '_') |>
  #         set_names() |>
  #         map_dfr(
  #           \(x)
  #           exactextractr::exact_extract(
  #             time_since_fire,
  #             # buffer points
  #             aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
  #             fun = 'mean'
  #           )
  #         ) |> 
  #         bind_cols(st_drop_geometry(aru_locs[,id_col]))
  #       
  #       fire_buff_out[[var]] <- fire_buff_extract
  #       
  #     }
  #     
  #     if(var.tmp == "fire_freq"){
  #       fire_buff_extract <-
  #         buff_size %>%
  #         # create column names
  #         str_c(var.tmp, ., sep = '_') |>
  #         set_names() |>
  #         map_dfr(
  #           \(x)
  #           exactextractr::exact_extract(
  #             fire_freq,
  #             # buffer points
  #             aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
  #             fun = 'mean'
  #           )
  #         ) |> 
  #         bind_cols(st_drop_geometry(aru_locs[,id_col]))
  #       
  #       fire_buff_out[[var]] <- fire_buff_extract
  #     }
  #     
  #     if(var.tmp == "fire_ret_int"){
  #       fire_buff_extract <-
  #         buff_size %>%
  #         # create column names
  #         str_c(var.tmp, ., sep = '_') |>
  #         set_names() |>
  #         map_dfr(
  #           \(x)
  #           exactextractr::exact_extract(
  #             fire_return_int,
  #             # buffer points
  #             aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
  #             fun = 'mean'
  #           )
  #         ) |> 
  #         bind_cols(st_drop_geometry(aru_locs[,id_col]))
  #       
  #       fire_buff_out[[var]] <- fire_buff_extract
  #     }
  #     
  #   }
  #   fire_buff_merge <- Reduce(function(x, y) merge(x, y, by = id_col), fire_buff_out) 
  # } else {
  #   fire_buff_merge <- NULL
  # }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Landscape Metrics
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Function should estimate specific landscape metrics for the
  ## interval-based raster data for CBI
  if(landscape_metrics){
    fire_lscp_out <- fire_lscp_fun(ras_int = cbi_int,
                                   ras_stack = cbi_stack,
                                   locs = aru_locs,
                                   years = survey_years,
                                   intervals = intervals,
                                   buff_size = buff_size,
                                   id_col = id_col,
                                   metrics = c("lsm_c_pland"))
  } else {
    fire_lscp_out <- NULL
  }

  return(list(FireSeverity = fire_sev,
              FireMetrics = fire_buff_merge,
              FireLscp = fire_lscp_out))

} ## function closure


buffSize <- 120

fire_sev21 <- aru_fire_prep(fire_prod = c("fire_severity"),
                            locs_from_cabio = TRUE,
                            survey_years = 2021:2024,
                            intervals = c("1-5", "6-10"), #only interested in recent fire
                            id_col = "deployment_name",
                            buff_size = buffSize,
                            landscape_metrics = TRUE
)

str(fire_sev21)
fire_sev21 <- fire_sev21$FireSeverity

## Write the robject
saveRDS(fire_sev21, file = here("./Data/Fire_ARU_21_24.RDS"))


## Histogram for the fire severity
nrow(fire_sev21$FireSeverity)
hist(fire_sev21$FireSeverity$Fire_Sev_mean_2016_2020)
hist(fire_sev21$FireSeverity$Fire_Sev_mean_2011_2015)
hist(fire_sev21$FireSeverity$Fire_Sev_mean_1986_2010)

## Bring in the dataframe from the ARU stuff
aru_meta <- read.csv(here("./Data/ARU_120m_New.csv")) |> 
  mutate(deployment_name = paste(group_id, visit_id, cell_id, unit_numbe, sep = "_"))

aru_meta <- aru_meta |> 
  filter(deployment_name %in% fire_sev21$FireSeverity$deployment_name) |> 
  left_join(fire_sev21$FireSeverity)

write.csv(fire_sev21, file = here("./Data/FireSeverity2021_MeanStDev_AllARUs.csv"))





