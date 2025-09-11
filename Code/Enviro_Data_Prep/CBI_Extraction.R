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
                          landscape_metrics = T
){
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CAbioacoustics Spatial
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(locs_from_cabio){
    aru_locs <- cabio_loc_query(years = survey_years)
  } else if (!is.null(custom_locs)) {
    aru_locs <- custom_locs
  } else {
    stop("No locations provided. Either set locs_from_cabio = TRUE or provide locations for custom_locs.")
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Create SF objects from the ARU coordinates
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Create SF objects from points
  if(!any(class(aru_locs) %in% "sf")){
    
    aru_locs <- st_as_sf(aru_locs, 
                         coords = c(x_col, 
                                    y_col), 
                         crs = .crs)
  }
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: CBI Categorical Extraction and prep
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Find the data for extraction
  #if(env_prod == "Fire_CBI"){
  
  ## Function bundle
  
  if(file.exists("C:/Users/srk252/Documents/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack.tif")){
    
    message("CBI rasters are fixed and exist in directory. Loading the fixed stack.")
    cbi_stack <- rast("C:/Users/srk252/Documents/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack.tif")
  
    } else {
    
    ## Find files
    cbi_path <- "C:/Users/srk252/Documents/data_for_spencer/cbi_sierra_cat_rasters/"
    cbi_files <- list.files(cbi_path, full.names = T, pattern = "(cbi_cat_)(\\d{4})(*.tif$)")
    
    ## Call function to stack if need be
    cbi_stack <- cbi_stack_fun(cbi_files)
    
    ## Add in zeros for the raster for downstream processing
    temp <- rast(nrows = nrow(cbi_stack[[1]]),
                 ncols = ncol(cbi_stack[[1]]),
                 xmin = xmin(cbi_stack[[1]]),
                 xmax = xmax(cbi_stack[[1]]),
                 ymin = ymin(cbi_stack[[1]]),
                 ymax = ymax(cbi_stack[[1]]),
                 crs = crs(cbi_stack[[1]]),
                 resolution = res(cbi_stack[[1]]),
                 vals = 0,
                 names = "template"
    )
    
    ## Reclassify the NAs to zero to fill in the rasters
    cbi_stack <- classify(cbi_stack, cbind(NA, 0))
    if(nlyr(cbi_stack) == length(cbi_files)){
      names(cbi_stack) <- str_extract(cbi_files, "\\d{4}")
    }
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Fire year interval stacking
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(!is.null(intervals)){
  cbi_int <- int_ras_fun(ras_stack = cbi_stack,
                         intervals = intervals,
                         sum_int = "max",
                         locations = aru_locs)
  } else {
    message("No intervals specified. Skipping interval binning.")
  }
  
  ## Project the points
  aru_locs <- st_transform(aru_locs,
                           crs = crs(cbi_stack))
  
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
  ## *************************************************************
  if("fire_severity" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals set. Fire severity calculated from summarized fire years.")
      cbi_sev <- cbi_int
    } else {cbi_sev <- cbi_stack}
    if(!is.null(buff_size)){
      aru_buffer <- st_buffer(aru_locs,
                              dist = buff_size)
      
      ## Check the number of unique years
      if(length(unique(aru_locs$survey_year)) == length(cbi_sev)){
        print("Matching fire intervals to the unique years in the data. Proceed.")
      } else { stop("Mismatching number of years for extraction.") }
      
      ## If the years are more than 1
      if(length(unique(aru_locs$survey_year)) > 1){
        fire_sev <- vector(mode = "list", length = length(unique(aru_locs$survey_year)))
        yrs_tmp <- unique(aru_locs$survey_year)
        
        for(i in 1:length(unique(aru_buffer$survey_year))){
          fire_sev[[i]] <- exactextractr::exact_extract(cbi_sev[[i]],
                                                        aru_buffer |> filter(survey_year == yrs_tmp[i]),
                                                        fun = c("mean", "stdev"),
                                                        append_cols = id_col)
        }
      }
    } else {
      
      fire_sev <- terra::extract(cbi_sev,
                                 aru_locs,
                                 bind = T) |> st_as_sf()
      
    }
    colnames(fire_sev)[colnames(fire_sev) != id_col] <- paste0("Fire_Sev_", gsub("[[:punct:]]", "_", colnames(fire_sev)[colnames(fire_sev) != id_col]))
    
  } else {
    fire_sev <- NULL
  }
  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Calculate fire variables
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ras_years <- as.numeric(names(cbi_stack))
  
  if("time_since_fire" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Time since fire doesn't accept intervals...output will be generated from single years.")
    }
    
    ## Time to most recent fire
    time_since_fire <- terra::app(cbi_stack, 
                                  fun = function(x) time_to_most_recent_fire(cell_values = x,
                                                                             years = ras_years))
  }
  
  if("fire_freq" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Fire frequency doesn't accept intervals...output will be generated from single years.")
    } 
    ## Fire frequency (num fires/total record length)
    fire_freq <- terra::app(cbi_stack, 
                            fun = function(x) fire_freq_calc(cell_values = x,
                                                             years = ras_years))
  }
  
  if("fire_ret_int" %in% fire_prod){
    if(!is.null(intervals)){
      message("Intervals are set. Fire return interval doesn't accept intervals.,,output will be generated from single years.")
    } 
    ## Fire return interval (mean of time between successive fires)
    fire_return_int <- terra::app(cbi_stack, 
                                  fun = function(x) fire_return_int(cell_values = x, 
                                                                    years = ras_years))
  }
  
  ## Report the status of s2 geometry for convenience
  if(!is.null(buff_size)){
    if(sf_use_s2() == T){
      message("s2 geometry enabled. Buff_size interpreted as meters.")
    } else {
      message("s2 disabled, if locations are geodetic (lat/lon) units interpretted as degrees.")
    }
  }
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  ## Subsection: Buffer Extraction from Custom Fire Variables
  ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(any(fire_prod %in% c("time_since_fire", "fire_freq", "fire_ret_int"))){
    variable_name <- fire_prod[!fire_prod == "fire_severity"]
    fire_buff_out <- vector(mode = "list", length = length(variable_name))
    
    for(var in 1:length(variable_name)){
      
      var.tmp <- variable_name[var]
      
      if(var.tmp == "time_since_fire"){
        fire_buff_extract <-
          buff_size %>%
          # create column names
          str_c(var.tmp, ., sep = '_') |>
          set_names() |>
          map_dfr(
            \(x)
            exactextractr::exact_extract(
              time_since_fire,
              # buffer points
              aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
              fun = 'mean'
            )
          ) |> 
          bind_cols(st_drop_geometry(aru_locs[,id_col]))
        
        fire_buff_out[[var]] <- fire_buff_extract
        
      }
      
      if(var.tmp == "fire_freq"){
        fire_buff_extract <-
          buff_size %>%
          # create column names
          str_c(var.tmp, ., sep = '_') |>
          set_names() |>
          map_dfr(
            \(x)
            exactextractr::exact_extract(
              fire_freq,
              # buffer points
              aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
              fun = 'mean'
            )
          ) |> 
          bind_cols(st_drop_geometry(aru_locs[,id_col]))
        
        fire_buff_out[[var]] <- fire_buff_extract
      }
      
      if(var.tmp == "fire_ret_int"){
        fire_buff_extract <-
          buff_size %>%
          # create column names
          str_c(var.tmp, ., sep = '_') |>
          set_names() |>
          map_dfr(
            \(x)
            exactextractr::exact_extract(
              fire_return_int,
              # buffer points
              aru_locs |> st_buffer(as.numeric(str_extract(x, "\\d"))),
              fun = 'mean'
            )
          ) |> 
          bind_cols(st_drop_geometry(aru_locs[,id_col]))
        
        fire_buff_out[[var]] <- fire_buff_extract
      }
      
    }
    fire_buff_merge <- Reduce(function(x, y) merge(x, y, by = id_col), fire_buff_out) 
  } else {
    fire_buff_merge <- NULL
  }
  
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
                            survey_years = c(2021, 2022, 2023, 2024),
                            intervals = c("1-5", "6-10"), #only interested in recent fire
                            id_col = "deployment_name",
                            buff_size = buffSize,
                            landscape_metrics = T
)

fire_sev21 <- fire_sev21$FireSeverity


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





