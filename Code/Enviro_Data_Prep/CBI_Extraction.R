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

## Testing the validity of supplying the same arus for all 4 years
locs <- st_read(here("Data/Spatial_Data/ARU_Locs_2021_2025.shp"))

## Small shifts in ARU locations make this a weird task
## Could estimate the shift per unit and accommodate with
## buffering...?
unique_units <- unique(locs$Cll_Unt)
table(locs$srvy_yr)

loc_option <- 2

if(loc_option == 1){
## Option 1
## Location are represented as an average 
locs_fixed <- locs |> 
  select(Cell_Unit = Cll_Unt, deployment_name = dplymn_, Long, Lat) |> 
  distinct() |> 
  arrange(Cell_Unit) |> 
  st_drop_geometry() |> 
  tidyr::crossing(survey_year = 2021:2025) |> 
  group_by(Cell_Unit, deployment_name, survey_year) |> 
  summarise(Lat = mean(Lat),
            Long = mean(Long)) |> 
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

locs <- locs_fixed
print("Using the centroid value for all observations and filling.")

} else if(loc_option == 2){

## Option 2
## Locations from the observed points are variable (accommodate shifts)
## but the missing years are the centroid value

## Observed ARU locations
locs_obs <- locs |> 
  select(Cell_Unit = Cll_Unt,
         deployment_name = dplymn_,
         survey_year = srvy_yr,
         geometry)

## Centroids to act as filler for the missing survey years
locs_centroids <- locs_obs |> 
  group_by(Cell_Unit, deployment_name) |> 
  summarise(geometry = st_centroid(st_union(geometry)), .groups = "drop")

## Expanded year set
years <- tibble(survey_year = 2021:2025)

site_year_grid <- locs_centroids |> 
  tidyr::crossing(years)

## Join to bring in the new sites
locs_cent_fill <- site_year_grid |> 
  left_join(
    locs_obs |> rename(geometry_obs = geometry),
    by = c("Cell_Unit", "deployment_name", "survey_year")
  )

## Replace the centroid with observed
locs_mixed <- locs_cent_fill |> 
  mutate(
    geometry = if_else(
      !st_is_empty(geometry_obs),
      geometry_obs,
      geometry
    )
  ) |> 
  select(Cell_Unit, deployment_name, survey_year, geometry) |> 
  st_as_sf()

locs <- locs_mixed
print("Filling coordinated with centroid position to deal with jitter.")

}


## Checking some distance based measurements
rep_d <- locs |> 
  st_transform(., crs = 3310) |> 
  group_by(deployment_name) |> 
  filter(n() > 1) |> 
  summarise(
    max_dist = {
      dist_mat <- st_distance(geometry)
      max(dist_mat)
    },
    mean_dist = {
      dist_mat <- st_distance(geometry)
      mean(dist_mat)
    },
    min_dist = {
      dist_mat <- st_distance(geometry)
      min(dist_mat)
    },
    num_dep = {
      length(geometry)
    },
    .groups = "drop"
  )

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
                          cbi_path = NULL,
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
  
  cbi_stack <- cbi_clean(cbi_path = NULL, 
                         cbi_stack_path = "D:/GIS_Data/CBI_Sierra/CBI_1985_2024_ZeroFilling_Stack_New.tif",
                         save_path = NULL)
  
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
                                   metrics = c("lsm_c_pland", "lsm_c_ed"))
  } else {
    fire_lscp_out <- NULL
  }

  return(list(FireSeverity = fire_sev,
              FireMetrics = fire_buff_merge,
              FireLscp = fire_lscp_out))

} ## function closure

## Set the buffer size
buffSize <- 500

fire_sev21 <- aru_fire_prep(fire_prod = c("fire_severity"),
                            locs_from_cabio = FALSE,
                            custom_locs = locs,
                            survey_years = c(2021, 2022, 2023, 2024, 2025),
                            intervals = c("1-10"), #only interested in recent fire
                            id_col = "deployment_name",
                            buff_size = buffSize,
                            landscape_metrics = TRUE
)

## Save the R object for later
saveRDS(fire_sev21, file = here("Data/FireMets_ARU_21_25_AllUnitsByYears.RDS"))



