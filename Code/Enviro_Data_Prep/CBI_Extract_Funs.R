## -------------------------------------------------------------
##
## Script name: CBI Extraction Functions
##
## Script purpose: Script housing the function for CBI_Extraction.R
##
## Author: Spencer R Keyser
##
## Date Created: 2025-03-20
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
library(ggplot2)
library(here)
library(exactextractr)
library(landscapemetrics)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Chunking functions
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Input validation
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

validate_inputs <- function(fire_prod, locs_from_cabio, custom_locs, survey_years, ras_stack, intervals, allow_gaps = TRUE) {
  valid_products <- c("fire_severity", "time_since_fire", "fire_freq", "fire_ret_int")
  
  ## Validate arguments for products
  if(!locs_from_cabio && is.null(custom_locs)) {
    stop("No locations provided. Either set locs_from_cabio = TRUE or provide locations for custom_locs.")
  }
  if(!all(fire_prod %in% valid_products)) {
    stop(sprintf("Invalid fire product. Must be one of: %s", 
                 paste(valid_products, collapse = ", ")))
  }
  if(!is.numeric(survey_years)) {
    stop("Invalid survey years provided. Survey_years should be numeric.")
  }
  
  ## Validate raster stack
  if(!is.null(ras_stack)) {
    # Check layer names exist
    if(is.null(names(ras_stack))) {
      stop("Raster stack layers must be named")
    }
    
    ## Extract years from layer names
    years <- stringr::str_extract(names(ras_stack), "\\d{4}")
    years <- as.numeric(years)
    
    if(any(is.na(years))) {
      stop("Could not extract valid years from all layer names")
    }
    
    ## Check temporal ordering
    if(!all(diff(years) >= 0)) {
      stop("Raster stack layers are not in chronological order. Reorder raster data chronologically.")
    }
    
    ## Check for temporal gaps
    if(!allow_gaps && !all(diff(years) == 1)) {
      stop("Gaps detected in temporal sequence. Gaps in raster data are not handled downstream. Set your intervals accordingly.")
    }
    
    ## Interval validation
    if(!is.null(intervals) & is.character(intervals)) {
      ## Sort intervals based on first number
      intervals <- intervals[order(sapply(strsplit(intervals, "-"), 
                                          function(x) as.numeric(x[1])))]
      
      ## Extract interval ranges and create full sequences
      int_ranges <- lapply(strsplit(intervals, "-"), function(x) {
        if(length(x) == 1) return(as.numeric(x))
        seq(as.numeric(x[1]), as.numeric(x[2]))
      })
      
      ## Check for gaps between intervals
      full_seq <- sort(unique(unlist(int_ranges)))
      gaps <- diff(full_seq)
      if(any(gaps > 1)) {
        gap_starts <- full_seq[which(gaps > 1)]
        gap_ends <- full_seq[which(gaps > 1) + 1]
        gap_msg <- paste(sprintf("Gap found between %d and %d", 
                                 gap_starts, gap_ends), 
                         collapse = "\n")
        
        ## User choice
        choice <- menu(c("Continue with gaps in the intervals", 
                         "Stop function to fix gaps"),
                       title = paste("Interval gaps detected:\n", gap_msg, 
                                     "\nWhat would you like to do?"))
        if(choice == 1) {
          message("Continuing with gaps in intervals.")
        } else {
          stop("Function stopped. Please fix interval gaps.")
        }
      }
      
      ## Check if intervals exceed data range
      max_year <- max(years)
      min_year <- min(years)
      max_interval <- max(full_seq)
      
      if(any(survey_years - max_interval < min_year)) {
        warning(sprintf(
          "Some intervals exceed available fire data range: %d-%d",
          min_year, max_year
        ))
      }
      
      ## Check if survey years exceed data range
      if(any(survey_years > max_year + 1)) {
        stop(sprintf(
          "Survey years exceed available fire data range: %d-%d",
          min_year, max_year
        ))
      }
    }
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Trim CBI Stack depending on the survey years
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cbi_trim <- function(cbi_stack, survey_years){
  max_survey_year <- max(survey_years)
  keep_indices <- which(as.numeric(names(cbi_stack)) < max_survey_year)
  return(terra::subset(cbi_stack, keep_indices))
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: CAbioacoustics function
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load the ROI
cabio_loc_query <- function(years = survey_years){
  
  ## Connect to CAbioacoustics DB
  CAbioacoustics::cb_connect_db()
  
  ## View the tables housed in the DB
  DBI::dbListTables(conn = conn)
  
  ## Retrieve the study area
  ## Breaking function...removing until I hear back from Jay
  #roi <- CAbioacoustics::cb_get_spatial(layer_name = "sierra_study_area")
  
  ## Grab the ARU location from the CAbioacoustics package
  ## Pull the spatial coordinates
  aru_df <-
    tbl(conn, "acoustic_field_visits") |>
    filter(survey_year %in% years & study_type == 'Sierra Monitoring' ) |>
    select(
      id,
      study_type,
      group_id,
      visit_id,
      cell_id,
      unit_number,
      deployment_name,
      survey_year,
      deploy_date,
      utm_zone,
      utme,
      utmn
    ) |>
    collect()
  
  aru_locs <- aru_df |>
    filter(!is.na(utm_zone)) |> 
    group_split(utm_zone) |>
    map_dfr(cb_make_aru_sf) |> 
    mutate(Cell_Unit = stringr::str_remove(deployment_name, "G[0-9]+_V[0-9]+_"))
  
  ## Remove connection
  CAbioacoustics::cb_disconnect_db()
  
  ## Return the object
  return(aru_locs)
  
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Prepare locations
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prep_locs <- function(locs_from_cabio,
                      custom_locs,
                      survey_years,
                      x_col,
                      y_col,
                      .crs){
  if(locs_from_cabio){
    aru_locs <- cabio_loc_query(years = survey_years)
  } else if (!is.null(custom_locs)) {
    aru_locs <- custom_locs
  } else {
    stop("No locations provided. Either set locs_from_cabio = TRUE or provide locations for custom_locs.")
  }
  
  ## Create SF objects from points
  if(!any(class(aru_locs) %in% "sf")){
    
    aru_locs <- st_as_sf(aru_locs, 
                         coords = c(x_col, 
                                    y_col), 
                         crs = .crs)
  }
  return(aru_locs)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: CBI Cleaning Function
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cbi_clean <- function(cbi_path = NULL, cbi_stack_path = NULL, save_path = NULL){
  
  if (!is.null(cbi_stack_path) && !file.exists(cbi_stack_path)) {
    stop("Supplied cbi_stack_path does not exist: ", cbi_stack_path)
  }
  if (!is.null(cbi_path) && !dir.exists(cbi_path)) {
    stop("Supplied cbi_path does not exist: ", cbi_path)
  }
  
  if(!is.null(cbi_stack_path)){
    
    message("CBI rasters are fixed and exist in directory. Loading the fixed stack.")
    
    cbi_stack <- terra::rast(cbi_stack_path)
    
  } else if(!is.null(cbi_path)) {
    
    message("No CBI stack path supplied. Proceeding with yearly rasters.")
    
    ## Find files
    cbi_files <- list.files(cbi_path, 
                            full.names = T, 
                            pattern = "(cbi_cat_)(\\d{4})(*.tif$)")
    
    if(length(cbi_files) == 0){
      stop("No files found in ", cbi_path, " matching pattern.")
    }
    
    ## Call function to stack if need be
    cbi_stack <- cbi_stack_fun(cbi_files)
    
    ## Reclassify the NAs to zero to fill in the rasters
    cbi_stack[is.na(cbi_stack)] <- 0
    
    if(nlyr(cbi_stack) == length(cbi_files)){
      names(cbi_stack) <- str_extract(cbi_files, "\\d{4}")
    }
    
    if(!is.null(save_path)){
      message("Saving processed stack to: ", save_path)
      terra::writeRaster(cbi_stack, filename = save_path)
    }
    
  } else { 
    
    stop("Both CBI data locations are null! Provide either cbi_stack_path or cbi_path.") 
  
  }
  
  return(cbi_stack)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: CBI Stacking Function
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cbi_stack_fun <- function(cbi_files){
  tryCatch({
    res <- terra::rast(cbi_files)
    return(res)
  }, error = function(e) {
    ## Error message
    message("Error:", e$message)
    
    ## if error is different spatial extents
    if(stringr::str_detect(e$message, "extents do not match")){
      
      ## Print a message
      cat("Extents do not match. Attempting to resolve. \n This could take a minute...")
      
      ## Fix extents
      cbi_list <- lapply(cbi_files, rast)
      ## Extents for all rasters
      l_extent <- lapply(cbi_list, ext)
      ## Make all extents bboxes
      cbi_bbox <- lapply(l_extent, function(x) vect(x))
      ## Union to find the largest extent
      big.ext <- Reduce(terra::union, cbi_bbox)
      ## Set the CRS
      crs(big.ext) <- crs(cbi_list[[1]])
      ## Extend all rasters to the same extent
      cbi_list <- lapply(cbi_list, function(x) extend(x, big.ext))
      ## Restack the rasters
      cbi_stack <- try({
        res <- do.call(c, cbi_list)
      }, silent = T)
      
      if(inherits(cbi_stack, "try-error")){
        stop("Check raster data. Could not stack.")
      } else {
        cat("Successful raster matching.")
        return(res)
      } #if error after alignment
      
    } #if error
    
  })}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire year interval stacking
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int_ras_fun <- function(ras_stack,
                        intervals,
                        sum_int,
                        locations){
  
  if(!is.null(intervals) & class(intervals) == "character"){
    intervals <- intervals[order(sapply(strsplit(intervals, "-"), 
                                        function(x) as.numeric(x[1])))]
    int_len <- length(intervals)
    seq_list <- vector(mode = "list", length = int_len)
    int <- 0
    while(int < int_len){
      int <- int + 1
      int.t <- strsplit(intervals[int], "-")[[1]]
      if(length(int.t) > 1){
        s_seq <- int.t[[1]]
        e_seq <- int.t[[2]]
        seq_list[[int]] <- seq(s_seq, e_seq, by = 1)
      } else {seq_list[[int]] <- as.numeric(int.t)}
    }
    
    ## Check if the seq list is given in dates or integer years
    if(!any(sort(unlist(seq_list)) %in% names(ras_stack))){
      aru_years <- unique(locations$survey_year)
      if(length(aru_years) > 1){
        
        ## Make sure all the years are integers
        years <- as.integer(aru_years)
        ## Find the number of years
        num_years <- length(years)
        ## Replicate the sequence list n times
        seq_list <- rep(list(seq_list), num_years)
        ## for each year replace the values to match the names of the fire data
        for(i in 1:num_years){
          seq_tmp <- seq_list[[i]]
          ## Back calculate the expected dates from the sequence of numbers
          ## Sort the unlisted seq_list vec
          seq_vec <- sort(unlist(seq_tmp))
          ## Find the year of interest
          years <- as.integer(aru_years[i])
          ## Find all years before the year of interest
          fire_years <- years - seq_vec
          ## Update the seq_list to the dates of the fires
          seq_list[[i]] <- lapply(seq_list[[i]], function(x) fire_years[x])
        }
      }
      
      if(length(aru_years) == 1){
        ## Back calculate the expected dates from the sequence of numbers
        ## Sort the unlisted seq_list vec
        seq_vec <- sort(unlist(seq_list))
        ## Find the year of interest
        years <- as.integer(aru_years)
        ## Find all years before the year of interest
        fire_years <- years - seq_vec
        ## Update the seq_list to the dates of the fires
        seq_list <- lapply(seq_list, function(x) fire_years[x])
      }
    } # fix year naming
    
    if(length(aru_years) == 1){
      ## Create the spatial products with the sequences for 1 year
      ras_list <- vector(mode = "list", length = int_len)
      for(i in 1:length(seq_list)){
        message("This is a single year workflow. Set by survey_years at the top-level.")
        message(paste("Finding", sum_int, "value for Year:", aru_years, "For interval", i, "of", int_len))
        seq_temp <- seq_list[[i]]
        ras_tmp <- ras_stack[[which(names(ras_stack) %in% seq_list[[i]])]]
        if(nlyr(ras_tmp) < 2){
          ras_list[[i]] <- ras_tmp
          names(ras_list[[i]]) <- intervals[i]  # Use interval name directly
        } else {
          ras_list[[i]] <- app(ras_tmp, fun = sum_int)
          names(ras_list[[i]]) <- intervals[i]  # Use interval name directly
        }
      }
      ras_int <- list(do.call(c, ras_list))
      names(ras_int[[1]]) <- intervals  # Name the combined stack layers
    }
    
    if(length(aru_years) > 1){
      ## Create the spatial products with the sequences for 1 year
      ras_list <- vector(mode = "list", length = int_len)
      ras_list <- rep(list(ras_list), num_years)
      for(i in 1:num_years){
        for(j in 1:int_len){
          print(paste("Finding", sum_int, "value for Year:", aru_years[i], "For interval", j, "of", int_len))
          seq_temp <- unlist(seq_list[[i]][j])
          ras_tmp <- ras_stack[[which(names(ras_stack) %in% seq_temp)]]
          if(nlyr(ras_tmp) < 2){
            ras_list[[i]][[j]] <- ras_tmp
            names(ras_list[[i]][[j]]) <- intervals[j]  # Use interval name directly
          } else {
            ras_list[[i]][[j]] <- app(ras_tmp, fun = sum_int)
            names(ras_list[[i]][[j]]) <- intervals[j]  # Use interval name directly
          }
        }
        
        ras_int <- vector(length(ras_list), mode = "list")
        for(k in 1:length(ras_list)){
          ras_int[[k]] <- do.call(c, ras_list[[k]])
        }
      }
    }
    
  } #FULL IF STATEMENT
  
  ## Return the stacked intervals
  return(ras_int)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Interval Process
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Process intervals
process_intervals <- function(cbi_stack, intervals, locations, int_fun = "max") {
  if(is.null(intervals)) {
    message("No intervals specified. Skipping interval binning. Returning 'cbi_stack'.")
    return(cbi_stack)
  } else {
    message(paste("Intervals specified. Binning using", int_fun, "function."))
    int_ras <- int_ras_fun(ras_stack = cbi_stack,
                           intervals = intervals,
                           sum_int = int_fun,
                           locations = locations)
    return(int_ras)
  }
}
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Extract fire severity
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_fire_sev <- function(cbi_sev, aru_locs, buff_size, id_col){
  #if("fire_severity" %in% fire_prod){
  
    if(!is.null(buff_size)){
    message(paste("Buffer size detected. Extracting from", buff_size, "buffer."))
    aru_buffer <- st_buffer(aru_locs,
                            dist = buff_size)
    
    ## Check the number of unique years
    if(length(unique(aru_locs$survey_year)) == length(cbi_sev)){
      message("Matching fire intervals to the unique years in the data. Proceed.")
    } else { stop("Mismatching number of years for extraction.") }
    
    ## If the years are more than 1
    #if(length(unique(aru_locs$survey_year)) > 1){
      fire_sev <- vector(mode = "list", length = length(unique(aru_locs$survey_year)))
      yrs_tmp <- unique(aru_locs$survey_year)
      
      for(i in 1:length(unique(aru_buffer$survey_year))){
        fire_sev[[i]] <- exactextractr::exact_extract(cbi_sev[[i]],
                                                      aru_buffer |> filter(survey_year == yrs_tmp[i]),
                                                      fun = c("mean", "stdev"),
                                                      append_cols = id_col) |> 
          mutate(SurvExtYr = unique(aru_buffer$survey_year)[i])
      }
      
      ## Bind the multi-year data together
      fire_sev <- do.call(rbind, fire_sev)
      
      ## Rename to more aesthetically pleasing names
      colnames(fire_sev)[colnames(fire_sev) != id_col] <- paste0("Fire_Sev_", gsub("[[:punct:]]", "_", colnames(fire_sev)[colnames(fire_sev) != id_col]))
    #}
  } else {
    ## This does point level extractions
    fire_sev <- terra::extract(cbi_sev,
                               aru_locs,
                               bind = T) |> st_as_sf()
    
  }
  
# } else {
#   fire_sev <- NULL
# }
  return(fire_sev)
} # end function

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Landscape metrics function
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Landscape Metrics
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function should estimate specific landscape metrics for the
## interval-based raster data for CBI
fire_lscp_fun <- function(ras_int, 
                          ras_stack, 
                          locs,
                          years,
                          intervals,
                          buff_size,
                          id_col,
                          metrics){
  
  for(i in 1:length(years)){
    
    print(paste("Year:", years[i]))
    
    ## Pull the ARUs for a specific year
    locs_tmp <- locs |> 
      filter(survey_year == years[i])
    
    ## Is there a buffer size?
    if(!is.null(buff_size)){ #| is.null(hex_size)){
      
      ## Check to make sure sf
      if(!is.null(id_col) & any(str_detect(class(locs_tmp), "sf"))){
        id.v <- as.vector(unlist(st_drop_geometry(locs_tmp[, id_col])))
      } else {
        id.v <- NULL
      }
      
      ## Individual years versus intervals for the fire data
      if(is.null(intervals)){ 
        message("Landscapemetrics done for individual years.")
        ras_lcpc <- ras_stack 
      } else { 
        ras_lcpc <- ras_int[[i]]
        message("Landscapemetrics done for intervals.")
      }
      
      ## Calculate LSCP
      lsm <- landscapemetrics::sample_lsm(landscape = ras_lcpc,
                                          y = locs_tmp,
                                          plot_id = id.v,
                                          shape = "circle",
                                          size = buff_size,
                                          what = metrics,
                                          #metric = "core",
                                          #level = "patch",
                                          return_raster = F,
                                          #type = "diversity metric",
                                          #classes_max = 3,
                                          verbose = FALSE,
                                          progress = T)
      
      ## Deal with some of the spatial inaccuracies by taking means
      lsm <- lsm |> 
        group_by(plot_id, class, layer, metric, level) |> 
        summarise(value = mean(value)) |> 
        ungroup()
      
      ## Fix layer names
      lyr_map <- data.frame(layer = 1:nlyr(ras_lcpc), lyr_name = names(ras_lcpc))
      
      ## Merge these
      lsm <- left_join(lsm, lyr_map)
      
      ## level of lsm
      lsm_lev <- unique(lsm$level)
      
      if(any(lsm_lev == "landscape")){
        lsm_l <- lsm |> 
          filter(level == "landscape") |> 
          select(plot_id, level, lyr_name, metric, value) |> 
          tidyr::pivot_wider(names_from = c(metric, lyr_name), 
                             values_from = value,
                             names_glue = "{lyr_name}_{metric}_l") |> 
          select(!level) |>
          mutate(Year = years[i])
        
      }
      
      if(any(lsm_lev == "class")){
        
        fire_class <- c("Unburned", "Low_Sev", "Mod_Sev", "High_Sev")
        
        lsm_c <- lsm |> 
          filter(level == "class") |> 
          ## This handles the shifting spatial points for each
          ## ARU across the years. A better way would be to filter
          ## the ARUs actually surveyed and then bridge the remaining
          ## ARUs unsampled in that year but sampled at another point
          ## for generating a square matrix for model fitting.
          group_by(plot_id, class, lyr_name, metric, level) |> 
          summarise(value = mean(value)) |> 
          ungroup() |> 
          ## Original functionality
          mutate(FireClass = case_when(class == 0 ~ "Unburned",
                                       class == 1 ~ "Low_Sev",
                                       class == 2 ~ "Mod_Sev",
                                       class == 3 ~ "High_Sev")) |> 
          select(plot_id, level, lyr_name, FireClass, metric, value) |> 
          tidyr::pivot_wider(names_from = c(metric, FireClass, lyr_name), 
                             values_from = value,
                             names_glue = "{FireClass}_{lyr_name}_{metric}_c") |> 
          select(!level) |> 
          mutate(Year = years[i])
        
      }
      
      if(exists("lsm_l") & exists("lsm_c")){ lsm <- full_join(lsm_l, lsm_c) } 
      if (exists("lsm_l") & !exists("lsm_c")) { lsm <- lsm_l }  
      if (!exists("lsm_l") & exists("lsm_c")) { lsm <- lsm_c }
    }
    
    ## Grow the DF for multiple years
    if(i == 1){
      lsm_out <- lsm
    } else {
      ## dplyr::bind_rows has preffered behavior to rbind
      lsm_out <- dplyr::bind_rows(lsm_out, lsm)
    }
    
  }
  ## Return
  return(lsm_out)
  
}

## -------------------------------------------------------------
##
## Begin Section: Fire Functions
##
## -------------------------------------------------------------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Year of most recent fire
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Potential improvements
## For each year assign the fires a simple year estimate
## 40 rasters with values (Year, NA)
## Trim the rasters based on the survey year max Year Stack < Survey Year
## DiffRasStack = Survey Year - Each image
## Min(DiffRasStacl) == Most Recent Fire rel survey year
## NA -> Min year 
## Repeat each survey year

## Fire variables
## Time since fire
time_to_most_recent_fire <- function(cell_values, years) {

  # Ensure years are numeric and sorted
  years <- sort(as.numeric(years))
  max_year <- max(years)
  
  # Check if there are any fires (values > 0)
  if (any(cell_values > 0, na.rm = TRUE)) {
    # Find the most recent fire year index
    fire_indices <- which(cell_values > 0)
    if(length(fire_indices) > 0) {
      most_recent_fire_index <- max(fire_indices)
      most_recent_fire_year <- years[most_recent_fire_index]
      time_since_fire <- max_year - most_recent_fire_year
      
      return(time_since_fire)
    }
  }
  
  # If no fires found, return the full time period
  return(max_year - min(years))
}

# ## Optimized time_to_most_recent_fire function
# create_last_fire_raster <- function(cbi_stack, temp.dir = NULL) {
#   
#   if(is.null(temp.dir)){
#     message("Custom fire metrics running without temp dir. \nFunction might fail without adequate RAM.")
#   }
#   
#   years <- as.numeric(names(cbi_stack))
#   
#   # Create year rasters (each layer contains its year value where fires occurred)
#   year_rasters <- list()
#   # for(i in 1:nlyr(cbi_stack)) {
#   #   year_rasters[[i]] <- terra::ifel(cbi_stack[[i]] > 0, years[i], NA)
#   # }
#   
#   for(i in 1:nlyr(cbi_stack)) {
#     ## Classification matrix
#     rcl_matrix <- matrix(c(
#       -Inf, 0, NA,      
#       0, Inf, years[i]
#     ), ncol = 3, byrow = TRUE)
#     
#     year_rasters[[i]] <- terra::classify(cbi_stack[[i]], rcl_matrix, include.lowest = FALSE)
#     
#     ## Clean in between iterations
#     gc()
#   }
#   
#   # Stack and find maximum (most recent) year
#   year_stack <- do.call(c, year_rasters)
#   #last_fire_year <- terra::app(year_stack, fun = max, na.rm = TRUE)
#   
#   if(!is.null(temp.dir)){
#     
#     ## Write the stack to avoid memory bog
#     print(paste0("Writing intermediate raster to temp directory specified: ", temp.dir))
#     writeRaster(year_stack, filename = paste0(temp.dir, "FiresAsYears.tif"), overwrite = TRUE)
#     
#     ## Remove and clean
#     remove(year_stack, year_rasters)
#   } else {
#     return(year_stack)
#   }
#   gc()
#   
# }

## Optimized time_to_most_recent_fire function
create_last_fire_raster <- function(cbi_stack, temp.dir = NULL, hsf_only = FALSE) {
  
  years <- as.numeric(names(cbi_stack))
  
  if(is.null(temp.dir)){
    message("Custom fire metrics running without temp dir. \nFunction might fail without adequate RAM.")
    
    # MEMORY-INTENSIVE APPROACH (your original)
    year_rasters <- list()
    
    for(i in 1:nlyr(cbi_stack)) {
      ## Classification matrix
      if(!hsf_only){
        if(i == 1){message("Using all CBI classes.")}
        rcl_matrix <- matrix(c(
          -Inf, 0, NA,      
          0, Inf, years[i]
        ), ncol = 3, byrow = TRUE)
      } else {
        if(i == 1){message("High-severity fire only.")}
        rcl_matrix <- matrix(c(
          -Inf, 2, NA,      
          2, Inf, years[i]
        ), ncol = 3, byrow = TRUE)
      }
      
      year_rasters[[i]] <- terra::classify(cbi_stack[[i]], rcl_matrix, include.lowest = FALSE)
      
      ## Clean in between iterations
      gc()
    }
    
    # Stack and find maximum (most recent) year
    year_stack <- do.call(c, year_rasters)
    return(year_stack)
    
  } else {
    message("Using file-based approach with temp directory: ", temp.dir)
    
    # Create directory if needed
    dir.create(temp.dir, showWarnings = FALSE, recursive = TRUE)
    
    # Pre-allocate file paths
    temp_files <- character(nlyr(cbi_stack))
    
    for(i in 1:nlyr(cbi_stack)) {
      message(paste("Processing layer", i, "of", nlyr(cbi_stack), "- Year:", years[i]))
      
      ## Classification matrix
      if(!hsf_only){
        if(i == 1){message("Using all CBI classes.")}
        rcl_matrix <- matrix(c(
          -Inf, 0, NA,      
          0, Inf, years[i]
        ), ncol = 3, byrow = TRUE)
      } else {
        if(i == 1){message("High-severity fire only.")}
        rcl_matrix <- matrix(c(
          -Inf, 2, NA,      
          2, Inf, years[i]
        ), ncol = 3, byrow = TRUE)
      }
      
      # Create file path for this layer
      temp_files[i] <- file.path(temp.dir, paste0("fire_year_", years[i], ".tif"))
      
      # Save directly to file (avoids memory buildup)
      terra::classify(cbi_stack[[i]], rcl_matrix, 
                      filename = temp_files[i],
                      include.lowest = FALSE, 
                      overwrite = TRUE)
    }
    return(temp_files)
  }
  
  gc()
}


time_since_fire <- function(last_fire_raster, survey_year) {
  if(is.character(class(last_fire_raster))){
    message("last_fire_raster variable is file paths. \nReading files as stack.")
    last_fire_raster <- rast(last_fire_raster)
  }
  ## Use trimming to restrict to the current survey year of interest
  message("Trimming")
  cbi_tsf <- cbi_trim(cbi_stack = last_fire_raster, survey_years = survey_year)
  
  ## Find the difference between survey year and all fires
  message("TSF for all years")
  cbi_tsf <- survey_year - cbi_tsf
  
  ## Find the minimum difference across all layers
  message("Calculating Mininmum TSF")
  cbi_tsf <- terra::app(cbi_tsf, fun = min, na.rm = TRUE)
  
  ## Handle areas with no fires - set to full time period -1 year to avoid confounding
  ## sites with older fire and sites with no fires 
  message("Calculating TSF")
  max_tsf <- survey_year - (min(as.numeric(names(last_fire_raster)))-1)
  
  message("Filling NAs with max year between survey date and last year of CBI.")
  time_since_fire_raster <- terra::ifel(is.na(cbi_tsf), max_tsf, cbi_tsf)
  
  return(time_since_fire_raster)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire frequency
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Fire freq
# Create the function to calculate time since the most recent event for each cell
fire_freq_calc <- function(cell_values, years) {
  
  # Check for fires, handling NAs
  has_fires <- any(cell_values > 0, na.rm = TRUE)
  
  if (is.na(has_fires)) {
    # All values are NA
    return(0)  # or return(NA) if you prefer
  }
  
  if (has_fires) {
    fire_years <- years[which(cell_values > 0)]
    num_fires <- length(fire_years)
    num_years <- length(years)
    
    fire_freq <- num_fires / num_years
    
  } else {
    fire_freq <- 0
  }
  
  return(fire_freq)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire Return Interval
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create the function to calculate fire return interval
## Return mean difference in time between fires (assumes yearly data)
## Returns NA under two conditions: Less than 2 fires produce NA
fire_return_int <- function(cell_values, years){
  # Check if there are any 1s in the cell
  if (any(cell_values > 0)) {
    
    # Find the most recent occurrence of 1 (event) - we look for the last occurrence of 1
    fire_years <- as.numeric(years[which(cell_values > 0)])
    
    if(length(fire_years) > 1){
      # Calculate the average time in-between fire events
      f_int <- diff(fire_years)
      ret_int <- mean(f_int)
    } else {
      ret_int <- NA
    }
  } else {
    # If no event occurred, return NA
    #ret_int <- length(years)
    ret_int <- NA
  }
  
  return(ret_int)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Wrapper for custom fire variables
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## WIP for fire freq and fire return int with multi-year data
## TSF is working

## Deal with issues in the shifting years
create_fire_metrics <- function(cbi_stack, 
                                fire_prod, 
                                survey_years, 
                                intervals = NULL, 
                                inter.dir = NULL,
                                hsf_only = FALSE) {
  
  if(!is.null(intervals)){
    message("Intervals are set. Time since fire, fire frequency, and fire return interval don't accept intervals...output will be generated from single years for all CBI years provided.")
  }
  
  if(length(fire_prod[!fire_prod %in% "fire_severity"]) == 0) {
    message("No custom fire metrics requested.")
    return(list())
  }
  
  ## Filter products
  custom_products <- fire_prod[!fire_prod %in% "fire_severity"]
  
  metrics_by_year <- list()
  
  ## Map functions except tsf for now
  metric_functions <- list(
    "fire_freq" = fire_freq_calc,
    "fire_ret_int" = fire_return_int
  )
  
  ## Pre-compute year stack
  last_fire_raster <- NULL
  if("time_since_fire" %in% custom_products) {
    message("Pre-computing last fire raster (one-time calculation)...")
    if(!is.null(inter.dir)){
      message("Creating intermediate files.")
      last_fire_raster <- create_last_fire_raster(cbi_stack, 
                                                  temp.dir = inter.dir,
                                                  hsf_only = hsf_only)
    } else {
      message("Proceeding with memory intensive process.")
      last_fire_raster <- create_last_fire_raster(cbi_stack, 
                                                  temp.dir = inter.dir,
                                                  hsf_only = hsf_only)
    }
  }
  
  ## Iterate through each survey year
  for(survey_year in survey_years) {
    message(paste("Creating fire metrics for survey year:", survey_year))
    
    year_metrics <- list()
    
    ## Handle time_since_fire with survey-year-specific approach
    if("time_since_fire" %in% custom_products) {
      message(paste("  - Calculating time_since_fire"))
      
      # ## Trim stack to only years up to survey year to avoid
      # ## misaligned years
      # all_years <- as.numeric(names(cbi_stack))
      # valid_years_idx <- which(all_years <= survey_year)
      # 
      # if(length(valid_years_idx) == 0) {
      #   warning(paste("No fire data available for years up to survey year:", survey_year))
      #   next
      # }
      # 
      # ## Create trimmed stack for this survey year
      # cbi_up_to_survey <- cbi_stack[[valid_years_idx]]
      
      ## Create last fire raster using only years up to survey year
      # last_fire_raster <- create_last_fire_raster(cbi_up_to_survey)
      
      year_metrics[["time_since_fire"]] <- time_since_fire(
        last_fire_raster = last_fire_raster, 
        survey_year = survey_year
      )
      names(year_metrics[["time_since_fire"]]) <- paste0("time_since_fire_", survey_year)
      
      ## Clean up
      gc()
    }
    
    ## Handle other metrics that need trimmed stacks and terra::app
    other_products <- custom_products[!custom_products %in% "time_since_fire"]
    
    if(length(other_products) > 0) {
      # Trim stack for this specific survey year (only for other metrics)
      cbi_trimmed <- cbi_trim(cbi_stack = cbi_stack, survey_years = survey_year)
      
      if(is.null(cbi_trimmed)) {
        warning(paste("No fire data available for survey year:", survey_year))
      } else {
        # Create each other requested metric for this survey year
        for(metric in other_products) {
          if(metric %in% names(metric_functions)) {
            
            message(paste("  - Calculating", metric))
            
            year_metrics[[metric]] <- terra::app(
              cbi_trimmed,
              fun = function(x) metric_functions[[metric]](
                cell_values = x,
                years = as.numeric(names(cbi_trimmed))
              ),
              cores = 1
            )
            
            # Name the raster layer
            names(year_metrics[[metric]]) <- paste0(metric, "_", survey_year)
            gc()
          }
        }
        
        # Clean up trimmed stack
        rm(cbi_trimmed)
        gc()
      }
    }
    
    # Store metrics for this survey year
    metrics_by_year[[as.character(survey_year)]] <- year_metrics
  }
  
  return(metrics_by_year)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire metrics with buffer
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extract_metrics <- function(metrics_list, aru_locs, survey_years, buff_size = NULL, id_col) {
  
  ## Check buff_size
  if (!is.null(buff_size)) {
    extraction_geom <- sf::st_buffer(aru_locs, dist = buff_size)
    fun_type <- "mean"
  } else {
    extraction_geom <- aru_locs
    fun_type <- "mean"
  }
  
  ## Initialize list to collect all metrics across years
  all_metrics <- list()
  
  ## Loop through years
  for(year in survey_years) {
    
    ## Get yearly metrics
    year_metrics <- metrics_list[[as.character(year)]]
    
    ## Filter geometry to this year's locations
    extraction_geom_year <- extraction_geom[extraction_geom$survey_year == year, ]
    
    ## Skip no data
    if(nrow(extraction_geom_year) == 0 || is.null(year_metrics)) next
    
    ## Extract each metric for this year
    for(metric_name in names(year_metrics)) {
      
      result <- exactextractr::exact_extract(
        year_metrics[[metric_name]],
        extraction_geom_year,
        fun = fun_type,
        append_cols = id_col
      )
      
      ## Create unique metric name with year
      final_metric_name <- paste0(metric_name, "_", year)
      names(result)[names(result) == fun_type] <- final_metric_name
      
      ## Store in collection
      all_metrics[[final_metric_name]] <- result
    }
  }
  
  ## Combine all metrics
  final_result <- Reduce(function(x, y) merge(x, y, by = id_col, all = TRUE), all_metrics)
  return(final_result)
}