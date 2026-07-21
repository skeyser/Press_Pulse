## ----------------------------------------------------------
##
## Script name: DCM Viz Functions (Updated for JJA/MAM flexibility)
##
## Script purpose: Dynamic community model visualization with seasonal climate support
##
## Author: Spencer R Keyser
##
## Date Created: 2026-04-06
## Updated: [Current Date] - Added JJA/MAM climate season flexibility
##
## Email: skeyser@wisc.edu
##
## Github: https://github.com/skeyser
##
## -----------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(data.table)
library(tibble)
library(coda)
library(MCMCvis)
library(tidyr)
library(tidybayes)

## -----------------------------------------------------------

# Helper function to get variable configuration based on climate season
get_var_config <- function(climate_season = "jja") {
  
  # Convert to lowercase for consistency
  season <- tolower(climate_season)
  
  if(!season %in% c("jja", "mam")) {
    stop("climate_season must be either 'jja' or 'mam'")
  }
  
  # Create season-specific variable names
  temp_anom <- paste0("tmax_anom_", season)
  precip_anom <- paste0("p_anom_", season)
  temp_trend <- paste0("tmax_trend_", season)
  tmin_trend <- paste0("tmin_trend_", season)
  precip_trend <- paste0("prec_trend_", season)
  
  var_config <- list(
    beta = list(
      vars = c("btmax_annual", "bprec_annual", "cc"),
      data_cols = c("btmax_annual", "bprec_annual", "cc"),
      scaling = c("z", "z", "z"),
      coef_pattern = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6"),
      has_quadratic = TRUE
    ),
    gamma = list(
      vars = c(temp_anom, precip_anom, "hsf_pland1_10", temp_trend, tmin_trend, precip_trend),
      data_cols = c(temp_anom, precip_anom, "hsf_pland1_10", temp_trend, tmin_trend, precip_trend),
      scaling = c("z_no_center", "z_no_center", "z", "raw", "raw", "raw"),
      coef_pattern = c("gamma1", "gamma2", "gamma3", "gamma6", "gamma7", "gamma8"),
      interactions = list(
        list(vars = c(temp_anom, "hsf_pland1_10"), coef = "gamma4"),
        list(vars = c(precip_anom, "hsf_pland1_10"), coef = "gamma5"),
        list(vars = c(temp_trend, "hsf_pland1_10"), coef = "gamma9"),
        list(vars = c(tmin_trend, "hsf_pland1_10"), coef = "gamma10"),
        list(vars = c(precip_trend, "hsf_pland1_10"), coef = "gamma11")
      ),
      has_quadratic = FALSE
    ),
    eps = list(
      vars = c(temp_anom, precip_anom, "hsf_pland1_10", temp_trend, tmin_trend, precip_trend),
      data_cols = c(temp_anom, precip_anom, "hsf_pland1_10", temp_trend, tmin_trend, precip_trend),
      scaling = c("z_no_center", "z_no_center", "z", "raw", "raw", "raw"),
      coef_pattern = c("eps1", "eps2", "eps3", "eps6", "eps7", "eps8"),
      interactions = list(
        list(vars = c(temp_anom, "hsf_pland1_10"), coef = "eps4"),
        list(vars = c(precip_anom, "hsf_pland1_10"), coef = "eps5"),
        list(vars = c(temp_trend, "hsf_pland1_10"), coef = "eps9"),
        list(vars = c(tmin_trend, "hsf_pland1_10"), coef = "eps10"),
        list(vars = c(precip_trend, "hsf_pland1_10"), coef = "eps11")
      ),
      has_quadratic = FALSE
    )
  )
  
  return(var_config)
}

## Updated marginal effects function
dynOcc_margeff <- function(var = "tmax_anom_jja",
                           model_component = "gamma",
                           param_compile = param.compile,
                           raw_data = bdata,
                           climate_season = "jja",  # NEW PARAMETER
                           level = "species",
                           pred_point = 100){
  
  # Get season-appropriate variable configuration
  var_config <- get_var_config(climate_season)
  
  # Get configuration
  config <- var_config[[model_component]]
  if(is.null(config)) {
    stop("Model component not recognized: ", model_component)
  }
  
  # Check if variable of interest is valid
  if(!var %in% config$vars) {
    stop("Variable '", var, "' not available for component '", model_component, 
         "' in season '", climate_season, "'. Available: ", paste(config$vars, collapse = ", "))
  }
  
  # Function to apply scaling transformations
  apply_scaling <- function(x, scale_type) {
    switch(scale_type,
           "z" = scale(x)[,1],
           "z_no_center" = scale(x, center=FALSE)[,1],
           "raw" = x,
           x)
  }
  
  # Process ALL variables - calculate their means and scaling
  all_var_data <- list()
  var_scaling_info <- list()
  original_sequence <- NULL  
  
  # Define which variables are dynamic (time-varying)
  season <- tolower(climate_season)
  dynamic_vars <- c(paste0("tmax_anom_", season), paste0("p_anom_", season), "hsf_pland1_10")
  
  for(i in 1:length(config$vars)) {
    var_name <- config$vars[i]
    data_col <- config$data_cols[i] 
    scale_type <- config$scaling[i]
    
    # Get raw data - check if it's a dynamic variable
    if(data_col %in% dynamic_vars) {
      if(is.list(raw_data[[data_col]])){
        raw_data[[data_col]] <- unlist(raw_data[[data_col]])
      }
      raw_values <- as.vector(raw_data[[data_col]])
    } else {
      raw_values <- raw_data[[data_col]]
    }
    
    # Store original stats for back-transformation if needed
    var_scaling_info[[var_name]] <- list(
      mean_orig = mean(raw_values, na.rm = TRUE),
      sd_orig = sd(raw_values, na.rm = TRUE),
      scale_type = scale_type
    )
    
    # Apply scaling transformation
    scaled_values <- apply_scaling(raw_values, scale_type)
    
    if(var_name == var) {
      # Variable of interest - create sequence across its range
      var_range <- range(scaled_values, na.rm = TRUE)
      all_var_data[[var_name]] <- seq(var_range[1], var_range[2], length.out = pred_point)
      
      # Also store original scale sequence for plotting
      orig_range <- range(raw_values, na.rm = TRUE)
      original_sequence <- seq(orig_range[1], orig_range[2], length.out = pred_point)
      
    } else {
      # All other variables - set to their means
      all_var_data[[var_name]] <- switch(scale_type,
                                         "z" = 0,
                                         "z_no_center" = 0,
                                         "raw" = mean(raw_values, na.rm = TRUE),
                                         mean(scaled_values, na.rm = TRUE))
    }
  }
  
  # Filter the dataframe to the components of interest
  mod.out <- param_compile |> 
    filter(str_detect(variable, model_component)) |> 
    select(sp.ind, variable, mean) |> 
    pivot_wider(names_from = variable, values_from = mean)
  
  # Get number of species
  nspec <- length(unique(mod.out$sp.ind))
  
  # Create the receiving structure
  out <- vector("list", nspec)
  
  # Calculate predictions for each species
  for(s in 1:nspec){
    if(model_component == "beta"){
      lin_pred <- mod.out$beta0[s] + 
        mod.out$beta1[s] * all_var_data$btmax_annual + 
        mod.out$beta2[s] * all_var_data$btmax_annual^2 + 
        mod.out$beta3[s] * all_var_data$bprec_annual + 
        mod.out$beta4[s] * all_var_data$bprec_annual^2 +
        mod.out$beta5[s] * all_var_data$cc +
        mod.out$beta6[s] * all_var_data$cc^2
      p <- plogis(lin_pred)
      p.df <- data.frame(Predicted = p, Value = all_var_data[[var]], 
                         Value_original = original_sequence,
                         Variable = var, SpInd = s)
      out[[s]] <- p.df
      
    } else if (model_component %in% c("gamma", "eps")){
      # Get the season-specific variable names
      temp_anom <- paste0("tmax_anom_", season)
      precip_anom <- paste0("p_anom_", season)
      temp_trend <- paste0("tmax_trend_", season)
      tmin_trend <- paste0("tmin_trend_", season)
      precip_trend <- paste0("prec_trend_", season)
      
      # Build linear predictor
      intercept_col <- paste0(model_component, "0")
      coef_prefix <- model_component
      
      lin_pred <- mod.out[[intercept_col]][s] + 
        mod.out[[paste0(coef_prefix, "1")]][s] * all_var_data[[temp_anom]] + 
        mod.out[[paste0(coef_prefix, "2")]][s] * all_var_data[[precip_anom]] + 
        mod.out[[paste0(coef_prefix, "3")]][s] * all_var_data[["hsf_pland1_10"]] + 
        mod.out[[paste0(coef_prefix, "4")]][s] * all_var_data[[temp_anom]] * all_var_data[["hsf_pland1_10"]] +
        mod.out[[paste0(coef_prefix, "5")]][s] * all_var_data[[precip_anom]] * all_var_data[["hsf_pland1_10"]] +
        mod.out[[paste0(coef_prefix, "6")]][s] * all_var_data[[temp_trend]] +
        mod.out[[paste0(coef_prefix, "7")]][s] * all_var_data[[tmin_trend]] +
        mod.out[[paste0(coef_prefix, "8")]][s] * all_var_data[[precip_trend]] +
        mod.out[[paste0(coef_prefix, "9")]][s] * all_var_data[[temp_trend]] * all_var_data[["hsf_pland1_10"]] +
        mod.out[[paste0(coef_prefix, "10")]][s] * all_var_data[[tmin_trend]] * all_var_data[["hsf_pland1_10"]] +
        mod.out[[paste0(coef_prefix, "11")]][s] * all_var_data[[precip_trend]] * all_var_data[["hsf_pland1_10"]]
      
      p <- plogis(lin_pred)
      p.df <- data.frame(Predicted = p, Value = all_var_data[[var]], 
                         Value_original = original_sequence,
                         Variable = var, SpInd = s)
      out[[s]] <- p.df
    }
  }
  return(out)
}

## Updated interaction function
dynOcc_interaction <- function(var1 = "tmax_anom_jja",
                               var2 = "hsf_pland1_10",
                               model_component = "gamma", 
                               param_compile = param.compile,
                               raw_data = bdata,
                               climate_season = "jja",  # NEW PARAMETER
                               pred_point = 50,
                               var2_levels = 5) {
  
  # Get season-appropriate configuration
  var_config <- get_var_config(climate_season)
  config <- var_config[[model_component]]
  
  if(is.null(config)) stop("Model component not recognized")
  if(!var1 %in% config$vars) stop("var1 not available")
  if(!var2 %in% config$vars) stop("var2 not available")
  
  # Apply scaling function
  apply_scaling <- function(x, scale_type) {
    switch(scale_type,
           "z" = scale(x)[,1],
           "z_no_center" = scale(x, center=FALSE)[,1],
           "raw" = x,
           x)
  }
  
  # Define dynamic variables for this season
  season <- tolower(climate_season)
  dynamic_vars <- c(paste0("tmax_anom_", season), paste0("p_anom_", season), "hsf_pland1_10")
  
  # Process variables
  other_var_means <- list()
  var1_seq_scaled <- NULL
  var2_seq_scaled <- NULL
  var1_seq_orig <- NULL
  var2_seq_orig <- NULL
  
  for(i in 1:length(config$vars)) {
    var_name <- config$vars[i]
    data_col <- config$data_cols[i] 
    scale_type <- config$scaling[i]
    
    # Get raw data
    if(data_col %in% dynamic_vars) {
      if(is.list(raw_data[[data_col]])){
        raw_data[[data_col]] <- unlist(raw_data[[data_col]])
      }
      raw_values <- as.vector(raw_data[[data_col]])
    } else {
      raw_values <- raw_data[[data_col]]
    }
    
    scaled_values <- apply_scaling(raw_values, scale_type)
    
    if(var_name == var1) {
      var_range <- range(scaled_values, na.rm = TRUE)
      var1_seq_scaled <- seq(var_range[1], var_range[2], length.out = pred_point)
      orig_range <- range(raw_values, na.rm = TRUE)
      var1_seq_orig <- seq(orig_range[1], orig_range[2], length.out = pred_point)
      
    } else if(var_name == var2) {
      var_range <- range(scaled_values, na.rm = TRUE)
      var2_seq_scaled <- seq(var_range[1], var_range[2], length.out = var2_levels)
      orig_range <- range(raw_values, na.rm = TRUE)
      var2_seq_orig <- seq(orig_range[1], orig_range[2], length.out = var2_levels)
      
    } else {
      other_var_means[[var_name]] <- switch(scale_type,
                                            "z" = 0,
                                            "z_no_center" = 0,
                                            "raw" = mean(raw_values, na.rm = TRUE),
                                            mean(scaled_values, na.rm = TRUE))
    }
  }
  
  # Get coefficients
  mod.out <- param_compile |> 
    filter(str_detect(variable, model_component)) |> 
    select(sp.ind, variable, mean) |> 
    pivot_wider(names_from = variable, values_from = mean)
  
  nspec <- length(unique(mod.out$sp.ind))
  results_list <- list()
  
  # Get season-specific variable names
  temp_anom <- paste0("tmax_anom_", season)
  precip_anom <- paste0("p_anom_", season)
  temp_trend <- paste0("tmax_trend_", season)
  tmin_trend <- paste0("tmin_trend_", season)
  precip_trend <- paste0("prec_trend_", season)
  
  for(s in 1:nspec) {
    species_data <- mod.out[s, ]
    
    species_results <- expand.grid(
      var1_scaled = var1_seq_scaled,
      var2_scaled = var2_seq_scaled,
      stringsAsFactors = FALSE
    )
    
    species_results$var1_original <- rep(var1_seq_orig, times = var2_levels)
    species_results$var2_original <- rep(var2_seq_orig, each = pred_point)
    
    predictions <- numeric(nrow(species_results))
    
    for(j in 1:nrow(species_results)) {
      all_vars <- other_var_means
      all_vars[[var1]] <- species_results$var1_scaled[j]
      all_vars[[var2]] <- species_results$var2_scaled[j]
      
      # Calculate linear predictor using season-specific variables
      intercept_col <- paste0(model_component, "0")
      coef_prefix <- model_component
      
      lin_pred <- species_data[[intercept_col]] + 
        species_data[[paste0(coef_prefix, "1")]] * all_vars[[temp_anom]] + 
        species_data[[paste0(coef_prefix, "2")]] * all_vars[[precip_anom]] + 
        species_data[[paste0(coef_prefix, "3")]] * all_vars[["hsf_pland1_10"]] + 
        species_data[[paste0(coef_prefix, "4")]] * all_vars[[temp_anom]] * all_vars[["hsf_pland1_10"]] +
        species_data[[paste0(coef_prefix, "5")]] * all_vars[[precip_anom]] * all_vars[["hsf_pland1_10"]] +
        species_data[[paste0(coef_prefix, "6")]] * all_vars[[temp_trend]] +
        species_data[[paste0(coef_prefix, "7")]] * all_vars[[tmin_trend]] +
        species_data[[paste0(coef_prefix, "8")]] * all_vars[[precip_trend]] +
        species_data[[paste0(coef_prefix, "9")]] * all_vars[[temp_trend]] * all_vars[["hsf_pland1_10"]] +
        species_data[[paste0(coef_prefix, "10")]] * all_vars[[tmin_trend]] * all_vars[["hsf_pland1_10"]] +
        species_data[[paste0(coef_prefix, "11")]] * all_vars[[precip_trend]] * all_vars[["hsf_pland1_10"]]
      
      predictions[j] <- plogis(lin_pred)
    }
    
    species_results$Predicted <- predictions
    species_results$SpInd <- s
    species_results$var1_name <- var1
    species_results$var2_name <- var2
    species_results$var2_level <- factor(rep(1:var2_levels, each = pred_point))
    species_results$var2_level_label <- factor(rep(round(var2_seq_orig, 2), each = pred_point))
    
    results_list[[s]] <- species_results
  }
  
  return(results_list)
}

## Updated significance checking function
check_specific_interaction_significance <- function(param_compile, 
                                                    var1, 
                                                    var2, 
                                                    model_component = "gamma",
                                                    climate_season = "jja") {  # NEW PARAMETER
  
  # Extract season from variable names if not provided explicitly
  season <- tolower(climate_season)
  
  # Create season-specific variable names for mapping
  temp_anom <- paste0("tmax_anom_", season)
  precip_anom <- paste0("p_anom_", season)
  temp_trend <- paste0("tmax_trend_", season)
  tmin_trend <- paste0("tmin_trend_", season)
  precip_trend <- paste0("prec_trend_", season)
  
  # Define the mapping from variable pairs to coefficient names
  interaction_mapping <- list(
    gamma = list(),
    eps = list()
  )
  
  # Populate mappings with season-specific variables
  for(component in c("gamma", "eps")) {
    interaction_mapping[[component]][[paste0(temp_anom, ":hsf_pland1_10")]] <- paste0(component, "4")
    interaction_mapping[[component]][[paste0(precip_anom, ":hsf_pland1_10")]] <- paste0(component, "5")
    interaction_mapping[[component]][[paste0(temp_trend, ":hsf_pland1_10")]] <- paste0(component, "9")
    interaction_mapping[[component]][[paste0(tmin_trend, ":hsf_pland1_10")]] <- paste0(component, "10")
    interaction_mapping[[component]][[paste0(precip_trend, ":hsf_pland1_10")]] <- paste0(component, "11")
  }
  
  # Create the key for this variable pair
  var1_tmp <- var1; var2_tmp <- var2
  
  if(stringr::str_detect(var1, "hsf")){
    var1 <- var2_tmp
    var2 <- var1_tmp
  }
  
  var_key <- paste(c(var1, var2), collapse = ":")
  
  # Find the coefficient name for this specific interaction
  coef_name <- interaction_mapping[[model_component]][[var_key]]
  
  if(is.null(coef_name)) {
    warning("No interaction coefficient found for ", var1, " x ", var2, " in ", model_component)
    return(data.frame(sp.ind = unique(param_compile$sp.ind), 
                      has_significant_interaction = FALSE))
  }
  
  # Check significance
  significance_results <- param_compile %>%
    filter(variable == coef_name) %>%
    mutate(
      has_significant_interaction = (lowci > 0 & hici > 0) | (lowci < 0 & hici < 0)
    ) %>%
    select(sp.ind, has_significant_interaction)
  
  return(significance_results)
}

# Update the plotting function to pass through climate_season
me_interaction_plots <- function(dat_list,
                                 plot_type = "lines",
                                 use_original_scale = TRUE,
                                 selected_species = NULL,
                                 facet_species = TRUE,
                                 species_names = NULL,
                                 palette = "viridis",
                                 x_label = NULL,
                                 y_label = "Predicted Probability",
                                 show_significance = TRUE,
                                 param_compile = NULL,
                                 model_component = "gamma",
                                 climate_season = "jja") {  # NEW PARAMETER
  
  # Combine list to dataframe
  dat <- do.call(rbind, dat_list)
  
  # Get the variable names from the data
  var1_name <- unique(dat$var1_name)[1]
  var2_name <- unique(dat$var2_name)[1]
  
  # Add species names if provided
  if(!is.null(species_names)) {
    dat$species_name <- species_names[dat$SpInd]
    facet_var <- "species_name"
  } else {
    facet_var <- "SpInd"
  }
  
  # Check for significance of THIS SPECIFIC interaction
  facet_labels <- NULL
  if(show_significance && !is.null(param_compile)) {
    
    cat("Checking significance for interaction:", var1_name, "x", var2_name, "in season", climate_season, "\n")
    
    significance_data <- check_specific_interaction_significance(
      param_compile, var1_name, var2_name, model_component, climate_season
    )
    
    # Create labels with asterisks for significant interactions
    if(!is.null(species_names)) {
      sig_vector <- significance_data$has_significant_interaction[
        match(unique(dat$SpInd), significance_data$sp.ind)
      ]
      sig_vector[is.na(sig_vector)] <- FALSE
      
      facet_labels <- ifelse(sig_vector, 
                             paste0(species_names[unique(dat$SpInd)], "*"),
                             species_names[unique(dat$SpInd)])
      names(facet_labels) <- species_names[unique(dat$SpInd)]
      
    } else {
      species_indices <- unique(dat$SpInd)
      sig_vector <- significance_data$has_significant_interaction[
        match(species_indices, significance_data$sp.ind)
      ]
      sig_vector[is.na(sig_vector)] <- FALSE
      
      facet_labels <- ifelse(sig_vector,
                             paste0("Species ", species_indices, "*"),
                             paste0("Species ", species_indices))
      names(facet_labels) <- species_indices
    }
    
    cat("Number of species with significant interactions:", sum(sig_vector, na.rm = TRUE), "\n")
  }
  
  # Filter species if specified
  if(!is.null(selected_species)) {
    dat <- dat %>% filter(SpInd %in% selected_species)
    
    if(!is.null(facet_labels)) {
      if(!is.null(species_names)) {
        keep_labels <- species_names[selected_species]
        facet_labels <- facet_labels[names(facet_labels) %in% keep_labels]
      } else {
        facet_labels <- facet_labels[names(facet_labels) %in% as.character(selected_species)]
      }
    }
  }
  
  # Choose scales
  x_var <- if(use_original_scale) "var1_original" else "var1_scaled"
  y_var_for_heatmap <- if(use_original_scale) "var2_original" else "var2_scaled"
  
  # Create the plot (same plotting code as before)
  if(plot_type == "lines") {
    p <- ggplot(dat, aes(x = !!sym(x_var), y = Predicted, 
                         color = var2_level_label)) +
      geom_line(size = 0.8, alpha = 0.8) +
      scale_color_viridis_d(option = palette, name = var2_name) +
      theme_bw() +
      labs(x = x_label %||% var1_name, y = y_label)
    
  } else if(plot_type == "heatmap") {
    p <- ggplot(dat, aes(x = !!sym(x_var), 
                         y = !!sym(y_var_for_heatmap), 
                         fill = Predicted)) +
      geom_tile() +
      scale_fill_viridis_c(option = palette, name = "Predicted\nProbability") +
      theme_bw() +
      labs(x = x_label %||% var1_name, 
           y = var2_name)
    
  } else if(plot_type == "ribbon") {
    dat_summary <- dat %>%
      group_by(!!sym(x_var), SpInd) %>%
      summarise(
        mean_pred = mean(Predicted),
        min_pred = min(Predicted),
        max_pred = max(Predicted),
        .groups = "drop"
      )
    
    p <- ggplot(dat_summary, aes(x = !!sym(x_var))) +
      geom_ribbon(aes(ymin = min_pred, ymax = max_pred), 
                  alpha = 0.3, fill = "steelblue") +
      geom_line(aes(y = mean_pred), color = "steelblue", size = 1) +
      theme_bw() +
      labs(x = x_label %||% var1_name, y = y_label)
  }
  
  # Add faceting with significance asterisks
  if(facet_species && length(unique(dat$SpInd)) > 1) {
    if(show_significance && !is.null(facet_labels)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), 
                          scales = "free_y",
                          labeller = labeller(.default = facet_labels))
    } else {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_y")
    }
  }
  
  return(p)
}

# Keep the existing plotting functions (me_plots) unchanged as they don't need modification
me_plots <- function(dat, 
                     x_var = "Value", 
                     y_var = "Predicted", 
                     color_var = "SpInd",
                     palette = "viridis",
                     group_only = FALSE,
                     line_color = "gray50",
                     line_alpha = 0.6,
                     line_size = 0.5,
                     x_label = NULL,
                     y_label = "Predicted Probability",
                     show_legend = TRUE) {
  
  p <- ggplot(data = dat, aes(x = !!sym(x_var), y = !!sym(y_var))) + 
    theme_bw() + 
    xlab(x_label %||% x_var) +
    ylab(y_label)
  
  if(group_only) {
    p <- p + 
      geom_line(aes(group = !!sym(color_var)), 
                color = line_color, 
                alpha = line_alpha, 
                size = line_size)
    
  } else {
    p <- p + 
      aes(color = factor(!!sym(color_var))) +
      geom_line(alpha = 0.8, size = 0.8) +
      labs(color = "Species")
    
    if(palette %in% c("viridis", "plasma", "inferno", "magma", "cividis", "turbo")) {
      p <- p + scale_color_viridis_d(option = palette)
    } else {
      p <- p + scale_color_viridis_d()
    }
    
    if(!show_legend) {
      p <- p + theme(legend.position = "none")
    }
  }
  
  return(p)
}
