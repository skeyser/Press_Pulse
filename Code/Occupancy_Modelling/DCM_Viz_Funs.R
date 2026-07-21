## ----------------------------------------------------------
##
## Script name: DCM Viz Functions
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-04-06
##
## Email: skeyser@wisc.edu
##
## Github: https://github.com/skeyser
##
## -----------------------------------------------------------
##
## Notes:
##
##
## -----------------------------------------------------------

## Defaults
options(scipen = 6, digits = 4)

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
## What we need is the variables held at their means for the effects and to
## focus on the predictor of interest

dynOcc_margeff <- function(var = "tmax_anom_jja",
                           model_component = "gamma",
                           param_compile = param.compile,
                           raw_data = bdata,
                           level = "species",
                           pred_point = 100){
  
  ## Define variables for each model component
  var_config <- list(
    beta = list(
      vars = c("btmax_annual", "bprec_annual", "cc"),
      data_cols = c("btmax_annual", "bprec_annual", "cc"),
      scaling = c("z", "z", "z"),
      coef_pattern = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6"),  # linear, quad, linear, quad, linear, quad
      has_quadratic = TRUE
    ),
    gamma = list(
      vars = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      data_cols = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      scaling = c("z_no_center", "z_no_center", "z", "raw", "raw", "raw"),
      coef_pattern = c("gamma1", "gamma2", "gamma3", "gamma6", "gamma7", "gamma8"),
      interactions = list(
        list(vars = c("tmax_anom_jja", "hsf_pland1_10"), coef = "gamma4"),
        list(vars = c("p_anom_jja", "hsf_pland1_10"), coef = "gamma5"),
        list(vars = c("tmax_trend_jja", "hsf_pland1_10"), coef = "gamma9"),
        list(vars = c("tmin_trend_jja", "hsf_pland1_10"), coef = "gamma10"),
        list(vars = c("prec_trend_jja", "hsf_pland1_10"), coef = "gamma11")
      ),
      has_quadratic = FALSE
    ),
    eps = list(
      vars = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      data_cols = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      scaling = c("z_no_center", "z_no_center", "z", "raw", "raw", "raw"),
      coef_pattern = c("eps1", "eps2", "eps3", "eps6", "eps7", "eps8"),
      interactions = list(
        list(vars = c("tmax_anom_jja", "hsf_pland1_10"), coef = "eps4"),
        list(vars = c("p_anom_jja", "hsf_pland1_10"), coef = "eps5"),
        list(vars = c("tmax_trend_jja", "hsf_pland1_10"), coef = "eps9"),
        list(vars = c("tmin_trend_jja", "hsf_pland1_10"), coef = "eps10"),
        list(vars = c("prec_trend_jja", "hsf_pland1_10"), coef = "eps11")
      ),
      has_quadratic = FALSE
    )
  )
  
  # Get configuration
  config <- var_config[[model_component]]
  if(is.null(config)) {
    stop("Model component not recognized: ", model_component)
  }
  
  # Check if variable of interest is valid
  if(!var %in% config$vars) {
    stop("Variable '", var, "' not available for component '", model_component, 
         "'. Available: ", paste(config$vars, collapse = ", "))
  }
  
  # Function to apply scaling transformations
  apply_scaling <- function(x, scale_type) {
    switch(scale_type,
           "z" = scale(x)[,1],                      # standardize (mean=0, sd=1)
           "z_no_center" = scale(x, center=FALSE)[,1],  # scale by SD only  
           "raw" = x,                               # no scaling
           x                                        # default
    )
  }
  
  # Process ALL variables - calculate their means and scaling
  all_var_data <- list()
  var_scaling_info <- list()
  original_sequence <- NULL  
  
  for(i in 1:length(config$vars)) {
    var_name <- config$vars[i]
    data_col <- config$data_cols[i] 
    scale_type <- config$scaling[i]
    
    # Get raw data
    if(data_col %in% c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10")) {
      if(is.list(raw_data[[data_col]])){
        raw_data[[data_col]] <- unlist(raw_data[[data_col]])
      }
      # Dynamic variable - flatten matrix across all years
      raw_values <- as.vector(raw_data[[data_col]])
    } else {
      # Static variable
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
      # All other variables - set to their means (which is 0 for scaled vars)
      all_var_data[[var_name]] <- switch(scale_type,
                                         "z" = 0,  # mean of standardized var is 0
                                         "z_no_center" = 0,  # mean when center=FALSE is also 0
                                         "raw" = mean(raw_values, na.rm = TRUE),
                                         mean(scaled_values, na.rm = TRUE)
      )
    }
  }
  
  
  ## Filter the dataframe to the components of interest
  mod.out <- param_compile |> 
    filter(str_detect(variable, model_component)) |> 
    select(sp.ind, variable, mean) |> 
    pivot_wider(names_from = variable, values_from = mean)
  
  # Get number of species
  nspec <- length(unique(mod.out$sp.ind))
  
  ## Create the receiving structure
  out <- vector("list", nspec)
  ## Formula for each
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
    } else if (model_component == "gamma"){
      lin_pred <- mod.out$gamma0[s] + 
        mod.out$gamma1[s] * all_var_data$tmax_anom_jja + 
        mod.out$gamma2[s] * all_var_data$p_anom_jja + 
        mod.out$gamma3[s] * all_var_data$hsf_pland1_10 + 
        mod.out$gamma4[s] * all_var_data$tmax_anom_jja * all_var_data$hsf_pland1_10 +
        mod.out$gamma5[s] * all_var_data$p_anom_jja * all_var_data$hsf_pland1_10 +
        mod.out$gamma6[s] * all_var_data$tmax_trend_jja +
        mod.out$gamma7[s] * all_var_data$tmin_trend_jja +
        mod.out$gamma8[s] * all_var_data$prec_trend_jja +
        mod.out$gamma9[s] * all_var_data$tmax_trend_jja * all_var_data$hsf_pland1_10 +
        mod.out$gamma10[s] * all_var_data$tmin_trend_jja * all_var_data$hsf_pland1_10 +
        mod.out$gamma11[s] * all_var_data$prec_trend_jja * all_var_data$hsf_pland1_10
      p <- plogis(lin_pred)
      p.df <- data.frame(Predicted = p, Value = all_var_data[[var]], 
                         Value_original = original_sequence,
                         Variable = var, SpInd = s)
      out[[s]] <- p.df
    } else if (model_component == "eps"){
      lin_pred <- mod.out$eps0[s] + 
        mod.out$eps1[s] * all_var_data$tmax_anom_jja + 
        mod.out$eps2[s] * all_var_data$p_anom_jja + 
        mod.out$eps3[s] * all_var_data$hsf_pland1_10 + 
        mod.out$eps4[s] * all_var_data$tmax_anom_jja * all_var_data$hsf_pland1_10 +
        mod.out$eps5[s] * all_var_data$p_anom_jja * all_var_data$hsf_pland1_10 +
        mod.out$eps6[s] * all_var_data$tmax_trend_jja +
        mod.out$eps7[s] * all_var_data$tmin_trend_jja +
        mod.out$eps8[s] * all_var_data$prec_trend_jja +
        mod.out$eps9[s] * all_var_data$tmax_trend_jja * all_var_data$hsf_pland1_10 +
        mod.out$eps10[s] * all_var_data$tmin_trend_jja * all_var_data$hsf_pland1_10 +
        mod.out$eps11[s] * all_var_data$prec_trend_jja * all_var_data$hsf_pland1_10
      p <- plogis(lin_pred)
      p.df <- data.frame(Predicted = p, Value = all_var_data[[var]], 
                         Value_original = original_sequence,
                         Variable = var, SpInd = s)
      out[[s]] <- p.df
    }
  }
  return(out)
}

me_plots <- function(dat, 
                     x_var = "Value", 
                     y_var = "Predicted", 
                     color_var = "SpInd",
                     palette = "viridis",
                     group_only = FALSE,  # NEW: Just group, don't color code
                     line_color = "gray50",  # Color for grouped lines
                     line_alpha = 0.6,     # Alpha for grouped lines
                     line_size = 0.5,      # Size for grouped lines
                     x_label = NULL,
                     y_label = "Predicted Probability",
                     show_legend = TRUE) {
  
  p <- ggplot(data = dat, aes(x = !!sym(x_var), y = !!sym(y_var))) + 
    theme_bw() + 
    xlab(x_label %||% x_var) +
    ylab(y_label)
  
  if(group_only) {
    # Group by species but don't color code - all lines same color
    p <- p + 
      geom_line(aes(group = !!sym(color_var)), 
                color = line_color, 
                alpha = line_alpha, 
                size = line_size)
    
  } else {
    # Original behavior - color code by species
    p <- p + 
      aes(color = factor(!!sym(color_var))) +
      geom_line(alpha = 0.8, size = 0.8) +
      labs(color = "Species")
    
    # Add color palette
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

## Interaction plots
dynOcc_interaction <- function(var1 = "tmax_anom_jja",
                               var2 = "hsf_pland1_10",
                               model_component = "gamma", 
                               param_compile = param.compile,
                               raw_data = bdata,
                               pred_point = 50,
                               var2_levels = 5) {
  
  # Your existing var_config...
  var_config <- list(
    gamma = list(
      vars = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      data_cols = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      scaling = c("z_no_center", "z_no_center", "z", "raw", "raw", "raw"),
      has_quadratic = FALSE
    ),
    eps = list(
      vars = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      data_cols = c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10", "tmax_trend_jja", "tmin_trend_jja", "prec_trend_jja"),
      scaling = c("z_no_center", "z_no_center", "z", "raw", "raw", "raw"),
      has_quadratic = FALSE
    )
  )
  
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
  
  # Process variables to get means for "other" variables
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
    if(data_col %in% c("tmax_anom_jja", "p_anom_jja", "hsf_pland1_10")) {
      if(is.list(raw_data[[data_col]])){
        raw_data[[data_col]] <- unlist(raw_data[[data_col]])
      }
      raw_values <- as.vector(raw_data[[data_col]])
    } else {
      raw_values <- raw_data[[data_col]]
    }
    
    scaled_values <- apply_scaling(raw_values, scale_type)
    
    if(var_name == var1) {
      # Variable 1 - create sequence
      var_range <- range(scaled_values, na.rm = TRUE)
      var1_seq_scaled <- seq(var_range[1], var_range[2], length.out = pred_point)
      
      orig_range <- range(raw_values, na.rm = TRUE)
      var1_seq_orig <- seq(orig_range[1], orig_range[2], length.out = pred_point)
      
    } else if(var_name == var2) {
      # Variable 2 - create levels
      var_range <- range(scaled_values, na.rm = TRUE)
      var2_seq_scaled <- seq(var_range[1], var_range[2], length.out = var2_levels)
      
      orig_range <- range(raw_values, na.rm = TRUE)
      var2_seq_orig <- seq(orig_range[1], orig_range[2], length.out = var2_levels)
      
    } else {
      # Other variables - set to means
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
  
  # Create results for each species
  results_list <- list()
  
  for(s in 1:nspec) {
    species_data <- mod.out[s, ]
    
    # Create dataframe for this species
    species_results <- expand.grid(
      var1_scaled = var1_seq_scaled,
      var2_scaled = var2_seq_scaled,
      stringsAsFactors = FALSE
    )
    
    # Add original scale values
    species_results$var1_original <- rep(var1_seq_orig, times = var2_levels)
    species_results$var2_original <- rep(var2_seq_orig, each = pred_point)
    
    # Calculate predictions for each combination
    predictions <- numeric(nrow(species_results))
    
    for(j in 1:nrow(species_results)) {
      # Create complete variable list for this prediction
      all_vars <- other_var_means  # Start with means for other variables
      all_vars[[var1]] <- species_results$var1_scaled[j]  # Set var1
      all_vars[[var2]] <- species_results$var2_scaled[j]  # Set var2
      
      # Calculate linear predictor
      if(model_component == "gamma") {
        lin_pred <- species_data$gamma0 + 
          species_data$gamma1 * all_vars$tmax_anom_jja + 
          species_data$gamma2 * all_vars$p_anom_jja + 
          species_data$gamma3 * all_vars$hsf_pland1_10 + 
          species_data$gamma4 * all_vars$tmax_anom_jja * all_vars$hsf_pland1_10 +
          species_data$gamma5 * all_vars$p_anom_jja * all_vars$hsf_pland1_10 +
          species_data$gamma6 * all_vars$tmax_trend_jja +
          species_data$gamma7 * all_vars$tmin_trend_jja +
          species_data$gamma8 * all_vars$prec_trend_jja +
          species_data$gamma9 * all_vars$tmax_trend_jja * all_vars$hsf_pland1_10 +
          species_data$gamma10 * all_vars$tmin_trend_jja * all_vars$hsf_pland1_10 +
          species_data$gamma11 * all_vars$prec_trend_jja * all_vars$hsf_pland1_10
        
      } else if(model_component == "eps") {
        lin_pred <- species_data$eps0 + 
          species_data$eps1 * all_vars$tmax_anom_jja + 
          species_data$eps2 * all_vars$p_anom_jja + 
          species_data$eps3 * all_vars$hsf_pland1_10 + 
          species_data$eps4 * all_vars$tmax_anom_jja * all_vars$hsf_pland1_10 +
          species_data$eps5 * all_vars$p_anom_jja * all_vars$hsf_pland1_10 +
          species_data$eps6 * all_vars$tmax_trend_jja +
          species_data$eps7 * all_vars$tmin_trend_jja +
          species_data$eps8 * all_vars$prec_trend_jja +
          species_data$eps9 * all_vars$tmax_trend_jja * all_vars$hsf_pland1_10 +
          species_data$eps10 * all_vars$tmin_trend_jja * all_vars$hsf_pland1_10 +
          species_data$eps11 * all_vars$prec_trend_jja * all_vars$hsf_pland1_10
      }
      
      predictions[j] <- plogis(lin_pred)
    }
    
    # Add predictions and metadata
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

me_interaction_plots <- function(dat_list,
                                 plot_type = "lines",
                                 use_original_scale = TRUE,
                                 selected_species = NULL,
                                 facet_species = TRUE,
                                 palette = "viridis",
                                 x_label = NULL,
                                 y_label = "Predicted Probability") {
  
  # Combine list to dataframe
  dat <- do.call(rbind, dat_list)
  
  # Filter species if specified
  if(!is.null(selected_species)) {
    dat <- dat %>% filter(SpInd %in% selected_species)
  }
  
  # Choose x variable scale
  x_var <- if(use_original_scale) "var1_original" else "var1_scaled"
  y_var_for_heatmap <- if(use_original_scale) "var2_original" else "var2_scaled"
  
  if(plot_type == "lines") {
    # Line plot - different lines for var2 levels
    p <- ggplot(dat, aes(x = !!sym(x_var), y = Predicted, 
                         color = var2_level_label)) +
      geom_line(size = 0.8, alpha = 0.8) +
      scale_color_viridis_d(option = palette, name = unique(dat$var2_name)) +
      theme_bw() +
      labs(x = x_label %||% unique(dat$var1_name), y = y_label)
    
  } else if(plot_type == "heatmap") {
    # Heatmap - FIXED
    p <- ggplot(dat, aes(x = !!sym(x_var), 
                         y = !!sym(y_var_for_heatmap), 
                         fill = Predicted)) +
      geom_tile() +
      scale_fill_viridis_c(option = palette, name = "Predicted\nProbability") +
      theme_bw() +
      labs(x = x_label %||% unique(dat$var1_name), 
           y = unique(dat$var2_name))
    
  } else if(plot_type == "ribbon") {
    # Ribbon showing range across var2
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
      labs(x = x_label %||% unique(dat$var1_name), y = y_label)
  }
  
  # Add faceting if requested
  if(facet_species && length(unique(dat$SpInd)) > 1) {
    p <- p + facet_wrap(~SpInd, scales = "free_y")
  }
  
  return(p)
}

# Function to check if interaction coefficients are "significant"
check_interaction_significance <- function(param_compile, model_component = "gamma") {
  
  # Define interaction coefficients for each model component
  interaction_coefs <- list(
    gamma = c("gamma4", "gamma5", "gamma9", "gamma10", "gamma11"),  # All the interaction terms
    eps = c("eps4", "eps5", "eps9", "eps10", "eps11")
  )
  
  coefs_to_check <- interaction_coefs[[model_component]]
  
  # Check significance for each species
  significance_results <- param_compile %>%
    filter(variable %in% coefs_to_check) %>%
    mutate(
      significant = case_when(
        lowci > 0 & hici > 0 ~ "positive",
        lowci < 0 & hici < 0 ~ "negative", 
        TRUE ~ "non_significant"
      ),
      any_significant = significant != "non_significant"
    ) %>%
    group_by(sp.ind) %>%
    summarise(
      has_significant_interaction = any(any_significant),
      n_significant = sum(any_significant),
      .groups = "drop"
    )
  
  return(significance_results)
}

me_interaction_plots <- function(dat_list,
                                 plot_type = "lines",
                                 use_original_scale = TRUE,
                                 selected_species = NULL,
                                 facet_species = TRUE,
                                 species_names = NULL,
                                 palette = "viridis",
                                 x_label = NULL,
                                 y_label = "Predicted Probability",
                                 show_significance = TRUE,  # NEW
                                 param_compile = NULL,      # NEW - for significance testing
                                 model_component = "gamma",  # NEW
                                 significance_alpha = 0.11) { # NEW - significance level (89% CI)
  
  # Combine list to dataframe
  dat <- do.call(rbind, dat_list)
  
  # Add species names if provided
  if(!is.null(species_names)) {
    dat$species_name <- species_names[dat$SpInd]
    facet_var <- "species_name"
    facet_labels <- species_names
  } else {
    facet_var <- "SpInd"
    facet_labels <- unique(dat$SpInd)
  }
  
  # Check for significance if requested
  significance_data <- NULL
  if(show_significance && !is.null(param_compile)) {
    significance_data <- check_interaction_significance(param_compile, model_component)
    
    # Create labels with asterisks
    if(!is.null(species_names)) {
      facet_labels <- ifelse(significance_data$has_significant_interaction[match(1:length(species_names), significance_data$sp.ind)],
                             paste0(species_names, "*"),
                             species_names)
      names(facet_labels) <- species_names
    } else {
      species_indices <- unique(dat$SpInd)
      facet_labels <- ifelse(significance_data$has_significant_interaction[match(species_indices, significance_data$sp.ind)],
                             paste0("Species ", species_indices, "*"),
                             paste0("Species ", species_indices))
      names(facet_labels) <- species_indices
    }
  }
  
  # Filter species if specified
  if(!is.null(selected_species)) {
    dat <- dat %>% filter(SpInd %in% selected_species)
  }
  
  # Choose scales
  x_var <- if(use_original_scale) "var1_original" else "var1_scaled"
  y_var_for_heatmap <- if(use_original_scale) "var2_original" else "var2_scaled"
  
  if(plot_type == "lines") {
    p <- ggplot(dat, aes(x = !!sym(x_var), y = Predicted, 
                         color = var2_level_label)) +
      geom_line(size = 0.8, alpha = 0.8) +
      scale_color_viridis_d(option = palette, name = unique(dat$var2_name)) +
      theme_bw() +
      labs(x = x_label %||% unique(dat$var1_name), y = y_label)
    
  } else if(plot_type == "heatmap") {
    p <- ggplot(dat, aes(x = !!sym(x_var), 
                         y = !!sym(y_var_for_heatmap), 
                         fill = Predicted)) +
      geom_tile() +
      scale_fill_viridis_c(option = palette, name = "Predicted\nProbability") +
      theme_bw() +
      labs(x = x_label %||% unique(dat$var1_name), 
           y = unique(dat$var2_name))
    
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
      labs(x = x_label %||% unique(dat$var1_name), y = y_label)
  }
  
  # Add faceting with significance asterisks
  if(facet_species && length(unique(dat$SpInd)) > 1) {
    if(show_significance && !is.null(significance_data)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), 
                          scales = "free_y",
                          labeller = labeller(.default = facet_labels))
    } else {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_y")
    }
  }
  
  return(p)
}

check_specific_interaction_significance <- function(param_compile, 
                                                    var1, 
                                                    var2, 
                                                    model_component = "gamma") {
  
  # Define the mapping from variable pairs to coefficient names
  interaction_mapping <- list(
    gamma = list(
      "tmax_anom_jja:hsf_pland1_10" = "gamma4",
      "p_anom_jja:hsf_pland1_10" = "gamma5", 
      "tmax_trend_jja:hsf_pland1_10" = "gamma9",
      "tmin_trend_jja:hsf_pland1_10" = "gamma10",
      "prec_trend_jja:hsf_pland1_10" = "gamma11"
    ),
    eps = list(
      "tmax_anom_jja:hsf_pland1_10" = "eps4",
      "p_anom_jja:hsf_pland1_10" = "eps5",
      "tmax_trend_jja:hsf_pland1_10" = "eps9", 
      "tmin_trend_jja:hsf_pland1_10" = "eps10",
      "prec_trend_jja:hsf_pland1_10" = "eps11"
    )
  )
  
  # Create the key for this variable pair (order matters for consistency)
  var_key <- paste(c(var1, var2), collapse = ":")
  
  # Find the coefficient name for this specific interaction
  coef_name <- interaction_mapping[[model_component]][[var_key]]
  
  if(is.null(coef_name)) {
    warning("No interaction coefficient found for ", var1, " x ", var2, " in ", model_component)
    return(data.frame(sp.ind = unique(param_compile$sp.ind), 
                      has_significant_interaction = FALSE))
  }
  
  # Check significance ONLY for this specific interaction coefficient
  significance_results <- param_compile %>%
    filter(variable == coef_name) %>%
    mutate(
      has_significant_interaction = (lowci > 0 & hici > 0) | (lowci < 0 & hici < 0)
    ) %>%
    select(sp.ind, has_significant_interaction)
  
  return(significance_results)
}

me_interaction_plots <- function(dat_list,
                                 plot_type = "lines",
                                 use_original_scale = TRUE,
                                 selected_species = NULL,
                                 facet_species = TRUE,
                                 species_names = NULL,
                                 palette = "viridis",
                                 x_label = NULL,
                                 y_label = "Predicted Probability",
                                 #legend_lab = NULL,
                                 show_significance = TRUE,
                                 param_compile = NULL,
                                 model_component = "gamma") {
  
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
    
    cat("Checking significance for interaction:", var1_name, "x", var2_name, "\n")
    
    significance_data <- check_specific_interaction_significance(
      param_compile, var1_name, var2_name, model_component
    )
    
    # Create labels with asterisks for significant interactions
    if(!is.null(species_names)) {
      # Match species indices to species names
      sig_vector <- significance_data$has_significant_interaction[
        match(unique(dat$SpInd), significance_data$sp.ind)
      ]
      sig_vector[is.na(sig_vector)] <- FALSE  # Handle any missing matches
      
      facet_labels <- ifelse(sig_vector, 
                             paste0(species_names[unique(dat$SpInd)], "*"),
                             species_names[unique(dat$SpInd)])
      names(facet_labels) <- species_names[unique(dat$SpInd)]
      
    } else {
      # Use species indices
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
    
    # Also filter facet labels if they exist
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
  
  # Create the plot
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
