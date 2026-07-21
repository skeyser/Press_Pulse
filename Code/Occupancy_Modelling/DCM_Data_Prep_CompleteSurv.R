## -------------------------------------------------------------
##
## Script name: Multi-season, multi-species Data Prep w/ only
## complete sites
##
## Script purpose: Prepare data for Bayesian multi-season model fit.
##
## Author: Spencer R Keyser
##
## Date Created: 2025-09-16
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
## -------------------------------------------------------------
##
## Notes:
## This script was written originally for prepping one year of
## data for a single season MSOM. I need to update the code to
## handle processing data for a multi-season or stacked static 
## occupancy modelling framework.
## Before I start editing decisions on the modelling packages
## is important.
## 1. spOccupancy (depends on the handling of multi-season; no colex in spOccupancy)
## 2. JAGS/Nimble (current implementation w/ FPs)
## 3. ubms (unlikely)
## 4. unmarked (potential)
## 5. flocker (FP, if I can get it to fit)
##
##
## -------------------------------------------------------------

## Defaults
options(scipen = 6, digits = 4)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(here)
library(abind)
library(lubridate)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Species Data
##
## -------------------------------------------------------------

## Load in the bird data for 2021-2025
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter_Thresh2024/2021_99Conf_OccSppList.RData"))
sp.det.list21 <- sp.det.list
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter_Thresh2024/2022_99Conf_OccSppList.RData"))
sp.det.list22 <- sp.det.list
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter_Thresh2024/2023_99Conf_OccSppList.RData"))
sp.det.list23 <- sp.det.list
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter_Thresh2024/2024_99Conf_OccSppList.RData"))
sp.det.list24 <- sp.det.list
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter_Thresh2024/2025_99Conf_OccSppList.RData"))
sp.det.list25 <- sp.det.list

names(sp.det.list)
str(sp.det.list21)

## Set the sampling interval for binning daily data apriori
bin <- 6

# ## Checking spatial locations
# coords <- st_read(here("./Data/Spatial_Data/ARU_Locs_2021_2025.shp"))
# coords <- coords |> filter(srvy_yr == 2021)
# 
# wbnu <- sp.det.list21$`White-breasted Nuthatch`
# result <- wbnu %>%
#   rowwise() %>%
#   mutate(Max_Occ = {
#     vals <- c_across(-Cell_Unit)
#     max_val <- max(vals, na.rm = TRUE)
#     ifelse(is.infinite(max_val), NA, max_val)
#   }) %>%
#   select(Cell_Unit, Max_Occ) |> 
#   mutate(Cell_Unit = ifelse(stringr::str_detect(Cell_Unit, "C[0-9]{4}_"), Cell_Unit, gsub("C", "C0", Cell_Unit)))
# 
# 
# ## Join by cell_unit
# wbnu <- coords |> 
#   rename("Cell_Unit" = "Cll_Unt") |> 
#   left_join(result)
# 
# ggplot(wbnu) +
#   geom_sf(aes(color = as.factor(Max_Occ)))

## Initial occupancy across species
init_occ <- function(x){
  in_occ <- rowSums(x[,-1], na.rm = T)
  in_occ <- sum(in_occ > 0)/sum(!is.na(in_occ))
  return(in_occ)
}

## Filter species by naive occupancy
sp.list.fun <- function(sp_list, naive_occ = 0.1){
  sp_inocc <- lapply(sp_list, init_occ)
  sp_inocc <- do.call(c, sp_inocc)
  species <- sp_inocc[sp_inocc >= naive_occ]
  species <- names(species)
  species <- species[!str_detect(species, "Kestrel|Hawk|Eagle|Falcon")]
  return(species)
}

## Take a species if naive occ > 10% on any year
species <- unique(c(
  sp.list.fun(sp.det.list21),
  sp.list.fun(sp.det.list22),
  sp.list.fun(sp.det.list23),
  sp.list.fun(sp.det.list24),
  sp.list.fun(sp.det.list25)))

## Extract these species from each list
sp.extract <- function(x, species){
  x <- x[names(x) %in% species]
  return(x)
}

## Apply across all the years
sp.det.list21 <- sp.extract(sp.det.list21, species = species)
sp.det.list22 <- sp.extract(sp.det.list22, species = species)
sp.det.list23 <- sp.extract(sp.det.list23, species = species)
sp.det.list24 <- sp.extract(sp.det.list24, species = species)
sp.det.list25 <- sp.extract(sp.det.list25, species = species)

## Group acoustically indistinguishable species
combine_species <- function(species_list, species_to_combine, new_name, date_cols) {
  # Start with first species
  combined_df <- species_list[[species_to_combine[1]]]
  
  # Sum values across species
  for(sp in species_to_combine[-1]) {
    combined_df[, date_cols] <- combined_df[, date_cols] + species_list[[sp]][, date_cols]
  }
  
  # Binarize: convert all values > 0 to 1
  combined_df <- combined_df %>%
    mutate(across(all_of(date_cols), ~as.integer(. > 0)))
  
  # Remove original species and add combined
  species_list[species_to_combine] <- NULL
  species_list[[new_name]] <- combined_df
  
  return(species_list)
}

## Combine Sapsuckers
sp.det.list21 <- combine_species(species_list = sp.det.list21, 
                                 species_to_combine = c("Red-breasted Sapsucker", "Red-naped Sapsucker", "Williamson's Sapsucker"), 
                                 new_name = "Sphyrapicus spp.",
                                 date_cols = colnames(sp.det.list21[[1]])[str_detect(colnames(sp.det.list21[[1]]), "Cell_Unit", negate = T)])
sp.det.list22 <- combine_species(species_list = sp.det.list22, 
                                 species_to_combine = c("Red-breasted Sapsucker", "Red-naped Sapsucker", "Williamson's Sapsucker"), 
                                 new_name = "Sphyrapicus spp.",
                                 date_cols = colnames(sp.det.list22[[1]])[str_detect(colnames(sp.det.list22[[1]]), "Cell_Unit", negate = T)])
sp.det.list23 <- combine_species(species_list = sp.det.list23, 
                                 species_to_combine = c("Red-breasted Sapsucker", "Red-naped Sapsucker", "Williamson's Sapsucker"), 
                                 new_name = "Sphyrapicus spp.",
                                 date_cols = colnames(sp.det.list23[[1]])[str_detect(colnames(sp.det.list23[[1]]), "Cell_Unit", negate = T)])
sp.det.list24 <- combine_species(species_list = sp.det.list24, 
                                 species_to_combine = c("Red-breasted Sapsucker", "Red-naped Sapsucker", "Williamson's Sapsucker"), 
                                 new_name = "Sphyrapicus spp.",
                                 date_cols = colnames(sp.det.list24[[1]])[str_detect(colnames(sp.det.list24[[1]]), "Cell_Unit", negate = T)])
sp.det.list25 <- combine_species(species_list = sp.det.list25, 
                                 species_to_combine = c("Red-breasted Sapsucker", "Red-naped Sapsucker", "Williamson's Sapsucker"), 
                                 new_name = "Sphyrapicus spp.",
                                 date_cols = colnames(sp.det.list25[[1]])[str_detect(colnames(sp.det.list25[[1]]), "Cell_Unit", negate = T)])

names(sp.det.list21)
names(sp.det.list25)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: 3-D array for site x rep x year
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Put this data into an array, but first we need to align the data
## so that we have an idea of the matching sites
## Add the padded 0 to ensure all cell units have C + 4 digits
padZero <- function(x){
  x$Cell_Unit <- ifelse(stringr::str_detect(x$Cell_Unit, "C[0-9]{4}_"), x$Cell_Unit, gsub("C", "C0", x$Cell_Unit))
  x <- x[order(x$Cell_Unit), ]
  return(x)
}

## Apply the zero padding across all data sets for alignment
sp.det.list21 <- lapply(sp.det.list21, padZero)
sp.det.list22 <- lapply(sp.det.list22, padZero)
sp.det.list23 <- lapply(sp.det.list23, padZero)
sp.det.list24 <- lapply(sp.det.list24, padZero)
sp.det.list25 <- lapply(sp.det.list25, padZero)

## Find all of the unique units across the 4 years
findUniq <- function(x){
  uniqUnits <- unique(x$Cell_Unit)
  return(uniqUnits)
}

## Get unique units for each species-year combination
uniq_units21 <- unique(unlist(lapply(sp.det.list21, findUniq)))
uniq_units22 <- unique(unlist(lapply(sp.det.list22, findUniq)))
uniq_units23 <- unique(unlist(lapply(sp.det.list23, findUniq)))
uniq_units24 <- unique(unlist(lapply(sp.det.list24, findUniq)))
uniq_units25 <- unique(unlist(lapply(sp.det.list25, findUniq)))

## Combine all unique units and get unique values
uniq_units <- unique(c(uniq_units21, uniq_units22, uniq_units23, uniq_units24, uniq_units25))

## These are the total unique units across the 5 survey years: 2095 ARUs
message(paste("Number of unique ARUs from 2021-2025:", length(uniq_units)))
head(uniq_units)
tail(uniq_units)

## We need to ensure that the number of unique units standardizes array for the occ model
template <- data.frame(Cell_Unit = uniq_units)

## Fill the missing ARUs with the complete sites
sp.det.fill <- function(x) {
  sp.det.fill <- template |> 
    left_join(x) |> 
    mutate(across(where(is.logical), as.numeric))
  return(sp.det.fill)
}

## Apply filling template across all species for each year
sp.det.std21 <- lapply(sp.det.list21, sp.det.fill)
sp.det.std22 <- lapply(sp.det.list22, sp.det.fill)
sp.det.std23 <- lapply(sp.det.list23, sp.det.fill)
sp.det.std24 <- lapply(sp.det.list24, sp.det.fill)
sp.det.std25 <- lapply(sp.det.list25, sp.det.fill)

## Double-check alignment with the rows 
nrow(sp.det.std21[[1]])
nrow(sp.det.std21[[1]]) == nrow(sp.det.std21[[2]])
nrow(sp.det.std21[[1]]) == nrow(sp.det.std22[[1]])
nrow(sp.det.std21[[1]]) == nrow(sp.det.std25[[1]])

## -------------------------------------------------------------
##
## End Section: 
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Format the detection covariates as site x rep x yr
##
## -------------------------------------------------------------

## Here is where we are making the decision to collapse the secondary 
## sampling units into 6-day non-overlapping bins.
## We do the same for the effort data later on
## Array for the data
## D1 (i) = Site, D2 (j) = Sampling Date, D3 (k) = species
samp.cols <- colnames(sp.det.list21[[1]])[str_detect(colnames(sp.det.list21[[1]]), "\\d")]
samp.cols <- as.Date(samp.cols, format = "%Y-%m-%d")

## Create different sampling periods based on bin
second_samp <- function(DAT, interval, id_col, eff = F, e.var = "Days"){
  tmp.cols <- colnames(DAT)[str_detect(colnames(DAT), "\\d")]
  tmp.splits <- split(tmp.cols, ceiling(seq_along(tmp.cols) / interval))
  tmp.new <- as.data.frame(cbind(
    DAT[, id_col], 
    matrix(data = NA, 
           nrow = nrow(DAT), 
           ncol = length(tmp.splits), 
           dimnames = list(NULL, paste0("J", seq(1:length(tmp.splits)))))
  ))
  for(i in 1:length(tmp.splits)){
    if(!eff){
      Jsum <- rowSums(DAT[, tmp.splits[[i]]], na.rm = T)
      Jsum <- ifelse(Jsum > 0, 1, 0)
    } 
    if(eff & e.var == "Days"){
      Jsum <- ifelse(DAT[, tmp.splits[[i]]] > 0, 1, 0)
      Jsum <- rowSums(Jsum, na.rm = T)
    }
    if(eff & e.var == "Hrs"){
      Jsum <- rowSums(DAT[, tmp.splits[[i]]], na.rm = T)
    }
    if(eff & e.var == "FirstJDay"){
      mjd <- as.Date(gsub("[[:punct:]]", "-", tmp.splits[[i]]), format = "%Y-%m-%d")
      mjd <- median(lubridate::yday(mjd))
      Jsum <- mjd
    }
    
    tmp.new[,i+1] <- Jsum
  }
  return(tmp.new)
}

## Apply the function
sp.det.list21.r <- lapply(sp.det.std21, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = F))
sp.det.list21.r <- lapply(sp.det.list21.r, function(x) x |> dplyr::select(matches("J[0-9]")))
sp.det.list22.r <- lapply(sp.det.std22, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = F))
sp.det.list22.r <- lapply(sp.det.list22.r, function(x) x |> dplyr::select(matches("J[0-9]")))
sp.det.list23.r <- lapply(sp.det.std23, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = F))
sp.det.list23.r <- lapply(sp.det.list23.r, function(x) x |> dplyr::select(matches("J[0-9]")))
sp.det.list24.r <- lapply(sp.det.std24, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = F))
sp.det.list24.r <- lapply(sp.det.list24.r, function(x) x |> dplyr::select(matches("J[0-9]")))
sp.det.list25.r <- lapply(sp.det.std25, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = F))
sp.det.list25.r <- lapply(sp.det.list25.r, function(x) x |> dplyr::select(matches("J[0-9]")))

## Now we want to merge the files across the years
## We need to check to ensure the matrices are square
sp.det.std <- vector(mode = "list", length = length(sp.det.list21.r))
for(i in 1:length(sp.det.std)){
  
  sp.det.list.tmp <- unlist(list(sp.det.list21.r[i], 
                                 sp.det.list22.r[i], 
                                 sp.det.list23.r[i],
                                 sp.det.list24.r[i],
                                 sp.det.list25.r[i]),
                            recursive = FALSE)
  
  years <- c(2021:2025)
  
  # Ensure all arrays have the same column names
  all_cols <- unique(unlist(lapply(sp.det.list.tmp, colnames)))
  sp.det.list.tmp <- lapply(sp.det.list.tmp, function(x) {
    missing_cols <- setdiff(all_cols, colnames(x))
    if(length(missing_cols) > 0) {
      x[missing_cols] <- NA
    }
    x[, all_cols]  # Reorder columns consistently
  })
  
  sp.det.std[[i]] <- abind(sp.det.list.tmp, along = 3)
  dimnames(sp.det.std[[i]]) <- list(NULL,
                                    colnames(sp.det.list.tmp[[1]]),
                                    paste(names(sp.det.list.tmp), years))
}

str(sp.det.std)
str(sp.det.std[[1]])
str(sp.det.std[[2]])

## Package the data up as a 4-D array for the y
dates <- dimnames(sp.det.std[[1]])[[2]]
years <- str_extract(dimnames(sp.det.std[[1]])[[3]], "\\d{4}")
species <- unique(unlist(lapply(sp.det.std, function(x) {
  trimws(str_extract(dimnames(x)[[3]], "[^0-9]+"))
})))

# Create 4-D array
sp.det.array <- abind(sp.det.std, along = 4)
str(sp.det.array)
# Add dimension names
dimnames(sp.det.array) <- list(
  Cell_Unit = uniq_units,
  Date = dates,
  Year = years,
  Species = species
)

str(sp.det.array)

y4d <- sp.det.array

## Make this into an array
nsite <- dim(y4d)[1]
nsurv <- dim(y4d)[2]
nyear <- dim(y4d)[3]
nspec <- dim(y4d)[4]

## Get the number of hours per survey for the detection covariate
## This needs to be edited to handle all years of data and 
## backfill in the sites that are missing with NA for NA covs
eff.dat <- list.files(here("Data/Occ_Data/Thresh_By_Species_NoDetFilter_Thresh2024/"), pattern = "Effort", full.names = T)

eff.dat <- lapply(eff.dat, read.csv)

eff.dat <- lapply(eff.dat, function(x) x[,-1])

eff.dat <- lapply(eff.dat, function(x){
  x <- x |> mutate(Cell_Unit = ifelse(stringr::str_detect(Cell_Unit, "C[0-9]{4}_"), 
                                      Cell_Unit, 
                                      gsub("C", "C0", Cell_Unit)))
  return(x)                 
})

eff.dat <- lapply(eff.dat, function(x){
  colnames(x) <- gsub("[[:punct:]]", "_", gsub("X", "", colnames(x)))
  return(x)
})

eff.days <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = T, e.var = "Days"))
eff.days <- lapply(eff.days, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  x <- x |> arrange(Cell_Unit)
  return(x)
})

eff.hrs <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = T, e.var = "Hrs"))
eff.hrs <- lapply(eff.hrs, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  x <- x |> arrange(Cell_Unit)
  return(x)
})

eff.jday <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = bin, id_col = "Cell_Unit", eff = T, e.var = "FirstJDay"))
eff.jday <- lapply(eff.jday, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  x <- x |> arrange(Cell_Unit)
  return(x)
})

## Now we need to bundle this up for all of the unique sites that exist
eff.days.std <- lapply(eff.days, function(x) {
  
  eff.fill <- template |> 
    left_join(x) |>
    mutate(across(where(is.logical), as.numeric))
  eff.fill[is.na(eff.fill)] <- 0
  return(eff.fill)
})

eff.hrs.std <- lapply(eff.hrs, function(x) {
  
  eff.fill <- template |> 
    left_join(x) |> 
    mutate(across(where(is.logical), as.numeric))
  eff.fill[is.na(eff.fill)] <- 0
  return(eff.fill)
})

eff.jday.std <- lapply(eff.jday, function(x) {
  
  eff.fill <- template |> 
    left_join(x) |> 
    mutate(across(where(is.logical), as.numeric))
  eff.fill[is.na(eff.fill)] <- 0
  return(eff.fill)
})

eff.days.std <- lapply(eff.days.std, function(x) x |> dplyr::select(matches("J[0-9]")))
eff.days.std <- array(unlist(eff.days.std), dim = c(nsite, nsurv, nyear),
                      dimnames = list(uniq_units,
                                      colnames(eff.days.std[[1]]),
                                      c(2021:2025)))

eff.hrs.std <- lapply(eff.hrs.std, function(x) x |> dplyr::select(matches("J[0-9]")))
eff.hrs.std <- array(unlist(eff.hrs.std), dim = c(nsite, nsurv, nyear),
                     dimnames = list(uniq_units,
                                     colnames(eff.hrs.std[[1]]),
                                     c(2021:2025)))

eff.jday.std <- lapply(eff.jday.std, function(x) x |> dplyr::select(matches("J[0-9]")))
eff.jday.std <- array(unlist(eff.jday.std), dim = c(nsite, nsurv, nyear),
                      dimnames = list(uniq_units,
                                      colnames(eff.jday.std[[1]]),
                                      c(2021:2025)))

## Mask the 0s to NAs since these values correspond to missing biological data
eff.hrs.std.na <- eff.hrs.std
eff.hrs.std.na[eff.hrs.std.na == 0] <- NA 
mean.hrs <- mean(eff.hrs.std.na, na.rm = T)
sd.hrs <- sd(eff.hrs.std.na, na.rm = T)

## Some statistics on the number of recording hrs
cat("Mean Hours:", round(mean.hrs, 2),
    "\nSD Hours:", round(sd.hrs, 2),
    "\nRange Hours:", round(range(eff.hrs.std.na, na.rm = T)))

## Scale by mean and sd
eff.hrs.std.sc <- (eff.hrs.std - mean.hrs) / sd.hrs

eff.jday.std.na <- eff.jday.std
eff.jday.std.na[eff.jday.std.na == 0] <- NA 
mean.jday <- mean(eff.jday.std.na, na.rm = T)
sd.jday <- sd(eff.jday.std.na, na.rm = T)

## JDAY info
cat("Mean Jday:", round(mean.jday, 2),
    "\nSD Jday:", round(sd.jday, 2),
    "\nRange Jday:", round(range(eff.jday.std.na, na.rm = T)))


## Scale by mean and sd
eff.jday.std.sc <- (eff.jday.std - mean.jday) / sd.jday


## Set y to NA for days without sampling
# Iterate through the matrix and set corresponding values in the 3D array to NA
# Iterate through the matrix and set corresponding values in the 3D array to NA for all layers
y4d.new <- y4d
for (i in 1:dim(eff.days.std)[1]) {
  for (j in 1:dim(eff.days.std)[2]) {
    for(t in 1:dim(eff.days.std)[3]) {
      for(k in 1:nspec){
        if (as.numeric(eff.days.std[i, j, t]) == 0) {
          # Set all layers (third dimension) of the array at position (i, j) to NA
          y4d.new[i, j, t, k] <- NA
        }
      }
    }
  }
}
## Print the modified 3D array
## This looks good...effort data matches biological data
head(y4d.new[,,1,])
head(y4d.new[,,2,])
head(y4d.new[1:5,,1,1])
head(eff.hrs.std[1:5,,1])

# Keep only Hermit warbler: dims become site x date x year (3D array)
ht <- y4d.new[, , , "Hermit Warbler"]
print(ht, max = 1000)

non_na_counts <- apply(ht, c(1, 3), function(x) sum(!is.na(x)))
keep_sites <- apply(non_na_counts > 0, 1, all)
filtered_ht <- ht[keep_sites, , ]


surveys <- non_na_counts # site x year matrix of survey counts

# Boolean: was site sampled in year?
sampled <- surveys > 0

summarise_dynamic_coverage <- function(sampled) {
  ## total sites
  n_total <- nrow(sampled)
  
  ## year names
  years <- colnames(sampled)
  
  ## ---- Single-year coverage ----
  single_year <- colSums(sampled)
  
  ## ---- Year-pair coverage ----
  year_pairs <- c()
  for (i in 1:(length(years) - 1)) {
    pairname <- paste0(years[i], "-", years[i+1])
    year_pairs[pairname] <- sum(sampled[, years[i]] & sampled[, years[i+1]])
  }
  
  ## ---- Multi-year consecutive coverage ----
  multi_sets_counts <- c(
    all_years = sum(apply(sampled, 1, all)),
    first3    = sum(sampled[, years[1]] & sampled[, years[2]] & sampled[, years[3]]),
    first4 = sum(sampled[, years[1]] & sampled[, years[2]] & sampled[, years[3]] & sampled[, years[4]]),
    last3     = sum(sampled[, years[2]] & sampled[, years[3]] & sampled[, years[4]]),
    last4 = sum(sampled[, years[2]] & sampled[, years[3]] & sampled[, years[4]] & sampled[, years[5]])
  )
  
  # Give readable names for multi-year sequences
  names(multi_sets_counts) <- c(
    paste(years, collapse = "-"),
    paste(years[1:3], collapse = "-"),
    paste(years[1:4], collapse = "-"),
    paste(years[2:4], collapse = "-"),
    paste(years[2:5], collapse = "-")
  )
  
  ## ---- Combine and calculate percentages ----
  df_single <- data.frame(
    Sequence = names(single_year),
    n_sites = as.integer(single_year),
    pct_sites = round(single_year / n_total * 100, 1)
  )
  
  df_pairs <- data.frame(
    Sequence = names(year_pairs),
    n_sites = as.integer(year_pairs),
    pct_sites = round(year_pairs / n_total * 100, 1)
  )
  
  df_multi <- data.frame(
    Sequence = names(multi_sets_counts),
    n_sites = as.integer(multi_sets_counts),
    pct_sites = round(multi_sets_counts / n_total * 100, 1)
  )
  
  coverage_summary <- rbind(df_single, df_pairs, df_multi)
  rownames(coverage_summary) <- NULL
  
  return(list(
    n_total_sites = n_total,
    coverage_summary = coverage_summary
  ))
}

result <- summarise_dynamic_coverage(sampled)

cat("Total sites:", result$n_total_sites, "\n")
print(result$coverage_summary)

## ***********************************************************
##
## Section Notes: For more robust estimation and faster
## model fitting we could consider doing the sites that
## have 5 year sample coverage 
##
## ***********************************************************
## Get the site names for sites that have complete survey coverage
aru.comp <- rownames(sampled[rowSums(sampled) == ncol(sampled),])



## -------------------------------------------------------------
##
## End Section:
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Wrangling covariations for DCM
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fire Covariates - Dynamic only
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ***********************************************************
##
## Section Notes: We have a series of decisions to make with
## the fire data.
## 1. Use continuous metrics for severity?
##  a. Proportion burned at low/mod vs high severity
##  b. Mean compositve burn index
## 2. Use dummy variables
##  a. >50% burned HSF or > 2.3 mean CBI
## 3. Add in more nuanced variables?
##  a. Time-since fire?
##  b. Fire return interval?
##
## Regardless of how this is done, ideally we can limit
## the amount of arbitrary binning that needs to happen.
## There are a million ways to do this. 
## Current approach uses pland of HSF 1-10 years ago, but this
## neglects low/mod severity fire which might be important for
## a some species (e.g., aerial insectivores, ground nesters, etc.)
##
## ***********************************************************

## Load the fire data
fire1_10 <- readRDS(here("Data/FireMets_ARU_21_25_AllUnitsByYears_1_10yr.RDS"))
fire_tsf <- readRDS(here("Data/FireMets_ARU_21_25_AllUnitsByYears_TSF.RDS"))
fire_tsf_hsf <- readRDS(here("Data/FireMets_ARU_21_25_AllUnitsByYears_TSF_HSFonly.RDS"))
fire_pland_severity_all <- readRDS(here("Data/FireMets_ARU_21_25_AllUnitsByYears_PLAND_FireSev_1985_2024.RDS")) 

## Build a full fire metrics dataset from different partitions
str(fire1_10)
str(fire_tsf)
str(fire_tsf_hsf)
str(fire_pland_severity_all)

## Add unique character strings for organization
## We want a wide structure going forward
fire1_10$FireSeverity <- fire1_10$FireSeverity |> 
  rename(Survey_Year = Fire_Sev_SurvExtYr) |> 
  mutate(Survey_Year = as.character(Survey_Year)) |> 
  rename_with(~ paste0("1-10_", .x, recycle0 = TRUE),
              starts_with("Fire"))

fire_pland_severity_all$FireSeverity <- fire_pland_severity_all$FireSeverity |> 
  rename(Survey_Year = Fire_Sev_SurvExtYr) |>
  mutate(Survey_Year = as.character(Survey_Year)) |>
  rename_with(~ paste0("1-40_", .x, recycle0 = TRUE),
              starts_with("Fire"))

fire_pland_severity_all$FireLscp <- fire_pland_severity_all$FireLscp |> 
  rename(Survey_Year = Year) |>
  mutate(Survey_Year = as.character(Survey_Year)) |>
  mutate(across(contains("pland"), ~ifelse(is.na(.x), 0, .x)))

fire1_10$FireLscp <- fire1_10$FireLscp |>
  rename(Survey_Year = Year) |>
  mutate(Survey_Year = as.character(Survey_Year)) |>
  mutate(across(contains("pland"), ~ifelse(is.na(.x), 0, .x))) |> 
  select(deployment_name, Survey_Year, contains("pland"))

fire_tsf <- fire_tsf |> 
  rename_with(~ paste0(.x, "_AllFire", recycle0 = TRUE),
              starts_with("time")) |> 
  tidyr::pivot_longer(cols = contains("time_since_fire"),
                      values_to = "TSF_AllFire",
                      names_to = "Survey_Year") |> 
  mutate(Survey_Year = str_extract(Survey_Year, "\\d+"))

fire_tsf_hsf <- fire_tsf_hsf |> 
  rename_with(~ paste0(.x, "_HSF", recycle0 = TRUE),
              starts_with("time")) |> 
  tidyr::pivot_longer(cols = contains("time_since_fire"),
                      values_to = "TSF_HSF",
                      names_to = "Survey_Year") |> 
  mutate(Survey_Year = str_extract(Survey_Year, "\\d+"))

## Merge everything into one dataframe for all fire data
fire_full <- fire1_10$FireSeverity |> 
  full_join(fire_pland_severity_all$FireSeverity) |>
  full_join(fire1_10$FireLscp) |>
  full_join(fire_pland_severity_all$FireLscp) |>
  full_join(fire_tsf) |> 
  full_join(fire_tsf_hsf) |> 
  select(deployment_name, Ref_Survey_Yr = Survey_Year, everything()) |>
  mutate(Cell_Unit = gsub("G\\d{3}_V\\d{1}_", "", deployment_name)) |> 
  group_by(Cell_Unit, Ref_Survey_Yr) |>
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) |> 
  pivot_wider(
    names_from = Ref_Survey_Yr,
    values_from = -c(Cell_Unit, Ref_Survey_Yr),  # Everything except these two
    names_glue = "{.name}"
  ) |> 
  tibble::column_to_rownames("Cell_Unit")

cor(fire_full) |> corrplot::corrplot()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Climate covariates - Dynamic and static
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load the climate data
clim <- readRDS(here("Data/Climate_DCM_Covs_2020_2025.RDS"))
str(clim)

static_clim <- clim[[1]] |> filter(Cell_Unit %in% uniq_units)

## Quick check on correlations amongst static climate vars
c.cm <- cor(static_clim[,2:21])
corrplot::corrplot(c.cm, method = "number")


dyn_clim <- lapply(clim[[2]], function(x) x[rownames(x) %in% uniq_units,])

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Habitat + Topography - Static only
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Forest data
cfo <- read.csv(here("Data/CFO_ARU_21_25.csv"))
cfo <- cfo |> 
  rename("deployment_name" = "dplymn_") |> 
  rename("Cell_Unit" = "Cll_Unt") |>
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  select(Cell_Unit, CanopyBaseHeight:nlcd) |> 
  group_by(Cell_Unit) |> 
  summarise(across(CanopyBaseHeight:nlcd, mean)) |> 
  tibble::column_to_rownames("Cell_Unit")

## Elevation data
topo <- read.csv(here("Data/Topo_Data_ARU_21_25.csv"))
topo <- topo |> 
  rename("deployment_name" = "dplymn_") |> 
  rename("Cell_Unit" = "Cll_Unt") |>
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  select(Cell_Unit, elevation) |> 
  group_by(Cell_Unit) |> 
  summarise(Elevation = mean(elevation)) |> 
  tibble::column_to_rownames("Cell_Unit") |> 
  as.matrix()

## -------------------------------------------------------------
##
## End Section:
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Subsetting based on complete survey history
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cfo.comp <- cfo[rownames(cfo) %in% aru.comp,]
hist(cfo.comp$nlcd)

topo.comp <- topo[rownames(topo) %in% aru.comp,]
hist(topo.comp)

sclim.comp <- static_clim[static_clim$Cell_Unit %in% aru.comp,]
hist(sclim.comp$Prcp_Base_Annual)
hist(sclim.comp$Tmax_Base_Annual)
hist(sclim.comp$Tmin_Base_Annual)
hist(sclim.comp$Trend_Tmax_MAM)
hist(sclim.comp$Trend_Tmin_MAM)
hist(sclim.comp$Trend_Prcp_MAM)
hist(sclim.comp$Trend_Tmax_JJA)
hist(sclim.comp$Trend_Tmin_JJA)
hist(sclim.comp$Trend_Prcp_JJA)

fire.comp <- fire_full[rownames(fire_full) %in% aru.comp,]
hist(fire.comp$`1-10_Fire_Sev_mean_2021`)
hist(fire.comp$`1-10_Fire_Sev_mean_2022`)
hist(fire.comp$`1-10_Fire_Sev_mean_2023`)
hist(fire.comp$`1-10_Fire_Sev_mean_2024`)
hist(fire.comp$`1-10_Fire_Sev_mean_2025`)
hist(fire.comp$TSF_HSF_2021)
hist(fire.comp$TSF_HSF_2022)
hist(fire.comp$TSF_HSF_2023)
hist(fire.comp$TSF_HSF_2024)
hist(fire.comp$TSF_HSF_2025)

fc <- cbind(sclim.comp, fire.comp)
fc |> 
  select(Trend_Tmax_JJA, contains("1-10_Fire_Sev_mean")) |> 
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "Year",
                      values_to = "Fire_Sev") |> 
  mutate(Year = gsub("1-10_Fire_Sev_mean_", "", Year)) |> 
  ggplot(aes(x = Fire_Sev, y = Trend_Tmax_JJA)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  facet_wrap(~Year) +
  labs(title = "Predictor Coverage Across Years")

fc |> 
  select(Trend_Prcp_JJA, contains("1-10_Fire_Sev_mean")) |> 
  tidyr::pivot_longer(cols = 2:6,
                      names_to = "Year",
                      values_to = "Fire_Sev") |> 
  mutate(Year = gsub("1-10_Fire_Sev_mean_", "", Year)) |> 
  ggplot(aes(x = Fire_Sev, y = Trend_Prcp_JJA)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  facet_wrap(~Year) +
  labs(title = "Predictor Coverage Across Years")

## -------------------------------------------------------------
##
## Begin Section: Wrap data up for NIMBLE
##
## -------------------------------------------------------------
# validSurv <- which(eff.hrs.std > 0, arr.ind = T)
# nValid <- nrow(validSurv)
# 
# win.data <- list(
#   y = y4d.new,
#   btemp = static_clim$Tmax_Base_JJA,
#   bprec = static_clim$Prcp_Base_JJA,
#   trendt = static_clim$Trend_Tmax,
#   trendp = static_clim$Trend_Prcp,
#   Lat = static_clim$Lat,
#   Long = static_clim$Long,
#   cc = cfo$CanopyCover,
#   fire_init = fire_init_dum,
#   tanom = dyn_clim$Tmax,
#   panom = dyn_clim$Prcp,
#   fire = fire_dyn_dum,
#   fire_pland = flscp_1_10_pland[,2:5],
#   eff.hrs.sc = eff.hrs.std.sc,
#   eff.jday.sc = eff.jday.std.sc,
#   eff.hrs = eff.hrs.std,
#   eff.jday = eff.jday.std,
#   nsites = nsite,
#   nyears = nyear,
#   nreps = nsurv,
#   nspec = nspec,
#   nValid = nValid,
#   valid_i = validSurv[,1],
#   valid_j = validSurv[,2],
#   valid_t = validSurv[,3]
# )
# 
# str(win.data)
# 
# saveRDS(win.data, file = here("Data/Occ_Data/MSOM_Multi_Wide_Preds.RDS"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Ragged form bc nimble doesn't like dynamic indexing
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Start by subsetting the ymat
y4d.sub <- y4d.new[dimnames(y4d.new)[[1]] %in% aru.comp,,,]
dim(y4d.sub)
hist(apply(apply(y4d.sub, c(1,4), function(x) max(x, na.rm = T)), 1, sum))
hist(apply(apply(y4d.sub, c(1,4), function(x) max(x, na.rm = T)), 2, sum))

## Subset the effort data
eff.hrs.std.sub <- eff.hrs.std[dimnames(eff.hrs.std)[[1]] %in% aru.comp,,]
dim(eff.hrs.std.sub)

eff.jday.std.sub <- eff.jday.std[dimnames(eff.jday.std)[[1]] %in% aru.comp,,]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Temporary removal of year 5
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# eff.hrs.std <- eff.hrs.std[,,1:4]
# eff.jday.std <- eff.jday.std[,,1:4]
# y4d.new <- y4d.new[,,1:4,]
# nyear <- 4

## -----------------------------------------------------------
##
## End Section: Temporary removal of year 5
##
## -----------------------------------------------------------

valid_idx <- which(eff.hrs.std.sub > 0, arr.ind = T)
valid_idx <- valid_idx[order(valid_idx[,1], valid_idx[,3], valid_idx[,2]), ]
max(valid_idx[,1])
length(unique(valid_idx[,1]))
#valid_idx <- valid_idx[order(valid_idx[,1], valid_idx[,2], valid_idx[,3]), ]

nObs <- nrow(valid_idx)
nsite <- length(unique(rownames(valid_idx)))
validARUs <- unique(rownames(valid_idx))
length(validARUs)
## Get the indexing
obs_site <- as.integer(valid_idx[,1])
obs_rep <- as.integer(valid_idx[,2])
obs_year <- as.integer(valid_idx[,3])

site_ids <- sort(unique(obs_site))
site_map <- setNames(seq_along(site_ids), site_ids)
obs_site <- unname(site_map[as.character(obs_site)])
max(obs_site)
nsite == length(unique(obs_site))

## Change the detection data into long format
y_long <- matrix(NA, nObs, nspec)
eff.hrs.long <- numeric(nObs)
eff.jday.long <- numeric(nObs)
y_wide <- y4d.sub
y_wide <- y_wide[dimnames(y_wide)[[1]] %in% validARUs,,,]
eff.hrs.std.sub <- eff.hrs.std.sub[dimnames(eff.hrs.std.sub)[[1]] %in% validARUs,,]
eff.jday.std.sub <- eff.jday.std.sub[dimnames(eff.jday.std.sub)[[1]] %in% validARUs,,]

for(v in 1:nObs){
  i <- obs_site[v]
  j <- obs_rep[v]
  t <- obs_year[v]
  y_long[v, ] <- y_wide[i,j,t,]
  
  if(eff.hrs.std.sub[i,j,t] == 0){stop("Zero Hit. Check indexing.")}
  else {eff.hrs.long[v] <- eff.hrs.std.sub[i,j,t]}
  
  if(eff.jday.std.sub[i,j,t] == 0){stop("Zero Hit. Check indexing.")}
  else {eff.jday.long[v] <- eff.jday.std.sub[i,j,t]}
}

any(is.na(y_long))
any(is.na(eff.hrs.long))
sum(eff.hrs.long == 0)/length(eff.hrs.long)
any(is.na(eff.jday.long))
sum(eff.jday.long == 0)/length(eff.jday.long)

glimpse(y_long)

## Need to match this for dnyamic predictors too?
static_clim <- static_clim |> filter(Cell_Unit %in% validARUs)
dyn_clim <- lapply(dyn_clim, function(x) x[rownames(x) %in% validARUs, ])
cancov <- cfo$CanopyCover[rownames(cfo) %in% validARUs]
canht <- cfo$CanopyHeight[rownames(cfo) %in% validARUs]
lf <- cfo$LadderFuelDensity[rownames(cfo) %in% validARUs]
nlcd.cc <- cfo$nlcd[rownames(cfo) %in% validARUs]
topo <- topo[rownames(topo) %in% validARUs,]
fire_full <- fire_full[rownames(fire_full) %in% validARUs, ]

win.rag <- list(
  y = y_long,
  y_wide = y_wide,
  btmax_jja = static_clim$Tmax_Base_JJA,
  btmin_jja = static_clim$Tmin_Base_JJA,
  btmax_annual = static_clim$Tmax_Base_Annual,
  btmin_annual = static_clim$Tmin_Base_Annual,
  bprec_jja = static_clim$Prcp_Base_JJA,
  bprec_annual = static_clim$Prcp_Base_Annual,
  tmax_trend_jja = static_clim$Trend_Tmax_JJA,
  tmin_trend_jja = static_clim$Trend_Tmin_JJA,
  prec_trend_jja = static_clim$Trend_Prcp_JJA,
  tmax_trend_mam = static_clim$Trend_Tmax_MAM,
  tmin_trend_mam = static_clim$Trend_Tmin_MAM,
  prec_trend_mam = static_clim$Trend_Prcp_MAM,
  Lat = static_clim$Lat,
  Long = static_clim$Long,
  cc = as.numeric(cancov),
  ch = as.numeric(canht),
  cc.nldc = as.numeric(nlcd.cc),
  lf = as.numeric(lf),
  ele = as.numeric(topo),
  tmax_anom_mam = dyn_clim$MAM_Tmax,
  p_anom_mam = dyn_clim$MAM_Prcp,
  tmax_anom_jja = dyn_clim$JJA_Tmax,
  p_anom_jja = dyn_clim$JJA_Prcp,
  tsf_hsf = fire_full |> select(contains("TSF_HSF")),
  tsf = fire_full |> select(contains("TSF_AllFire")),
  cbi1_10 = fire_full |> select(contains("1-10_Fire_Sev_mean")),
  cbi1_40 = fire_full |> select(contains("1-40_Fire_Sev_mean")),
  hsf_pland1_10 = fire_full |>  select(contains("High_Sev_1-10_pland")),
  hsf_pland1_40 = fire_full |>  select(contains("High_Sev_1-40_pland")),
  eff.hrs = eff.hrs.long,
  eff.jday = eff.jday.long,
  nsites = nsite,
  nyears = nyear,
  nreps = nsurv,
  nspec = nspec,
  nObs = nObs,
  site_obs = obs_site,
  year_obs = obs_year,
  rep_obs = obs_rep,
  cell_id = as.integer(factor(gsub("_U\\d{1}", "", validARUs))),
  n_cells = max(as.integer(factor(gsub("_U\\d{1}", "", validARUs))))
)


str(win.rag)

saveRDS(win.rag, file = here("Data/Occ_Data/DCM_Ragged_Complete_Survey_Only.RDS"))

## -----------------------------------------------------------
##
## Begin Section: Old code
##
## -----------------------------------------------------------

## Checking data
win.rag <- readRDS(file = here("Data/Occ_Data/DCM_Ragged_Full_2021_2025_FireCont_Thresh24.RDS"))

nspec <- win.rag$nspec

zst <- apply(win.rag$y_wide, c(1,3,4), max, na.rm = T)
zst[is.infinite(zst)] <- NA
zst[is.na(zst)] <- 0

zst[1,2,]
win.rag$y[win.rag$site_obs[1], win.rag$year_obs[2]]

errors <- 0
for(v in 1:win.rag$nObs){
  i <- win.rag$site_obs[v]
  j <- win.rag$rep_obs[v]
  t <- win.rag$year_obs[v]
  for(k in 1:nspec) {
    if(win.rag$y[v, k] != win.rag$y_wide[i, j, t, k]) {
      errors <- errors + 1
      cat("Mismatch at v =", v, "site", i, "rep", j, "year", t, "\n")
    }
  }
}
cat(errors, "mismatches found\n")

# Initialise zst with NA
zst <- array(NA, dim = c(win.rag$nsites, win.rag$nyears, win.rag$nspec))
dim(zst)
max(win.rag$site_obs)
max(win.rag$site_obs)

# Loop over ragged rows and mark detections
for(v in 1:win.rag$nObs){
  i <- win.rag$site_obs[v]
  t <- win.rag$year_obs[v]
  for(k in 1:win.rag$nspec){
    y_val <- win.rag$y[v, k]
    if(!is.na(y_val) && y_val == 1){
      zst[i, t, k] <- 1
    } else if(!is.na(y_val) && y_val == 0) {
      # Keep as 0 if not already 1
      if(is.na(zst[i, t, k])) zst[i, t, k] <- 0
    }
  }
}

# Fill remaining NA (unsurveyed site-years) with 0
zst[is.na(zst)] <- 0

check_z_alignment <- function(zst, y_long, site_obs, year_obs, nspec) {
  mismatches <- 0L
  mismatch_list <- list()
  
  for (v in seq_len(nrow(y_long))) {
    i <- site_obs[v]
    t <- year_obs[v]
    
    for (k in seq_len(nspec)) {
      y_val <- y_long[v, k]
      z_val <- zst[i, t, k]
      
      # Check surveyed occasions only
      if (!is.na(y_val)) {
        # Main consistency rule:
        # z must be 1 if detection occurred, can be 0 if no detection
        if (y_val == 1 && z_val != 1) {
          mismatches <- mismatches + 1L
          mismatch_list[[length(mismatch_list) + 1L]] <- 
            sprintf("Mismatch at v=%d (site=%d, year=%d, species=%d): y=1, z=%d", v, i, t, k, z_val)
        }
      }
    }
  }
  
  cat("Total mismatches:", mismatches, "\n")
  if (mismatches > 0) {
    cat("First few mismatches:\n")
    print(head(mismatch_list, 10))
  } else {
    cat("All surveyed occasions align with zst.\n")
  }
  
  invisible(mismatch_list)
}

# Assuming win.rag as your list with ragged data:
check_z_alignment(
  zst        = zst,
  y_long     = win.rag$y,
  site_obs   = win.rag$site_obs,
  year_obs   = win.rag$year_obs,
  nspec      = win.rag$nspec
)

## Checking raw relationships
# Get dims for clarity
dimnames(win.rag$y_wide)  # check ordering: site, rep, year, species
nsites   <- dim(win.rag$y_wide)[1]
nreps    <- dim(win.rag$y_wide)[2]
nyears   <- dim(win.rag$y_wide)[3]
nspecies <- dim(win.rag$y_wide)[4]

## Get the position of the NAs to check MNAR
yy <- win.rag$y_wide[,,,1]
yy.sum <- apply(yy, MARGIN = c(1, 3), FUN = function(x) sum(!is.na(x)))

par(mfrow = c(2,3))

hist(yy.sum[,1])
hist(yy.sum[,2])
hist(yy.sum[,3])
hist(yy.sum[,4])
hist(yy.sum[,5])

# Year index for "year 1" in your occupancy model (often t=1)
year1_index <- 1

# Logical detection/non-detection per site/species in year 1
raw_det_matrix <- matrix(NA, nrow = nsites, ncol = nspecies)

for (sp in 1:nspecies) {
  for (site in 1:nsites) {
    # collapse across all reps for year 1, species sp
    this_year_reps <- win.rag$y_wide[site, , year1_index, sp]
    raw_det_matrix[site, sp] <- as.integer(any(this_year_reps == 1, na.rm = TRUE))
  }
}

colnames(raw_det_matrix) <- paste("Species", dimnames(win.rag$y_wide)[[4]])

raw_occ_df <- data.frame(
  site = 1:nsites,
  cc_raw = win.rag$cc,
  cc   = poly(win.rag$cc, 2)[,1],
  cc2 = poly(win.rag$cc, 2)[,2],
  cc.nlcd = win.rag$cc.nldc,
  ch = win.rag$ch,
  lf = win.rag$lf,
  temp = win.rag$btmax_jja,
  tmin = win.rag$btmin_jja,
  prec = win.rag$bprec_jja,
  Lat = win.rag$Lat,
  Ele = win.rag$ele,
  TSF_HSF = win.rag$tsf_hsf[,1],
  TSF = win.rag$tsf[,1]
)

# Add columns for each species
raw_occ_df <- cbind(raw_occ_df, raw_det_matrix)
raw_occ_long <- raw_occ_df %>%
  pivot_longer(
    cols = starts_with("Species"),
    names_to = "species",
    values_to = "detected"
  ) |> 
  mutate(species = gsub("Species ", "", species))

## Plots
ggplot(raw_occ_long, aes(x = cc_raw, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", 
              method.args = list(family = binomial), 
              formula = y ~ poly(x, 2),
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Canopy Cover",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = cc.nlcd, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", 
              method.args = list(family = binomial), 
              #formula = y ~ I(x^2),
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Canopy Cover: NLCD",
       y = "Naïve Occupancy (≥1 detection in Year 1)")


ggplot(raw_occ_long, aes(x = cc_raw, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Canopy Cover",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = cc_raw, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), 
              formula = y ~ x + I(x^2),
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Canopy Cover",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = cc2, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Canopy Cover",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = scale(ch), y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Canopy Height",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = scale(lf), y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Ladder Fuel Density",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = Ele, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Elevation",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = Ele^2, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Elevation Sq.",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = Ele, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", 
              method.args = list(family = binomial), 
              formula = y ~ poly(x, 2),
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Elevation Sq.",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = temp, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Temperature JJA",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = temp, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", 
              method.args = list(family = binomial), 
              formula = y ~ poly(x, 2),
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Temperature JJA",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = tmin, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Min Temperature JJA",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = tmin, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", 
              method.args = list(family = binomial), 
              formula = y ~ poly(x, 2),
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Min Temperature Sq. JJA",
       y = "Naïve Occupancy (≥1 detection in Year 1)")


ggplot(raw_occ_long, aes(x = temp^2, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Temperature Sq. JJA",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = prec, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Precip JJA",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = prec, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", 
              method.args = list(family = binomial), 
              formula = y ~ poly(x, 2),
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Precip JJA",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = atemp, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Temperature Annnual",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = atemp^2, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Temperature Sq. Annua;",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = atmin, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Min Temperature Annnual",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = aprec, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Precip Annual",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = TrendT, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Temp Trend",
       y = "Naïve Occupancy (≥1 detection in Year 1)")


ggplot(raw_occ_long, aes(x = TrendTmin, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Min Temp Trend",
       y = "Naïve Occupancy (≥1 detection in Year 1)")


ggplot(raw_occ_long, aes(x = TrendP, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Prec Trend",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = TSF_HSF, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Time Since HSF",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = TSF_HSF, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), 
              formula = y ~ poly(x, 2), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Time Since HSF",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = TSF, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), 
              #formula = y ~ poly(x, 2), 
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Time Since Fire",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long |> filter(Sev > 0), aes(x = Sev, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), 
              #formula = y ~ poly(x, 2), 
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Fire Severity 1-10 years",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long, aes(x = pland_HMSF, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), 
              #formula = y ~ poly(x, 2), 
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Prop. High/Mod SF 1-10 years",
       y = "Naïve Occupancy (≥1 detection in Year 1)")

ggplot(raw_occ_long |> filter(Sev > 0), aes(x = Sev, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), 
              formula = y ~ poly(x, 2), 
              se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Fire Severity 1-40 years",
       y = "Naïve Occupancy (≥1 detection in Year 1)")