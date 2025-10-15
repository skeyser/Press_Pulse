## -------------------------------------------------------------
##
## Script name: Multi-season, multi-species Data Prep
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
## 2. JAGS/Nimble (potential)
## 3. ubms (maybe for initial)
## 4. unmarked (unlikely)
## 5. flocker (FP, if I can get it to fit)
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
library(abind)
library(lubridate)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Species Data
##
## -------------------------------------------------------------

## Load in the bird data for 2021

## We start by using 1 species of interest (Lazuli's bunting)
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter/2021_99Conf_OccSppList.RData"))
sp.det.list21 <- sp.det.list
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter/2022_99Conf_OccSppList.RData"))
sp.det.list22 <- sp.det.list
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter/2023_99Conf_OccSppList.RData"))
sp.det.list23 <- sp.det.list
load(here("./Data/Occ_Data/Thresh_By_Species_NoDetFilter/2024_99Conf_OccSppList.RData"))
sp.det.list24 <- sp.det.list

## Species to use
species <- c("Lazuli Bunting",
             "Hermit Warbler",
             "Western Tanager",
             "Western Bluebird",
             "Hermit Thrush",
             "Mountain Chickadee",
             "MacGillivray's Warbler")

sp.extract <- function(x, species){
  x <- x[names(x) %in% species]
  return(x)
}

sp.det.list21 <- sp.extract(sp.det.list21, species = species)
sp.det.list22 <- sp.extract(sp.det.list22, species = species)
sp.det.list23 <- sp.extract(sp.det.list23, species = species)
sp.det.list24 <- sp.extract(sp.det.list24, species = species)

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

sp.det.list21 <- lapply(sp.det.list21, padZero)
sp.det.list22 <- lapply(sp.det.list22, padZero)
sp.det.list23 <- lapply(sp.det.list23, padZero)
sp.det.list24 <- lapply(sp.det.list24, padZero)

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

## Combine all unique units and get unique values
uniq_units <- unique(c(uniq_units21, uniq_units22, uniq_units23, uniq_units24))

## These are the total unique units across the 4 survey years: 2095 ARUs
length(uniq_units)
head(uniq_units)

## We need to ensure that the number of unique units standardizes array for the occ model
template <- data.frame(Cell_Unit = uniq_units)

## Fill the missing ARUs with the complete sites
sp.det.fill <- function(x) {
  
  sp.det.fill <- template |> 
    left_join(x) |> 
    #select(-Cell_Unit) |> 
    mutate(across(where(is.logical), as.numeric))
  return(sp.det.fill)
}

## Apply across all species for each year
sp.det.std21 <- lapply(sp.det.list21, sp.det.fill)
sp.det.std22 <- lapply(sp.det.list22, sp.det.fill)
sp.det.std23 <- lapply(sp.det.list23, sp.det.fill)
sp.det.std24 <- lapply(sp.det.list24, sp.det.fill)

## Double-check
nrow(sp.det.std21[[1]]) == nrow(sp.det.std21[[2]])
nrow(sp.det.std21[[1]]) == nrow(sp.det.std22[[1]])

# ## Now we want to merge the files across the years
# sp.det.std <- vector(mode = "list", length = length(sp.det.std21))
# for(i in 1:length(sp.det.std21)){
#   
#   sp.det.list.tmp <- unlist(list(sp.det.std21[i], 
#                           sp.det.std22[i], 
#                           sp.det.std23[i],
#                           sp.det.std24[i]),
#                           recursive = FALSE)
#   
#   for(j in 1:length(sp.det.list.tmp)){
#     cnames <- colnames(sp.det.list.tmp[[j]])
#     cnames <- paste0(names(sp.det.list.tmp[j]), " ", unique(str_extract(cnames[str_detect(cnames, "\\d{4}")], "\\d{4}")))
#     names(sp.det.list.tmp)[j] <- cnames
#   }
#   
#   sp.det.list.tmp <- lapply(sp.det.list.tmp, function(x) {
#     colnames(x) <- gsub("\\d{4}-", "", colnames(x))
#     return(x)
#   }) 
#   
#   # Ensure all arrays have the same column names
#   all_cols <- unique(unlist(lapply(sp.det.list.tmp, colnames)))
#   sp.det.list.tmp <- lapply(sp.det.list.tmp, function(x) {
#     missing_cols <- setdiff(all_cols, colnames(x))
#     if(length(missing_cols) > 0) {
#       x[missing_cols] <- NA
#     }
#     x[, all_cols]  # Reorder columns consistently
#   })
#   
#   ## Store the Cell_Unit
#   Cell_Unit <- sp.det.list.tmp[[1]]$Cell_Unit
#   
#   sp.det.list.tmp <- lapply(sp.det.list.tmp, function(x) x |> select(-Cell_Unit))
#   
#   sp.det.std[[i]] <- abind(sp.det.list.tmp, along = 3)
#   dimnames(sp.det.std[[i]]) <- list(Cell_Unit,
#                                     colnames(sp.det.list.tmp[[1]]),
#                                     names(sp.det.list.tmp))
# }
# 
# str(sp.det.std)
# 
# 
# ## Package the data up as a 4-D array for the y
# cell_units <- dimnames(sp.det.std[[1]])[[1]]
# dates <- dimnames(sp.det.std[[1]])[[2]]
# years <- str_extract(dimnames(sp.det.std[[1]])[[3]], "\\d{4}")
# species <- unique(unlist(lapply(sp.det.std, function(x) {
#   trimws(str_extract(dimnames(x)[[3]], "[^0-9]+"))
# })))
# 
# # Create 4-D array
# sp.det.array <- abind(sp.det.std, along = 4)
# 
# # Add dimension names
# dimnames(sp.det.array) <- list(
#   Cell_Unit = cell_units,
#   Date = dates,
#   Year = years,
#   Species = species
# )
# 
# str(sp.det.array)
# 
# y4d <- sp.det.array



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

## Array for the data
## D1 (i) = Site, D2 (j) = Sampling Date, D3 (k) = species
samp.cols <- colnames(sp.det.list21[[1]])[str_detect(colnames(sp.det.list21[[1]]), "\\d")]
samp.cols <- as.Date(samp.cols, format = "%Y-%m-%d")

## Create different sampling periods
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
sp.det.list21.r <- lapply(sp.det.std21, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = F))
sp.det.list21.r <- lapply(sp.det.list21.r, function(x) x |> dplyr::select(matches("J[0-9]")))
sp.det.list22.r <- lapply(sp.det.std22, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = F))
sp.det.list22.r <- lapply(sp.det.list22.r, function(x) x |> dplyr::select(matches("J[0-9]")))
sp.det.list23.r <- lapply(sp.det.std23, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = F))
sp.det.list23.r <- lapply(sp.det.list23.r, function(x) x |> dplyr::select(matches("J[0-9]")))
sp.det.list24.r <- lapply(sp.det.std24, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = F))
sp.det.list24.r <- lapply(sp.det.list24.r, function(x) x |> dplyr::select(matches("J[0-9]")))

## I need to insert the code above here
## Now we want to merge the files across the years
sp.det.std <- vector(mode = "list", length = length(sp.det.list21.r))
for(i in 1:length(sp.det.std)){
  
  sp.det.list.tmp <- unlist(list(sp.det.list21.r[i], 
                                 sp.det.list22.r[i], 
                                 sp.det.list23.r[i],
                                 sp.det.list24.r[i]),
                            recursive = FALSE)
  
  years <- c(2021:2024)
  
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


## Package the data up as a 4-D array for the y
dates <- dimnames(sp.det.std[[1]])[[2]]
years <- str_extract(dimnames(sp.det.std[[1]])[[3]], "\\d{4}")
species <- unique(unlist(lapply(sp.det.std, function(x) {
  trimws(str_extract(dimnames(x)[[3]], "[^0-9]+"))
})))

# Create 4-D array
sp.det.array <- abind(sp.det.std, along = 4)

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

## Let's get a sense for the number of surveys that have been repeated
# table(nyear <- apply(y, c(1,2), function(x) sum(!is.na(x))))
# nvisits <- apply(y, c(1,3), function(x) sum(!is.na(x)))
# nvisits[nvisits > 0] <- 1
# table(rowSums(nvisits))
# table(nsites <- apply(y, c(2,3), function(x) sum(!is.na(x))))

## Get the number of hours per survey for the detection covariate
## This needs to be edited to handle all years of data and 
## backfill in the sites that are missing with NA for NA covs
eff.dat <- list.files(here("Data/Occ_Data/Thresh_By_Species_NoDetFilter/"), pattern = "Effort", full.names = T)

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

eff.days <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Days"))
eff.days <- lapply(eff.days, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  x <- x |> arrange(Cell_Unit)
  #x <- x[,-1]
  #colnames(x) <- NULL
  return(x)
})

eff.hrs <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Hrs"))
eff.hrs <- lapply(eff.hrs, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  x <- x |> arrange(Cell_Unit)
  #x <- x[,-1]
  #colnames(x) <- NULL
  return(x)
})

eff.jday <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "FirstJDay"))
eff.jday <- lapply(eff.jday, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  x <- x |> arrange(Cell_Unit)
  #x <- x[,-1]
  #colnames(x) <- NULL
  return(x)
})

## Now we need to bundle this up for all of the unique sites that exist
eff.days.std <- lapply(eff.days, function(x) {
  
  eff.fill <- template |> 
    left_join(x) |> 
    #select(-Cell_Unit) |> 
    mutate(across(where(is.logical), as.numeric))
  eff.fill[is.na(eff.fill)] <- 0
  return(eff.fill)
})

eff.hrs.std <- lapply(eff.hrs, function(x) {
  
  eff.fill <- template |> 
    left_join(x) |> 
    #select(-Cell_Unit) |> 
    mutate(across(where(is.logical), as.numeric))
  eff.fill[is.na(eff.fill)] <- 0
  return(eff.fill)
})

eff.jday.std <- lapply(eff.jday, function(x) {
  
  eff.fill <- template |> 
    left_join(x) |> 
    #select(-Cell_Unit) |> 
    mutate(across(where(is.logical), as.numeric))
  eff.fill[is.na(eff.fill)] <- 0
  return(eff.fill)
})

eff.days.std <- lapply(eff.days.std, function(x) x |> dplyr::select(matches("J[0-9]")))
eff.days.std <- array(unlist(eff.days.std), dim = c(nsite, nsurv, nyear),
                      dimnames = list(uniq_units,
                                      colnames(eff.days.std[[1]]),
                                      c(2021:2024)))

eff.hrs.std <- lapply(eff.hrs.std, function(x) x |> dplyr::select(matches("J[0-9]")))
eff.hrs.std <- array(unlist(eff.hrs.std), dim = c(nsite, nsurv, nyear),
                     dimnames = list(uniq_units,
                                     colnames(eff.hrs.std[[1]]),
                                     c(2021:2024)))

eff.jday.std <- lapply(eff.jday.std, function(x) x |> dplyr::select(matches("J[0-9]")))
eff.jday.std <- array(unlist(eff.jday.std), dim = c(nsite, nsurv, nyear),
                      dimnames = list(uniq_units,
                                      colnames(eff.jday.std[[1]]),
                                      c(2021:2024)))

## Mask the 0s to NAs since these values correspond to missing biological data
eff.hrs.std.na <- eff.hrs.std
eff.hrs.std.na[eff.hrs.std.na == 0] <- NA 
mean.hrs <- mean(eff.hrs.std.na, na.rm = T)
sd.hrs <- sd(eff.hrs.std.na, na.rm = T)

## Scale by mean and sd
eff.hrs.std.sc <- (eff.hrs.std - mean.hrs) / sd.hrs

eff.jday.std.na <- eff.jday.std
eff.jday.std.na[eff.jday.std.na == 0] <- NA 
mean.jday <- mean(eff.jday.std.na, na.rm = T)
sd.jday <- sd(eff.jday.std.na, na.rm = T)

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
# Print the modified 3D array
print(y4d.new[,,1,])
print(y4d.new[,,2,])

# Keep only Hermit Thrush: dims become site x date x year (3D array)
ht <- y4d.new[, , , "Hermit Warbler"]

print(ht, max = 1000)

non_na_counts <- apply(ht, c(1, 3), function(x) sum(!is.na(x)))
keep_sites <- apply(non_na_counts > 0, 1, all)
filtered_ht <- ht[keep_sites, , ]


surveys <- non_na_counts # site x year matrix of survey counts

# Boolean: was site sampled in year?
sampled <- surveys > 0

summarise_dynamic_coverage <- function(sampled) {
  # total sites
  n_total <- nrow(sampled)
  
  # year names
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
    last3     = sum(sampled[, years[2]] & sampled[, years[3]] & sampled[, years[4]])
  )
  
  # Give readable names for multi-year sequences
  names(multi_sets_counts) <- c(
    paste(years, collapse = "-"),
    paste(years[1:3], collapse = "-"),
    paste(years[2:4], collapse = "-")
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

## -------------------------------------------------------------
##
## End Section:
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Format occupancy, persistence, and colonization
## following site x year
##
## -------------------------------------------------------------

## Load the fire data
fire <- readRDS(here("Data/Fire_ARU_21_24_AllUnitsByYears.RDS"))
fsev <- fire

## Combining output before reclassifying
fsev <- fsev |>
  mutate(Cell_Unit = gsub("G\\d{3}_V\\d{1}_", "", deployment_name)) |> 
  group_by(Cell_Unit, Fire_Sev_SurvExtYr) |>
  select(Cell_Unit, Fire_Sev_SurvExtYr, Fire_Sev_mean_1_10) |> 
  summarise(Fire_Sev_mean_1_10 = mean(Fire_Sev_mean_1_10)) |> 
  arrange(Cell_Unit) |>
  mutate(FCat1_10 = case_when(Fire_Sev_mean_1_10 == 0 ~ "Unburned",
                              Fire_Sev_mean_1_10 > 0 & Fire_Sev_mean_1_10 < 2.25 ~ "Low/Mod",
                              Fire_Sev_mean_1_10 >= 2.25 ~ "High"))



## Load the climate data
clim <- readRDS(here("Data/Climate_DCM_Covs.RDS"))
str(clim)

static_clim <- clim[[1]] |> filter(Cell_Unit %in% uniq_units)
dyn_clim <- lapply(clim[[2]], function(x) x[rownames(x) %in% uniq_units,])

# ## Merge the climate and fire data
# fclim <- clim |> 
#   rename("deployment_name" = "dplymn_") |> 
#   #left_join(fsev, by = c("deployment_name", "srvy_yr" = "Fire_Sev_SurvExtYr")) |> 
#   rename("Cell_Unit" = "Cll_Unt") |> 
#   arrange(Cell_Unit) |> 
#   #mutate(LagYr = srvy_yr - Year) |>
#   #filter(LagYr == 1) |> 
#   filter(Year %in% c(2020:2024)) |> 
#   select(Cell_Unit, srvy_yr, Year, bclim, Anomaly, Lat, Long)

## Let's make this into an array to match the data
# fclim.arr <- fclim |> 
#   filter(Cell_Unit %in% uniq_units) |>
#   filter(Year != 2020) |> 
#   arrange(Cell_Unit) |> 
#   tidyr::pivot_wider(names_from = Year,
#                      values_from = Anomaly,
#                      id_cols = "Cell_Unit",
#                      values_fn = mean) |>
#   tibble::column_to_rownames("Cell_Unit") |> 
#   as.matrix()

## Array for the fire data
## I should remedy some of the issue by creating a unique set of 
##ARUs and extracting only from those for each year
## Duplications are a result of small shifts in ARU placement
fire.arr <- fsev |>
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  tidyr::pivot_wider(names_from = Fire_Sev_SurvExtYr,
                     values_from = FCat1_10,
                     id_cols = "Cell_Unit") |>
  tibble::column_to_rownames("Cell_Unit") |>
  as.matrix()

fire_init <- fire.arr[,1]
fire_init <- as.factor(fire_init)
levels(fire_init) <- c("Unburned", "Low/Mod", "High")
str(fire_init)

fire_init_dum <- model.matrix(~fire_init)

# Create dynamic dummy variables across years
# First, create an array to store the dummy variables
# Create dynamic dummy variables maintaining yearly structure
fire_dyn <- fire.arr[,2:4]  # Keep as matrix
fire_dyn_dum <- array(NA, dim=c(nrow(fire_dyn), ncol(fire_dyn), 2))  # sites x years x (categories-1)

# Create dummy variables for each year
for(t in 1:ncol(fire_dyn)) {
  temp_dum <- model.matrix(~factor(fire_dyn[,t], levels=c("Unburned", "Low/Mod", "High")))[,-1]
  fire_dyn_dum[,t,] <- temp_dum
}

str(fire_dyn_dum)

## Take the baseline climate for initial occupancy
# bclim <- fclim |> 
#   filter(Cell_Unit %in% uniq_units) |> 
#   select(Cell_Unit, Year, bclim, Lat, Long) |> 
#   distinct() |> 
#   tidyr::pivot_wider(names_from = Year,
#                      values_from = bclim,
#                      id_cols = "Cell_Unit",
#                      values_fn = mean) |>
#   tibble::column_to_rownames("Cell_Unit") |> 
#   as.matrix()

## Forest data
cfo <- read.csv(here("Data/CFO_ARU_21_24.csv"))
cfo <- cfo |> 
  rename("deployment_name" = "dplymn_") |> 
  rename("Cell_Unit" = "Cll_Unt") |>
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  select(Cell_Unit, CanopyBaseHeight:SurfaceFuels) |> 
  group_by(Cell_Unit) |> 
  summarise(across(CanopyBaseHeight:SurfaceFuels, mean)) |> 
  tibble::column_to_rownames("Cell_Unit")

## Elevation data
topo <- read.csv(here("Data/Topo_Data_ARU_21_24.csv"))
topo <- topo |> 
  rename("deployment_name" = "dplymn_") |> 
  rename("Cell_Unit" = "Cll_Unt") |>
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  select(Cell_Unit, Elevation) |> 
  group_by(Cell_Unit) |> 
  summarise(Elevation = mean(Elevation)) |> 
  tibble::column_to_rownames("Cell_Unit") |> 
  as.matrix()

## -------------------------------------------------------------
##
## End Section:
##
## -------------------------------------------------------------


## -------------------------------------------------------------
##
## Begin Section: Wrap data up for NIMBLE
##
## -------------------------------------------------------------
validSurv <- which(eff.hrs.std > 0, arr.ind = T)
nValid <- nrow(validSurv)

win.data <- list(
  y = y4d.new,
  btemp = static_clim$Tmax_Base_JJA,
  bprec = static_clim$Prcp_Base_JJA,
  trendt = static_clim$Trend_Tmax,
  trendp = static_clim$Trend_Prcp,
  Lat = static_clim$Lat,
  Long = static_clim$Long,
  cc = cancov,
  fire_init = fire_init_dum,
  tanom = dyn_clim$Tmax,
  panom = dyn_clim$Prcp,
  fire = fire_dyn_dum,
  eff.hrs.sc = eff.hrs.std.sc,
  eff.jday.sc = eff.jday.std.sc,
  eff.hrs = eff.hrs.std,
  eff.jday = eff.jday.std,
  nsites = nsite,
  nyears = nyear,
  nreps = nsurv,
  nspec = nspec,
  nValid = nValid,
  valid_i = validSurv[,1],
  valid_j = validSurv[,2],
  valid_t = validSurv[,3]
)

str(win.data)

saveRDS(win.data, file = here("Data/Occ_Data/MSOM_Multi_Test_Categorical_Skip.RDS"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Ragged form bc nimble doesn't like dynamic indexing
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
valid_idx <- which(eff.hrs.std > 0, arr.ind = T)
valid_idx <- valid_idx[order(valid_idx[,1], valid_idx[,3], valid_idx[,2]), ]
max(valid_idx[,1])
length(unique(valid_idx[,1]))
#valid_idx <- valid_idx[order(valid_idx[,1], valid_idx[,2], valid_idx[,3]), ]

nObs <- nrow(valid_idx)
nsite <- length(unique(rownames(valid_idx)))
validARUs <- unique(rownames(valid_idx))

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
y_wide <- y4d.new
y_wide <- y_wide[dimnames(y_wide)[[1]] %in% validARUs,,,]

for(v in 1:nObs){
  i <- obs_site[v]
  j <- obs_rep[v]
  t <- obs_year[v]
  y_long[v, ] <- y_wide[i,j,t,]
  eff.hrs.long[v] <- eff.hrs.std[i,j,t]
  eff.jday.long[v] <- eff.jday.std[i,j,t]
}

any(is.na(y_long))
any(is.na(eff.hrs.long))
any(is.na(eff.jday.long))

glimpse(y_long)

## Need to match this for dnyamic predictors too?
static_clim <- static_clim |> filter(Cell_Unit %in% validARUs)
dyn_clim <- lapply(dyn_clim, function(x) x[rownames(x) %in% validARUs, ])
cancov <- cfo$CanopyCover[rownames(cfo) %in% validARUs]
canht <- cfo$CanopyHeight[rownames(cfo) %in% validARUs]
lf <- cfo$LadderFuelDensity[rownames(cfo) %in% validARUs]
topo <- topo[rownames(topo) %in% validARUs,]
fire_dyn <- fire.arr[,2:4]  # Keep as matrix
fire_dyn <- fire_dyn[rownames(fire_dyn) %in% validARUs, ]
fire_dyn_dum <- array(NA, dim=c(nrow(fire_dyn), ncol(fire_dyn), 2))  # sites x years x (categories-1)

# Create dummy variables for each year
for(t in 1:ncol(fire_dyn)) {
  temp_dum <- model.matrix(~factor(fire_dyn[,t], levels=c("Unburned", "Low/Mod", "High")))[,-1]
  fire_dyn_dum[,t,] <- temp_dum
}

str(fire_dyn_dum)


win.rag <- list(
  y = y_long,
  y_wide = y_wide,
  btmax = static_clim$Tmax_Base_JJA,
  btmin = static_clim$Tmin_Base_JJA,
  atmax = static_clim$Tmax_Base_Annual,
  atmin = static_clim$Tmin_Base_Annual,
  bprec = static_clim$Prcp_Base_JJA,
  aprec = static_clim$Prcp_Base_Annual,
  trendt = static_clim$Trend_Tmax,
  trendtmin = static_clim$Trend_Tmin,
  trendp = static_clim$Trend_Prcp,
  Lat = static_clim$Lat,
  Long = static_clim$Long,
  cc = as.numeric(cancov),
  ch = as.numeric(canht),
  lf = as.numeric(lf),
  ele = as.numeric(topo),
  tanom = dyn_clim$Tmax,
  panom = dyn_clim$Prcp,
  fire = fire_dyn_dum,
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

saveRDS(win.rag, file = here("Data/Occ_Data/MSOM_Multi_Test_Categorical_Ragged_Full.RDS"))


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
  cc   = win.rag$cc,
  ch = win.rag$ch,
  lf = win.rag$lf,
  temp = win.rag$btmax,
  tmin = win.rag$btmin,
  prec = win.rag$bprec,
  atemp = win.rag$atmax,
  atmin = win.rag$atmin,
  aprec = win.rag$aprec,
  Lat = win.rag$Lat,
  Ele = win.rag$ele,
  TrendT = win.rag$trendt,
  TrendTmin = win.rag$trendtmin,
  TrendP = win.rag$trendp,
  Fire = fsev |> 
    filter(Cell_Unit %in% validARUs) |> 
    filter(Fire_Sev_SurvExtYr == 2021) |> 
    pull(Fire_Sev_mean_1_10)
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
ggplot(raw_occ_long, aes(x = scale(cc), y = detected)) +
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


ggplot(raw_occ_long, aes(x = temp, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
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

ggplot(raw_occ_long, aes(x = Fire, y = detected)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Fire Severity 1-10 Post",
       y = "Naïve Occupancy (≥1 detection in Year 1)")
