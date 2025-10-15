## -------------------------------------------------------------
##
## Script name: Multi-season Data Prep
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
load(here("./Data/Occ_Data/Thresh_By_Species/2021_99Conf_OccSppList.RData"))
sp.det.list21 <- sp.det.list$`Lazuli Bunting`
load(here("./Data/Occ_Data/Thresh_By_Species/2022_99Conf_OccSppList.RData"))
sp.det.list22 <- sp.det.list$`Lazuli Bunting`
load(here("./Data/Occ_Data/Thresh_By_Species/2023_99Conf_OccSppList.RData"))
sp.det.list23 <- sp.det.list$`Lazuli Bunting`
load(here("./Data/Occ_Data/Thresh_By_Species/2024_99Conf_OccSppList.RData"))
sp.det.list24 <- sp.det.list$`Lazuli Bunting`

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: 3-D array for site x rep x year
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Put this data into an array, but first we need to align the data
## so that we have an idea of the matching sites
## Add the padded 0 to ensure all cell units have C + 4 digits
sp.det.list21$Cell_Unit <- ifelse(stringr::str_detect(sp.det.list21$Cell_Unit, "C[0-9]{4}_"), sp.det.list21$Cell_Unit, gsub("C", "C0", sp.det.list21$Cell_Unit))
sp.det.list21 <- sp.det.list21[order(sp.det.list21$Cell_Unit), ]
sp.det.list22 <- sp.det.list22[order(sp.det.list22$Cell_Unit), ]
sp.det.list23 <- sp.det.list23[order(sp.det.list23$Cell_Unit), ]
sp.det.list24 <- sp.det.list24[order(sp.det.list24$Cell_Unit), ]

## Find all of the unique units across the 4 years
uniq_units <- c(unique(sp.det.list21$Cell_Unit), unique(sp.det.list22$Cell_Unit), unique(sp.det.list23$Cell_Unit), unique(sp.det.list24$Cell_Unit))

## These are the total unique units across the 4 survey years: 2095 ARUs
uniq_units <- unique(uniq_units)
length(uniq_units)

## We need to ensure that the number of unique units standardizes array for the occ model
template <- data.frame(Cell_Unit = uniq_units)

## Merge the dataframes for each year with this
sp.det.list <- list(sp.det.list21,
                    sp.det.list22,
                    sp.det.list23,
                    sp.det.list24)

sp.det.std <- lapply(sp.det.list, function(x) {
  
  sp.det.fill <- template |> 
    left_join(x) |> 
    #select(-Cell_Unit) |> 
    mutate(across(where(is.logical), as.numeric))
  return(sp.det.fill)
})

# ## Make this into an array
# nsite <- length(uniq_units)
# nsurv <- sum(str_detect(colnames(sp.det.list21), "^\\d{4}-"))
# nyear <- length(sp.det.std)
# y <- array(unlist(sp.det.std), dim = c(nsite, nsurv, nyear),
#            dimnames = list(uniq_units,
#                            gsub("^\\d{4}-", "", colnames(sp.det.list21)[colnames(sp.det.list21) != "Cell_Unit"]),
#                            c(2021:2024)))
# 
# ## Let's get a sense for the number of surveys that have been repeated
# table(nyear <- apply(y, c(1,2), function(x) sum(!is.na(x))))
# nvisits <- apply(y, c(1,3), function(x) sum(!is.na(x)))
# nvisits[nvisits > 0] <- 1
# table(rowSums(nvisits))
# table(nsites <- apply(y, c(2,3), function(x) sum(!is.na(x))))

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
samp.cols <- colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "\\d")]
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
sp.det.list.r <- lapply(sp.det.std, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = F))
sp.det.list.r <- lapply(sp.det.list.r, function(x) x |> dplyr::select(matches("J[0-9]")))

## Make this into an array
nsite <- length(uniq_units)
nsurv <- sum(str_detect(colnames(sp.det.list.r[[1]]), "J\\d{1}"))
nyear <- length(sp.det.list.r)
y <- array(unlist(sp.det.list.r), dim = c(nsite, nsurv, nyear),
           dimnames = list(uniq_units,
                           colnames(sp.det.list.r[[1]]),
                           c(2021:2024)))

## Let's get a sense for the number of surveys that have been repeated
# table(nyear <- apply(y, c(1,2), function(x) sum(!is.na(x))))
# nvisits <- apply(y, c(1,3), function(x) sum(!is.na(x)))
# nvisits[nvisits > 0] <- 1
# table(rowSums(nvisits))
# table(nsites <- apply(y, c(2,3), function(x) sum(!is.na(x))))


## Get the number of hours per survey for the detection covariate
## This needs to be edited to handle all years of data and 
## backfill in the sites that are missing with NA for NA covs
eff.dat <- list.files(here("Data/Occ_Data/Thresh_By_Species/"), pattern = "Effort", full.names = T)

eff.dat <- lapply(eff.dat, read.csv)

eff.dat <- lapply(eff.dat, function(x) x[,-1])

eff.dat <- lapply(eff.dat, function(x){
  colnames(x) <- gsub("[[:punct:]]", "_", gsub("X", "", colnames(x)))
  return(x)
})

eff.days <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Days"))
eff.days <- lapply(eff.days, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  #x <- x[,-1]
  #colnames(x) <- NULL
  return(x)
})

eff.hrs <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Hrs"))
eff.hrs <- lapply(eff.hrs, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
  #x <- x[,-1]
  #colnames(x) <- NULL
  return(x)
})

eff.jday <- lapply(eff.dat, function(x) second_samp(DAT = x, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "FirstJDay"))
eff.jday <- lapply(eff.jday, function(x){
  colnames(x)[colnames(x) == "V1"] <- "Cell_Unit"
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
y.new <- y
for (i in 1:dim(eff.days.std)[1]) {
  for (j in 1:dim(eff.days.std)[2]) {
    for(t in 1:dim(eff.days.std)[3]) {
      if (as.numeric(eff.days.std[i, j, t]) == 0) {
        # Set all layers (third dimension) of the array at position (i, j) to NA
        y.new[i, j, t] <- NA
      }
    }
  }
}
# Print the modified 3D array
print(y.new[,,1])
print(y.new[,,2])

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
#fsev <- fire$FireSeverity
fsev <- fire

## Combining output before reclassifying
# fsev <- fsev |> 
#   
# 
# fsev <- fsev |> 
#   mutate(FCat1_5 = case_when(Fire_Sev_mean_1_5 == 0 ~ "Unburned",
#                              Fire_Sev_mean_1_5 > 0 & Fire_Sev_mean_1_5 <= 2.25 ~ "Low/Mod",
#                              Fire_Sev_mean_1_5 > 2.25 ~ "High")) |> 
#   mutate(FCat6_10 = case_when(Fire_Sev_mean_6_10 == 0 ~ "Unburned",
#                               Fire_Sev_mean_6_10 > 0 & Fire_Sev_mean_6_10 <= 2.25 ~ "Low/Mod",
#                               Fire_Sev_mean_6_10 > 2.25 ~ "High")) |>
#   mutate(FCat1_10 = case_when(Fire_Sev_mean_6_10 == 0 & Fire_Sev_mean_1_5 == 0 ~ "Unburned",
#                               Fire_Sev_mean_6_10 > 2.25 | Fire_Sev_mean_1_5 > 2.25 ~ "High",
#                               Fire_Sev_mean_6_10 > 0 | Fire_Sev_mean_1_5 > 0 ~ "Low/Mod"))



## Load the anomaly data
clim <- read.csv(here("Data/Tmax_Anomaly_ARU_21_24.csv"))

## Merge the climate and fire data
fclim <- clim |> 
  rename("deployment_name" = "dplymn_") |> 
  #left_join(fsev, by = c("deployment_name", "srvy_yr" = "Fire_Sev_SurvExtYr")) |> 
  rename("Cell_Unit" = "Cll_Unt") |> 
  arrange(Cell_Unit) |> 
  #mutate(LagYr = srvy_yr - Year) |>
  #filter(LagYr == 1) |> 
  filter(Year %in% c(2020:2023)) |> 
  select(Cell_Unit, srvy_yr, Year, bclim, Anomaly, Lat, Long)

## Let's make this into an array to match the data
fclim.arr <- fclim |> 
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  tidyr::pivot_wider(names_from = Year,
                     values_from = Anomaly,
                     id_cols = "Cell_Unit",
                     values_fn = mean) |>
  tibble::column_to_rownames("Cell_Unit") |> 
  as.matrix()

## Array for the fire data
## I should remedy some of the issue by creating a unique set of 
##ARUs and extracting only from those for each year
## Duplications are a result of small shifts in ARU placement
fire.arr <- fsev |> 
  mutate(Cell_Unit = gsub("G\\d{3}_V\\d{1}_", "", deployment_name)) |> 
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  tidyr::pivot_wider(names_from = Fire_Sev_SurvExtYr,
                     values_from = Fire_Sev_mean_1_10,
                     id_cols = "Cell_Unit",
                     values_fn = mean) |>
  tibble::column_to_rownames("Cell_Unit") |>
  as.matrix()


## Take the baseline climate for initial occupancy
bclim <- fclim |> 
  filter(Cell_Unit %in% uniq_units) |> 
  select(Cell_Unit, Year, bclim, Lat, Long) |> 
  distinct() |> 
  tidyr::pivot_wider(names_from = Year,
                     values_from = bclim,
                     id_cols = "Cell_Unit",
                     values_fn = mean) |>
  tibble::column_to_rownames("Cell_Unit") |> 
  as.matrix()

## Forest data
cancov <- read.csv(here("Data/Cancov_ARU_21_24.csv"))
cancov <- cancov |> 
  rename("deployment_name" = "dplymn_") |> 
  rename("Cell_Unit" = "Cll_Unt") |>
  filter(Cell_Unit %in% uniq_units) |> 
  arrange(Cell_Unit) |> 
  select(Cell_Unit, cancov) |> 
  group_by(Cell_Unit) |> 
  summarise(cancov = mean(cancov)) |> 
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

win.data <- list(
  y = y.new,
  bclim = bclim,
  cc = cancov,
  anom = fclim.arr,
  fire = fire.arr,
  eff.hrs = eff.hrs.std.sc,
  eff.jday = eff.jday.std.sc,
  nsites = nsite,
  nyears = nyear,
  nreps = nsurv
)

str(win.data)

saveRDS(win.data, file = here("Data/Occ_Data/LABU_Multi_Test.RDS"))

## -------------------------------------------------------------
##
## Begin Section: Reformat and save data for Unmarked/UBMS model
##
## -------------------------------------------------------------

## Reshape detection data to matrix format (sites x years*visits)
y_ubms <- matrix(NA, nrow = nsite, ncol = nyear * nsurv)
for(t in 1:nyear) {
  col_start <- (t-1) * nsurv + 1
  col_end <- t * nsurv
  y_ubms[, col_start:col_end] <- y.new[,,t]
}

## Create site covariates dataframe
siteCovs <- data.frame(
  init_temp = bclim[,1],
  init_fire = fire.arr[,1]# Using first year as example, adjust as needed
)

## Create yearly site covariates (time-varying)
yearlySiteCovs <- list(
  temp_anomaly = fclim.arr[,1:ncol(fclim.arr), drop=FALSE],
  fire = fire.arr[,1:ncol(fclim.arr), drop=FALSE]
)


## Create observation covariates
## Need to reshape into matrices with sites x (years*visits)
effort <- matrix(NA, nrow = nsite, ncol = nyear * nsurv)
jday <- matrix(NA, nrow = nsite, ncol = nyear * nsurv)

for(t in 1:nyear) {
  col_start <- (t-1) * nsurv + 1
  col_end <- t * nsurv
  effort[, col_start:col_end] <- eff.hrs.std[,,t]
  jday[, col_start:col_end] <- eff.jday.std[,,t]
}

obsCovs <- list(
  effort = effort,
  jday = jday
)

## Create unmarkedMultFrame object
library(unmarked)
umf <- unmarkedMultFrame(
  y = y_ubms,
  siteCovs = siteCovs,
  yearlySiteCovs = yearlySiteCovs,
  obsCovs = obsCovs,
  numPrimary = nyear
)

# Save the unmarkedMultiFrame object
saveRDS(umf, file = here("Data/Occ_Data/LABU_Multi_UBMS.RDS"))

## Unmarked model fit
## Formula
(fm <- colext(psiformula = ~ init_temp + init_fire,
              gammaformula = ~ temp_anomaly + fire + temp_anomaly*fire,
              epsilonformula = ~ temp_anomaly + fire + temp_anomaly*fire,
              pformula = ~ effort + jday,
  umf))

smoothed(fm)

nd <- expand.grid(temp_anomaly = seq(-3, 6, length = 100), fire = seq(0, 3, length = 100))
eps <- predict(fm, type = "ext", newdata = nd, appendData = TRUE)
gmm <- predict(fm, type = "col", newdata = nd, appendData = TRUE)

ggplot(data = eps, aes(y = fire, x = temp_anomaly, fill = Predicted)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  xlab("Max Temp. Anomaly") + 
  ylab("Fire Severity") + 
  labs(fill = "Extinction Prob.")

ggplot(data = gmm, aes(y = fire, x = temp_anomaly, fill = Predicted)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  xlab("Max Temp. Anomaly") + 
  ylab("Fire Severity") + 
  labs(fill = "Colonization Prob.")



## UBMS model fit
fit_labu <- ubms::stan_colext(psiformula = ~ init_temp,
                              gammaformula = ~ temp_anomaly,
                              epsilonformula = ~ temp_anomaly,
                              pformula = ~ effort + jday,
                              data = umf,
                              chains = 3,
                              iter = 3000)
