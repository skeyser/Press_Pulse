## -------------------------------------------------------------
##
## Script name: Flocker Modelling
##
## Script purpose: Fitting occupancy models with Flocker for
## Sierra Nevada dataset.
##
## Author: Spencer R Keyser
##
## Date Created: 2025-04-01
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
library(stringr)

## This is the alternative false-positive model structure
#renv::install("jsocolar/flocker@alt-fp")
library(flocker)

## -------------------------------------------------------------

## *************************************************************
##
## Section Notes: What data does flocker need? What decisions
## need to be made?
## 1. Observation data
##  1.1. FP Model - 0 for non-det and probablity score for det
##  1.2. Matrix -> ncol = J sampling events; nrow = N sites * K species
##  1.3. FP Model input per Wood and Socolar et al., 2022
##  Steps: Categorical 
##
## *************************************************************

## Check out the data format
d <- simulate_flocker_data(n_rep = 6,
                           n_pt = 200,
                           n_sp = 30,
                           n_season = 1,
                           rep_constant = F)
str(d)

fd_rep_varying <- make_flocker_data(
  obs = d$obs, 
  unit_covs = d$unit_covs, 
  event_covs = d$event_covs,
  type = "single"
)




## -------------------------------------------------------------
##
## Begin Section: Species data processing
##
## -------------------------------------------------------------

## Load in the real Sierra data
## Remember the number of sites is conditional on the settings from
## occ_det_gen_file() function...
## setting no_det = 2 and eff_filter = 10 reduces the number of sites
## when combined with a month of sampling in this model
load(here("./Data/Generated_DFs/Occ_Mod_Data/Flocker/2021_90Conf_OccSppList.RData"))

## ADDED removing species and name change for PACFLY
sp.det.list <- sp.det.list[which(!names(sp.det.list) %in% c("Red-tailed Hawk",
                                                            "Osprey",
                                                            "Red-shouldered Hawk",
                                                            "Clark's Nutcracker",
                                                            "American Kestrel"))]

names(sp.det.list)[which(names(sp.det.list) == "Pacific-slope Flycatcher")] <- "Western Flycatcher"

## Cell Unit mapping file
cu.map <- data.frame(ID = 1:length(sp.det.list[[1]]$Cell_Unit), Cell_Unit = sp.det.list[[1]]$Cell_Unit)
cu.map <- cu.map |>
  mutate(Cell_Unit = ifelse(
    stringr::str_detect(string = Cell_Unit, pattern = "C[0-9]{3}"),
    gsub(pattern = "(C)([0-9]{3})(_U[0-9]+)$", replacement = "\\10\\2\\3", x = Cell_Unit),
    Cell_Unit
  ))

## Reduce the number of dates for sampling
## Check the duration of the sampling period
## 60 sampling periods
length(colnames(sp.det.list[[1]])) - 1

## Array for the data
## D1 (i) = Site, D2 (j) = Sampling Date, D3 (k) = species
samp.cols <- colnames(sp.det.list[[1]])[str_detect(colnames(sp.det.list[[1]]), "\\d")]
samp.cols <- as.Date(samp.cols, format = "%Y-%m-%d")
#samp.period <- interval(ymd("2021-06-01"), ymd("2021-06-30"))
#samp.ind <- as.character(samp.cols[which(samp.cols %within% samp.period)])
#sp.det.list <- lapply(sp.det.list, function(x) x |> select(1, all_of(samp.ind)))

## Create different sampling periods
second_samp <- function(DAT, 
                        interval, 
                        id_col,
                        binary = T,
                        eff = F, 
                        e.var = "Days"){
  tmp.cols <- colnames(DAT)[str_detect(colnames(DAT), "\\d")]
  tmp.splits <- split(tmp.cols, ceiling(seq_along(tmp.cols) / interval))
  tmp.new <- as.data.frame(cbind(DAT[, id_col], 
                                 matrix(data = NA, nrow = nrow(DAT), ncol = length(tmp.splits), 
                                        dimnames = list(NULL, paste0("J", seq(1:length(tmp.splits)))))))
  for(i in 1:length(tmp.splits)){
    if(!eff & binary){
      print("Binary Species Concatenation by any detection.")
      Jsum <- rowSums(DAT[, tmp.splits[[i]]], na.rm = T)
      Jsum <- ifelse(Jsum > 0, 1, 0)
    }  else if(!eff & !binary){
      print("Raw BirdNet Score Concatenation by max score.")
      Jsum <- apply(DAT[, tmp.splits[[i]]], 1, function(x) {
        if(all(is.na(x))) NA else max(x, na.rm = TRUE)
      })
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
sp.det.list.r <- lapply(sp.det.list, function(x) second_samp(DAT = x, interval = 6, binary = F, id_col = "Cell_Unit", eff = F))

samp.cols <- colnames(sp.det.list.r[[1]])[str_detect(colnames(sp.det.list.r[[1]]), "J")]
nsite <- nrow(sp.det.list.r[[1]])
nrep <- ncol(sp.det.list.r[[1]]) - 1
nspec <- length(names(sp.det.list.r))
sp.det <- lapply(sp.det.list.r, function(x) x |> select(all_of(samp.cols)))
y <- array(unlist(lapply(sp.det, as.matrix)), dim = c(nsite, nrep, nspec))
dimnames(y) <- list(NULL, NULL, names(sp.det.list))

## No NAs present in the data
table(nsurveys <- apply(y[,,1], 1, function(x) sum(!is.na(x))))

## Species with 0 occurrences
## 9 species
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == -Inf] <- NA

## Naive Occupancy Estimates
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE) / nrow(tmp))

drop.sp <- which(obs.occ == 0)
y <- y[,,-drop.sp]
sp.names <- dimnames(y)
sp.df <- data.frame(Index = 1:length(sp.names[[3]]),
                    Species = sp.names[[3]])
#write.csv(sp.df, file = here::here("Code/Occupancy_Modeling/SpeciesIndex_Filtered.csv"))

## Redefine nspec
nspec <- dim(y)[3]
# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum)) # Compute and print sorted species counts

## -------------------------------------------------------------
##
## End Section: Species data processing
##
## -------------------------------------------------------------


## -------------------------------------------------------------
##
## Begin Section: Model Covariates
##
## -------------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: EC (Event covariates) prep
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Morphological Characteristics
morph <- read.csv("Data/SierraBirds_Mass_Beak_PCA.csv")
morph <- morph |> 
  filter(Com_Name %in% sp.names[[3]]) |> 
  arrange(match(Com_Name, sp.names[[3]]))

bm <- morph$Mass
blc <- morph$Beak_Length_Culmen
bd <- morph$Beak_Depth
pc1 <- morph$PCA1
pc2 <- morph$PCA2

## Get the number of hours per survey for the detection covariate
eff.dat <- read.csv(here("./Data/Generated_DFs/Occ_Mod_Data/Flocker/2021_90Conf_OccEffortFileSubset.csv"))
eff.dat <- eff.dat[,-1]
colnames(eff.dat) <- gsub("[[:punct:]]", "_", gsub("X", "", colnames(eff.dat)))

## Subset the rows of the effort file to the rows of the species data
eff.dat <- eff.dat |> 
  filter(Cell_Unit %in% sp.det.list[[1]]$Cell_Unit)

## Summarize the effort data at the same temporal interval
eff.days <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Days")
colnames(eff.days)[colnames(eff.days) == "V1"] <- "Cell_Unit"
eff.days <- eff.days[,-1]
colnames(eff.days) <- NULL

## Summarize the total number of hours surveyed per sampling unit
eff.hrs <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "Hrs")
colnames(eff.hrs)[colnames(eff.hrs) == "V1"] <- "Cell_Unit"
#eff.hrs[eff.hrs == 0] <- NA
eff.hrs <- eff.hrs[,-1]
colnames(eff.hrs) <- NULL


## Summarize the median Jdate
eff.jday <- second_samp(DAT = eff.dat, interval = 6, id_col = "Cell_Unit", eff = T, e.var = "FirstJDay")
colnames(eff.jday)[colnames(eff.jday) == "V1"] <- "Cell_Unit"
#eff.jday[is.na(eff.days)] <- NA
eff.jday <- eff.jday[,-1]
colnames(eff.jday) <- NULL

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Filter obs by effort
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Set y to NA for days without sampling
# Iterate through the matrix and set corresponding values in the 3D array to NA
# Iterate through the matrix and set corresponding values in the 3D array to NA for all layers
for (i in 1:nrow(eff.days)) {
  for (j in 1:ncol(eff.days)) {
    if (as.numeric(eff.days[i, j]) == 0) {
      # Set all layers (third dimension) of the array at position (i, j) to NA
      y[i, j, ] <- NA
    }
  }
}

## Compile all the species into a single long DF for flocker
bind_3d_matrix <- function(array_3d) {
  result <- NULL
  for(i in 1:dim(array_3d)[3]) {
    slice <- array_3d[,,i]
    slice_df <- as.data.frame(slice)
    slice_df$Species <- dimnames(array_3d)[[3]][i]  # Add slice identifier
    result <- rbind(result, slice_df)
  }
  colnames(result) <- c(paste0("Obs_J", 1:ncol(array_3d)), "Species")
  return(result)
}

obs <- bind_3d_matrix(y)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Data checking
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Double check the different cutoffs for each species
print(y)

acwo <- y[,,1]
acwo <- apply(acwo, 1, function(x) max(x, na.rm = T))
acwo <- data.frame(Det = acwo, Cell_Unit = cu.map$Cell_Unit)
nrow(acwo[acwo$Det > 0,])/nrow(acwo)

## Load in the spatial coordinates
source(here("./Code/Acoustic_Data_Prep/Sierra_functions.R"))
aru_sf <- aru_sf_query(years = c(2021))
aru_sf <- aru_sf |> 
  filter(Cell_Unit %in% acwo$Cell_Unit) |>
  left_join(acwo) |> 
  mutate(Det_Cat = as.factor(case_when(Det < 0.161828 ~ "No Det",
                                       Det >= 0.161828 & Det < 0.329982 ~ "90",
                                       Det >= 0.329982 & Det < 0.795825 ~ "95",
                                       Det >= 0.795825 ~ "99"))) |> 
  mutate(Det_Cat2 = as.factor(case_when(Det < 0.9 ~ "No Det",
                                        Det >= 0.9 & Det < 0.95 ~ "90-95",
                                        Det >= 0.95 & Det < 0.975 ~ "95-97.5",
                                        Det >= 0.975 & Det < 0.99 ~ "97.5-99",
                                        Det >= 0.99 & Det < 1 ~ "99-100",
                                        Det == 1 ~ "100")))

mapview::mapview(aru_sf, zcol = "Det_Cat")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Shuffling NAs
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shuffle NAs to the back
# Assuming your dataframe is called 'df'
obs_Species <- unique(obs$Species)
obs <- obs[!colnames(obs) == "Species"]
colnames(obs) <- NULL
rownames(obs) <- NULL

## Moving all the NAs to the end
obs_shuffled <- t(apply(obs, 1, function(x) {
  c(x[!is.na(x)], x[is.na(x)])
}))

## Replicate with the event covariates
eff.hrs

eff_hrs_shuf <- t(apply(eff.hrs, 1, function(x) {
  c(x[x > 0], x[x == 0])
}))

eff_days_shuf <- t(apply(eff.days, 1, function(x) {
  c(x[x > 0], x[x == 0])
})) 

ec <- list(ec1 = eff_hrs_shuf, ec2 = eff_days_shuf)

## Make the
# If your dataframe is called 'df'
ec$ec1 <- do.call(rbind, replicate(length(obs_Species), ec$ec1, simplify = FALSE))
ec$ec2 <- do.call(rbind, replicate(length(obs_Species), ec$ec2, simplify = FALSE))

nrow(ec$ec1) == nrow(obs)

## Make all 0 effort NA in flocker
ec$ec1[ec$ec1 == 0] <- NA
ec$ec2[ec$ec2 == 0] <- NA

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: UC (unit covariate) prep
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aru_meta <- read.csv(here("./Data/ARU_120m.csv"))
aru_meta <- aru_meta |> 
  filter(cell_unit %in% cu.map$Cell_Unit) |> 
  arrange(match(cell_unit, cu.map$Cell_Unit)) |>
  select(Cell_Unit = cell_unit,
         Long = X, 
         Lat = Y,
         ele = topo_elev, 
         tpi = topo_tpi,
         tmx = tmx_bcm_mn, 
         ppt = ppt_bcm_mn,
         fire1yr_cbi_mn, 
         fire2_5yr_cbi_mn,
         fire6_10yr_cbi_mn, 
         fire11_35yr_cbi_mn,
         stage = standage_f3_mn, 
         cc = cpycovr_f3_mn) |> 
  mutate(fire1_5yr_cbi_mn = (fire1yr_cbi_mn + fire2_5yr_cbi_mn) / 2) |> 
  select(-fire1yr_cbi_mn)

## Scale the preds of interest
uc <- data.frame(
  ID = aru_meta$Cell_Unit,
  Long = scale(aru_meta$Long),
  Lat = scale(aru_meta$Lat),
  ele = scale(aru_meta$ele),
  ppt = scale(aru_meta$ppt),
  tmx = scale(aru_meta$tmx),
  cbi1_5 = scale(aru_meta$fire1_5yr_cbi_mn),
  cbi6_10 = scale(aru_meta$fire6_10yr_cbi_mn),
  cbi11_35 = scale(aru_meta$fire11_35yr_cbi_mn),
  stage = scale(aru_meta$stage), 
  cc = scale(aru_meta$cc)
)

uc_rep <- do.call(rbind, replicate(length(obs_Species), uc, simplify = FALSE))
uc_rep$species <- rep(obs_Species, each = nrow(uc))

## -------------------------------------------------------------
##
## Begin Section: Create the blocking for categorical models
##
## -------------------------------------------------------------
## *************************************************************
##
## Section Notes: Making bins per species-specific thresholds
## this is only done for the detections
## 1. 0.9 - 0.95
## 2. 0.95 - 0.975
## 3. 0.975 - 0.99
## 4. 1
## *************************************************************
## Load the threshold file
thresh <- read.csv(here("./Data/Thresholds_2021_20230309.csv"))
thresh <- thresh |> 
  select(species, cutoff90.r_conf,
         cutoff95.r_conf, cutoff99.r_conf) |> 
  mutate(species = ifelse(species == "Pacific-slope Flycatcher", "Western Flycatcher", species))

mod.block <- obs
mod.block$block <- uc_rep$species
colnames(mod.block) <- c(paste0("J", 1:(ncol(mod.block)-1)), "block")


mod.block <- mod.block |>
  tidyr::pivot_longer(cols = contains("J"),
                      names_to = "Event",
                      values_to = "y") |> 
  group_by(block) |> 
  filter(!is.na(y) & y > 0) |> 
  left_join(thresh, by = c("block" = "species")) |> 
  mutate(y_bin = case_when(y < cutoff95.r_conf ~ 1,
                           y >= cutoff95.r_conf & y < cutoff99.r_conf ~ 2,
                           y >= cutoff99.r_conf & y < 1 ~ 3,
                           y == 1 ~ 4)) |> 
  mutate(block = as.factor(block))

table(mod.block$y_bin)
table(mod.block$block)
levels(mod.block$block)
str(mod.block)
quantile(table(mod.block$block))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: BRMS model for species binning
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Subsampling to speed model fits while retaining
## full information for species with low SS
## Define threshold for "low" sample size
low_threshold <- 310  # adjust as needed

# Create two groups and then bind them
low_n_species <- mod.block |>
  group_by(block) |>
  mutate(n_obs = n()) |>
  filter(n_obs <= low_threshold) # keep all low n species

high_n_species <- mod.block |>
  group_by(block) |>
  mutate(n_obs = n()) |>
  filter(n_obs > low_threshold) |>
  group_by(block) |>
  slice_sample(n = 500) # subsample high n species

# Combine the datasets
mod.block.sub <- bind_rows(low_n_species, high_n_species)
table(mod.block.sub$block)
str(mod.block.sub)

## Fit the model
bin_mod <- brm(
  y_bin ~ 1 + (1|block), 
  family = categorical(),
  backend = "cmdstanr", 
  cores = 4, 
  data = mod.block.sub
)

bin_mod
#saveRDS(bin_mod, here("./Data/Flocker/Model_Output/CategoricalModOut_Subsampled.RDS"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Prepare the FP Data for flocker
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load the model fit from the subsampled binned data
out <- readRDS(here("./Data/Flocker/Model_Output/CategoricalModOut_Subsampled.RDS"))

## Make an empty DF
sp.vals <- data.frame(
  species = obs_Species
) |> mutate(
  block = NA,
  frac = NA,
  frac1 = NA,
  pp_1 = NA,
  pp_2 = NA,
  pp_3 = NA,
  pp_4 = NA,
  qq_1 = NA,
  qq_2 = NA,
  qq_3 = NA,
  qq_4 = NA
)

## Get species mean confidence scores
for(i in seq_len(nrow(sp.vals))){
  a <- unique(mod.block$block[mod.block$block == sp.vals$species[i]])
  assertthat::assert_that(length(a) == 1)
  sp.vals$block[i] <- a
  sp.vals$frac[i] <- mean(mod.block$y[mod.block$block == sp.vals$species[i]])
}

# Get species-specific predictions
newdata <- data.frame(block = unique(mod.block.sub$block)) |>
  mutate(block = as.factor(block))

# Get posterior predictions
epreds <- posterior_epred(
  out, 
  newdata = newdata,
  re_formula = NULL  # This ensures random effects are included
)

# Calculate mean probabilities per species
species_probs <- apply(epreds, c(2,3), mean)

# Create readable output
species_summary <- data.frame(
  species = newdata$block,
  prob_bin1 = species_probs[,1],
  prob_bin2 = species_probs[,2],
  prob_bin3 = species_probs[,3],
  prob_bin4 = species_probs[,4]
)

# View results
head(species_summary)

## Generate midpoint values for probabilities above/below scores
midpoints <- c(0.925, 0.97, 0.995, 1)
for(i in seq_len(nrow(sp.vals))){
  sp.vals$frac1[i] <- species_probs[i,] %*% midpoints
  pp <- qq <- rep(NA, length(midpoints))
  for(j in 1:length(midpoints)){
    pp[j] <- species_probs[i, j] * midpoints[j]
    qq[j] <- species_probs[i, j] * (1 - midpoints[j])
  }
  sp.vals$pp_1[i] <- pp[1]/sum(pp)
  sp.vals$qq_1[i] <- qq[1]/sum(qq)
  sp.vals$pp_2[i] <- pp[2]/sum(pp)
  sp.vals$qq_2[i] <- qq[2]/sum(qq)
  sp.vals$pp_3[i] <- pp[3]/sum(pp)
  sp.vals$qq_3[i] <- qq[3]/sum(qq)
  sp.vals$pp_4[i] <- pp[4]/sum(pp)
  sp.vals$qq_4[i] <- qq[4]/sum(qq)
}

## Adding 
frac <- pp <- qq <- flocker:::new_array(obs_shuffled, -99)
for(i in 1:dim(obs_shuffled)[1]){
  species <- uc_rep$species[i]
  tmp.thresh <- thresh[thresh$species == species,]
  for(j in 1:dim(obs_shuffled)[2]){
    frac[i,j] <- sp.vals$frac1[sp.vals$species == species]
    if(is.na(obs_shuffled[i,j])){
      next
    }
    if(obs_shuffled[i,j] < tmp.thresh$cutoff90.r_conf){
      obs_shuffled[i,j] <- 0
    } else {
      pp[i,j] <- sp.vals$pp_1[sp.vals$species == species]
      qq[i,j] <- sp.vals$qq_1[sp.vals$species == species]
    }
    if(obs_shuffled[i,j] > tmp.thresh$cutoff95.r_conf){
      pp[i,j] <- sp.vals$pp_2[sp.vals$species == species]
      qq[i,j] <- sp.vals$qq_2[sp.vals$species == species]
    }
    if(obs_shuffled[i,j] > tmp.thresh$cutoff99.r_conf){
      pp[i,j] <- sp.vals$pp_3[sp.vals$species == species]
      qq[i,j] <- sp.vals$qq_3[sp.vals$species == species]
    }
    if(obs_shuffled[i,j] == 1){
      pp[i,j] <- sp.vals$pp_4[sp.vals$species == species]
      qq[i,j] <- sp.vals$qq_4[sp.vals$species == species]
    }
  }
}

block <- flocker:::new_array(obs)
for(i in 1:dim(obs)[1]){
  block[i,] <- sp.vals$block[sp.vals$species == uc_rep$species[i]]
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Make flocker data object
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Bundle this up
fp_data <- list(
  block = block,
  frac_true = frac,
  P = pp,
  QQ = qq
)

#uc_rep$species <- as.character(as.factor(uc_rep$species))

## Try to package
fd <- make_flocker_data(obs = obs_shuffled,
                        unit_covs = uc_rep,
                        event_covs = ec,
                        type = "single",
                        fp = TRUE,
                        fp_data = fp_data
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Fit flocker model
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fp_occ_mod <- flock(
  f_occ = ~ 0 + Intercept + (1 | g1 | species),
  f_det = ~ 0 + Intercept + ec1 + ec2 + (1 | g1 | species),
  fp = TRUE,
  flocker_data = fd,
  cores = 4, 
  backend = "cmdstanr"
)