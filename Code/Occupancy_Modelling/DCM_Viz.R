## ----------------------------------------------------------
##
## Script name: DCM Model Summaries and Plots
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-10-16
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
library(here)
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(coda)
library(MCMCvis)
library(tidyr)
library(tidybayes)
library(stringr)
library(patchwork)

## -----------------------------------------------------------

source(here("./Code/Occupancy_Modelling/DCM_Viz_Funs.R"))
source(here("./Code/Occupancy_Modelling/DCM_Viz_Funs_Seasons.R"))

## Model object with info
mod.obj <- readRDS(here("Data/Occ_Data/DCM_Ragged_Full_2021_2025_FireCont_Thresh24.RDS")) 
mod.obj <- readRDS(here("Data/Occ_Data/DCM_Ragged_Full_2021_2024_RevisedCovs_Tester4yrOnly.RDS"))
mod.obj <- readRDS(here("Data/Occ_Data/DCM_Ragged_Complete_Survey_Only.RDS"))

specs <- dimnames(mod.obj$y_wide)$Species

## Read in saved objects
# samples <- readRDS("D:/DCM_Samples/DCMmodel_mcmc_output.rds")
# samplesFP <- readRDS("D:/DCM_Samples/DCMmodel_mcmc_output_FP_21_25_Thresh24_NoTraits_V2_Simplified.rds")
# samplesFP <- readRDS("D:/DCM_Samples/DCMmodel_mcmc_output_FP_21_25_Thresh24_NoTraits_V3_Parallel.rds")
# samplesFP <- readRDS("D:/DCM_Samples/DCMmodel_mcmc_output_FP_21_24_Spring_Climate.rds")
samplesFP <- readRDS("D:/DCM_Samples/DCMmodel_mcmc_output_FP_21_25_Summer_Climate_CompSurv.rds")
samplesFP <- readRDS("D:/DCM_Samples/DCMmodel_mcmc_output_FP_21_25_Spring_Climate_CompSurv.rds")

mcmc_list_fp <- mcmc.list(
  mcmc(samplesFP$chain1), 
  mcmc(samplesFP$chain2), 
  mcmc(samplesFP$chain3)
)

rm(samplesFP)
gc()

any(MCMCvis::MCMCsummary(mcmc_list_fp,
                     params = 'gamma',
                     round = 2,
                     exact = F,
                     probs = c(0.055, 0.5, 0.945))$Rhat > 1.1)
gamma_summary <- MCMCvis::MCMCsummary(mcmc_list_fp,
                                      params = 'gamma',
                                      round = 2,
                                      exact = F,
                                      probs = c(0.055, 0.5, 0.945))
gsum <- gamma_summary[gamma_summary$Rhat > 1.1, ]

any(MCMCvis::MCMCsummary(mcmc_list_fp,
                         params = 'eps',
                         round = 2,
                         exact = F,
                         probs = c(0.055, 0.5, 0.945))$Rhat > 1.1)
eps_summary <- MCMCvis::MCMCsummary(mcmc_list_fp,
                                    params = 'eps',
                                    round = 2,
                                    exact = F,
                                    probs = c(0.055, 0.5, 0.945))

eps_summary[eps_summary$Rhat > 1.1, ]

MCMCvis::MCMCtrace(mcmc_list_fp, params = "eps10")

any(MCMCvis::MCMCsummary(mcmc_list_fp,
                         params = 'beta',
                         round = 2,
                         exact = F,
                         probs = c(0.055, 0.5, 0.945))$Rhat > 1.1)
beta_summary <- MCMCvis::MCMCsummary(mcmc_list_fp,
                                     params = 'beta',
                                     round = 2,
                                     exact = F,
                                     probs = c(0.055, 0.5, 0.945))

beta_summary[beta_summary$Rhat > 1.1, ]

any(MCMCvis::MCMCsummary(mcmc_list_fp,
                         params = 'alpha',
                         round = 2,
                         exact = F,
                         probs = c(0.055, 0.5, 0.945))$Rhat > 1.1)
alpha_summary <- MCMCvis::MCMCsummary(mcmc_list_fp,
                                      params = 'alpha',
                                      round = 2,
                                      exact = F,
                                      probs = c(0.055, 0.5, 0.945))
alpha_summary[alpha_summary$Rhat > 1.1, ]
plogis(alpha_summary[grep("^alpha0\\[", rownames(alpha_summary)),]$mean)

any(MCMCvis::MCMCsummary(mcmc_list_fp,
                         params = 'mu.',
                         round = 2,
                         exact = F,
                         probs = c(0.055, 0.5, 0.945))$Rhat > 1.1)
mu.summary <- MCMCvis::MCMCsummary(mcmc_list_fp,
                                   params = 'mu.',
                                   round = 2,
                                   exact = F,
                                   probs = c(0.055, 0.5, 0.945))
mu.summary[mu.summary$Rhat > 1.1, ]

any(MCMCvis::MCMCsummary(mcmc_list_fp,
                         params = 'alphaFP',
                         round = 2,
                         exact = F,
                         probs = c(0.055, 0.5, 0.945))$Rhat > 1.1)

fp.summary <- MCMCvis::MCMCsummary(mcmc_list_fp,
                                   params = 'alphaFP',
                                   round = 2,
                                   exact = F,
                                   probs = c(0.055, 0.5, 0.945))


any(MCMCvis::MCMCsummary(mcmc_list_fp,
                         params = c('alpha_year'),
                         round = 2,
                         exact = F,
                         probs = c(0.055, 0.5, 0.945))$Rhat > 1.1)

# Extract alpha_year posterior means
alpha_year_summary <- MCMCvis::MCMCsummary(mcmc_list_fp, 
                                           params = 'alpha_year',
                                           round = 3) %>%
  rownames_to_column("parameter") %>%
  # Parse the parameter names to extract species and year
  tidyr::extract(parameter, 
                 into = c("species", "year"), 
                 regex = "alpha_year\\[(\\d+), (\\d+)\\]") %>%
  mutate(species = as.numeric(species),
         year = as.numeric(year)) %>%
  select(species, year, mean, `2.5%`, `97.5%`)

# Create the line plot
ggplot(alpha_year_summary, aes(x = year, y = mean, color = factor(species))) +
  geom_line(alpha = 0.7) +
  geom_point(size = 0.5, alpha = 0.7) +
  labs(x = "Year", 
       y = "Detection Effect (log-odds scale)",
       title = "Species-Specific Yearly Detection Effects",
       subtitle = "Each line represents one species across years") +
  theme_minimal() +
  theme(legend.position = "none") +  # Too many species for useful legend
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(breaks = 1:5)  # Assuming you have 4 years

# Alternative: Add confidence intervals for a few select species
select_species <- c(7, 18, 40)  # The ones with big drops you mentioned

ggplot(alpha_year_summary, aes(x = year, y = mean)) +
  geom_line(aes(group = species), alpha = 0.3, color = "gray60") +
  geom_line(data = filter(alpha_year_summary, species %in% select_species),
            aes(color = factor(species)), linewidth = 1) +
  geom_ribbon(data = filter(alpha_year_summary, species %in% select_species),
              aes(ymin = `2.5%`, ymax = `97.5%`, fill = factor(species)), 
              alpha = 0.2) +
  labs(x = "Year", 
       y = "Detection Effect (log-odds scale)",
       title = "Yearly Detection Effects: All Species (gray) vs Selected Species (colored)",
       color = "Species", fill = "Species") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(breaks = 1:5)


## Yearly differences in detection prob by year and species seems to be important
## based on estimates of the parameter
spec <- 11
MCMCvis::MCMCsummary(mcmc_list_fp,
                     params = paste0('alpha_year\\', "[", spec, ", ", 1:5, "\\]"),
                     #params = "alpha_year",
                     round = 2,
                     exact = F,
                     ISB = F,
                     probs = c(0.055, 0.5, 0.945))

## ACF plots
coda::autocorr.plot(mcmc_list_fp[, "gamma0[1]"])
coda::autocorr.plot(mcmc_list_fp[, "delta0.eps0"])
coda::autocorr.plot(mcmc_list_fp[, "alphaFP"])




## Community Summaries
MCMCvis::MCMCplot(mcmc_list_fp,
                     params = "mu.gamma",
                     exact = F,
                  ref_ovl = T,
                  ref = 0,
                  ci = c(50, 89),
                  ylab = "Parameter")

## Traits
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "delta1.eps",
                  exact = F,
                  ci = c(50, 89),
                  rank = F,
                  ref = 0,
                  ref_ovl = T,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "delta2.eps",
                  exact = F,
                  ci = c(50, 89),
                  rank = F,
                  ref = 0,
                  ref_ovl = T,
                  guide_lines = T)

## Detection estimates
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "alpha0",
                  exact = F,
                  ci = c(50, 95),
                  rank = F,
                  ref = 0,
                  ref_ovl = T,
                  ylab = "Detection Probability (logit)")

## Occupancy
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "beta0",
                  excl = "mu.beta0",
                  labels = specs,
                  sz_labels = 0.5,
                  exact = F,
                  ci = c(50, 89),
                  rank = T,
                  ylab = "Occupancy Probability (logit)")

## Colonization prob. mean
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "^gamma0",
                  exact = F,
                  ci = c(50, 89),
                  labels = specs,
                  sz_labels = 0.5,
                  rank = T,
                  ylab = "Colonization Probability (logit)")

## Extinction prob. mean
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "^eps0",
                  exact = F,
                  ci = c(50, 89),
                  labels = specs,
                  sz_labels = 0.5,
                  rank = T,
                  ylab = "Extinction Probability (logit)")

## False-positive rates
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "alphaFP",
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "FP Rate",
                  ref = 0.05)

MCMCvis::MCMCtrace(mcmc_list_fp,
                   params = "alphaFP",
                   pdf = F)

MCMCvis::MCMCtrace(mcmc_list_fp,
                   params = "mu.beta0")

## Year effects
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = "tau.alpha_year",
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Variance of Year Effect",
                  ref = 0.05)

## Initial occupancy
par(mfrow = c(1,2))
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'beta1',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Baseline T",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'beta2',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Baseline T sq",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'beta3',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Baseline P",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'beta4',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Baseline P sq",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'beta5',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Canopy Cover",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'beta6',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Canopy Cover Sq",
                  ref = 0)

## Colonization coefficients
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma1',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "T Anomaly",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma2',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "P Anomaly",
                  ref = 0)

par(mfrow = c(1,2))
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma3',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "TSF HSF",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma4',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Tanom x Fire (Col)",
                  ref = 0,
                  guide_lines = T)

par(mfrow = c(1,2))
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma5',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "P Anom x Fire",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma6',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Trend Tmax",
                  ref = 0,
                  guide_lines = T)


par(mfrow = c(1,2))
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma7',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Trend Tmin",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma8',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Trend Precip (Col)",
                  ref = 0,
                  guide_lines = T)

dev.off()
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma9',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "T Max Trend x Fire",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma10',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "T Min Trend x Fire",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'gamma11',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Precip Trend x Fire",
                  ref = 0,
                  guide_lines = T)

## Extinction coefficients
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps1',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "T Anomaly (Ext)",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps2',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "P Anomaly (Ext)",
                  ref = 0)

par(mfrow = c(1,2))
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps3',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "HSF PLAND (Ext)",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps4',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Tanom x Fire (Ext)",
                  ref = 0,
                  guide_lines = T)

par(mfrow = c(1,2))
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps5',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "P Anom * Fire (Ext)",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps6',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "T Trend (Ext)",
                  ref = 0,
                  guide_lines = T)


par(mfrow = c(1,2))
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps7',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Tmin Trend (Ext)",
                  ref = 0,
                  guide_lines = T)


MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps8',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Trend Precip (Ext)",
                  ref = 0,
                  guide_lines = T)


MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps9',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "T Max Trend x Fire",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps10',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "T Min Trend x Fire",
                  ref = 0,
                  guide_lines = T)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'eps11',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.75,
                  ylab = "Species",
                  xlab = "Precip Trend x Fire",
                  ref = 0,
                  guide_lines = T)

## Detection effects
MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'alpha1',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Effort Hrs",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'alpha2',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "J Day",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'alpha3',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = T,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "J Day^2",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'alpha4',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  labels = specs,
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Canopy Cover (Detection)",
                  ref = 0)

MCMCvis::MCMCplot(mcmc_list_fp,
                  params = 'sd.cell_det',
                  exact = T,
                  ISB = T,
                  ci = c(50, 89),
                  ref_ovl = T,
                  rank = F,
                  sz_labels = 0.5,
                  ylab = "Species",
                  xlab = "Cell-level Ranef (Detection)",
                  ref = 0)

## Better plotting from the posteriors
## ***********************************************************
##
## Section Notes:
## To-Do
## 1. Community responses in detection, occ, col/ext (Half-eyes)
## 2. Species-specific responses (Eff plots)
## 3. Species-specific responses (Spaghetti plots)
## 4. Species Ext vs. Col rates (differences for "declining" species)
## 5. Predictive surfaces (maps)
##
## ***********************************************************

## Line plots for the species
## Generate prediction coefficients

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Bundle up the draws for parameters of interest
## to keep it light we will only take occ, col/ex, and detection
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Make a big dataframe with the output
n_chains <- 3
n_iter <- nrow(mcmc_list_fp[[1]])
str(mcmc_list_fp)
samples_df <- mcmc_list_fp |> 
  as.matrix() |> 
  as_tibble() |> 
  select(contains(c("beta", "gamma", "alpha", "eps"))) |> 
  mutate(.chain = rep(1:n_chains, each = n_iter),
         .iteration = rep(1:n_iter, n_chains)) |> 
  pivot_longer(cols = -c(.chain, .iteration),
               names_to = "variable",
               values_to = "value")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Community Means
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Mapping variable names
var_names <- c("gamma", "eps", "beta", "alpha")
var_map <- data.frame(var_names, param = 0:11)
var_map <- var_map |>
  expand(param, var_names) |> 
  arrange(var_names) |> 
  mutate(var_types = case_when(var_names == "alpha" ~ "Detection",
                               var_names == "beta" ~ "Occupancy",
                               var_names == "gamma" ~ "Colonization",
                               var_names == "eps" ~ "Extinction"))

## Add the predictor names
## Change equals to str_detect to accommodate the flexible
## indexing added to the end of non-community vectors
var_map <- var_map |> 
  mutate(PredPretty = case_when(var_names == "alpha" & param == 0 ~ "Intercept",
                                var_names == "alpha" & param == 1 ~ "Effort Hours",
                                var_names == "alpha" & param == 2 ~ "Ordinal Date",
                                var_names == "alpha" & param == 3 ~ "Ordinal Date^2",
                                var_names == "alpha" & param == 4 ~ "Canopy Cover",
                                var_names == "beta" & param == 0 ~ "Intercept",
                                var_names == "beta" & param == 1 ~ "Annual Tmax",
                                var_names == "beta" & param == 2 ~ "Annual Tmax^2",
                                var_names == "beta" & param == 3 ~ "Annual Precip",
                                var_names == "beta" & param == 4 ~ "Annual Precip^2",
                                var_names == "beta" & param == 5 ~ "Canopy Cover",
                                var_names == "beta" & param == 6 ~ "Canopy Cover^2",
                                var_names == "gamma" & param == 0 ~ "Intercept",
                                var_names == "gamma" & param == 1 ~ "Tmax Anomaly",
                                var_names == "gamma" & param == 2 ~ "Ppt Anomaly",
                                var_names == "gamma" & param == 3 ~ "HSF",
                                var_names == "gamma" & param == 4 ~ "Tmax Anom x HSF",
                                var_names == "gamma" & param == 5 ~ "Ppt Anom x HSF",
                                var_names == "gamma" & param == 6 ~ "Tmax Trend",
                                var_names == "gamma" & param == 7 ~ "Tmin Trend",
                                var_names == "gamma" & param == 8 ~ "Ppt Trend",
                                var_names == "gamma" & param == 9 ~ "Tmax Trend x HSF",
                                var_names == "gamma" & param == 10 ~ "Tmin Trend x HSF",
                                var_names == "gamma" & param == 11 ~ "Ppt Trend x HSF",
                                var_names == "eps" & param == 0 ~ "Intercept",
                                var_names == "eps" & param == 1 ~ "Tmax Anomaly",
                                var_names == "eps" & param == 2 ~ "Ppt Anomaly",
                                var_names == "eps" & param == 3 ~ "HSF",
                                var_names == "eps" & param == 4 ~ "Tmax Anom x HSF",
                                var_names == "eps" & param == 5 ~ "Ppt Anom x HSF",
                                var_names == "eps" & param == 6 ~ "Tmax Trend",
                                var_names == "eps" & param == 7 ~ "Tmin Trend",
                                var_names == "eps" & param == 8 ~ "Ppt Trend",
                                var_names == "eps" & param == 9 ~ "Tmax Trend x HSF",
                                var_names == "eps" & param == 10 ~ "Tmin Trend x HSF",
                                var_names == "eps" & param == 11 ~ "Ppt Trend x HSF",
                                )) |> 
  drop_na() |> 
  mutate(variable = paste0(var_names, param),
         cvariable = paste0("mu.", var_names, param))

## Sample only the mean responses
com.means <- samples_df |> 
  filter(str_detect(variable, "mu")) |> 
  left_join(var_map |> select(-variable), by = c("variable" = "cvariable"))

com.means |>
  filter(str_detect(variable, "gamma|eps")) |> 
  filter(PredPretty != "Intercept") |> 
  ggplot(aes(x = value, y = PredPretty)) +
  stat_halfeye(
    .width = c(0.50, 0.89),
    point_interval = median_qi,
    fill = "steelblue",
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    x = "Parameter Value",
    y = "Parameter"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 14)
  ) + 
  facet_wrap(~var_types, ncol = 2)

com.means |>
  filter(str_detect(variable, "beta")) |> 
  ggplot(aes(x = value, y = PredPretty)) +
  stat_halfeye(
    .width = c(0.50, 0.89),
    point_interval = median_qi,
    fill = "steelblue",
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Parameter Value",
    y = "Parameter"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 14)
  )

com.means |>
  filter(str_detect(variable, "alpha")) |> 
  ggplot(aes(x = value, y = PredPretty)) +
  stat_halfeye(
    .width = c(0.50, 0.89),
    point_interval = median_qi,
    fill = "steelblue",
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Parameter Value",
    y = "Parameter"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 14)
  )

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Species-specific effects plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sp.id <- data.frame(species = specs, sp.ind = as.character(1:length(specs)))

## Sample only the mean responses
col.means <- samples_df |> 
  filter(str_detect(variable, "mu|sd", negate = TRUE)) |>
  filter(str_detect(variable, "gamma[\\d+]")) |>
  separate(col = variable, into = c("variable", "sp.ind"), sep = "\\[") |> 
  mutate(sp.ind = gsub("\\]", "", sp.ind)) |>
  left_join(var_map |> select(-cvariable), by = c("variable" = "variable")) |> 
  left_join(sp.id, by = "sp.ind")

p.col <- col.means |> 
  filter(variable == "gamma0") |> 
  ggplot(aes(x = plogis(value), y = species)) +
  stat_halfeye(
    .width = c(0.50, 0.89),
    point_interval = median_qi,
    fill = "darkgreen",
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Parameter Estimate",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 14)
  )

## Sample only the mean responses
ext.means <- samples_df |> 
  filter(str_detect(variable, "mu|sd", negate = TRUE)) |>
  filter(str_detect(variable, "eps[\\d+]")) |>
  separate(col = variable, into = c("variable", "sp.ind"), sep = "\\[") |> 
  mutate(sp.ind = gsub("\\]", "", sp.ind)) |>
  left_join(var_map |> select(-cvariable), by = c("variable" = "variable")) |> 
  left_join(sp.id, by = "sp.ind")

p.ext <- ext.means |> 
  filter(variable == "eps0") |> 
  ggplot(aes(x = plogis(value), y = species)) +
  stat_halfeye(
    .width = c(0.50, 0.89),
    point_interval = median_qi,
    fill = "darkred",
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Parameter Estimate",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 14)
  )

p.col + p.ext

## Plot col/ext together
dyn.mean <- rbind(col.means, ext.means) 

dyn.mean |> 
  filter(PredPretty == "Intercept") |> 
  ggplot(aes(x = plogis(value), y = species, fill = var_types)) +
  stat_halfeye(
    .width = c(0.50, 0.89),
    point_interval = median_qi,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("darkgreen", "darkred")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(
    x = "Parameter Estimate",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 14)
  )


## Spaghetti plots to help visualize
bdata <- readRDS(here("Data/Occ_Data/DCM_Ragged_Full_2021_2025_FireCont_Thresh24.RDS"))  
bdata <- readRDS(here("Data/Occ_Data/DCM_Ragged_Full_2021_2024_RevisedCovs_Tester4yrOnly.RDS"))  
bdata <- readRDS(here("Data/Occ_Data/DCM_Ragged_Complete_Survey_Only.RDS"))  


# Modified parameter compilation to include uncertainty
param.compile <- samples_df |> 
  group_by(variable) |> 
  summarise(
    mean = mean(value),
    median = median(value),
    lowci = quantile(value, 0.055),
    hici = quantile(value, 0.945),
    lower_80 = quantile(value, 0.1),
    upper_80 = quantile(value, 0.9),
    .groups = 'drop'
  ) |> 
  filter(str_detect(variable, "^beta|^gamma|^eps")) |> 
  separate(col = variable, into = c("variable", "sp.ind"), sep = "\\[") |> 
  mutate(sp.ind = as.numeric(gsub("\\]", "", sp.ind))) |> 
  arrange(variable, sp.ind)

# Community-level compilation with uncertainty
com.compile <- samples_df |> 
  group_by(variable) |> 
  summarise(
    mean = mean(value),
    median = median(value),
    lowci = quantile(value, 0.055),
    hici = quantile(value, 0.945),
    lower_80 = quantile(value, 0.1),
    upper_80 = quantile(value, 0.9),
    .groups = 'drop'
  ) |> 
  filter(str_detect(variable, "^mu.eps"))


## Testing the function
me <- dynOcc_margeff(var = "tmax_trend_mam",
               model_component = "eps",
               param_compile = param.compile,
               raw_data = bdata,
               level = "species",
               pred_point = 100,
               climate_season = "mam")

## Unlist and stack
me <- do.call(rbind, me)

me_plots(me, 
         x_var = "Value_original", 
         group_only = F,
         line_color = "grey",
         line_alpha = 0.5,
         line_size = 0.8,
         x_label = "Tmax",
         y_label = expression(sigma),
         show_legend = F)

## Interactions
dyn_type <- "gamma"
dyn_lab <- ifelse(dyn_type == "eps", "Ext. Prob", "Col. Prob")

interaction_data <- dynOcc_interaction(
  var1 = "hsf_pland1_10",
  var2 = "tmin_trend_mam",
  model_component = dyn_type,
  param_compile = param.compile,
  raw_data = bdata,
  var2_levels = 3,
  climate_season = "mam"
)

# MCMCvis::MCMCplot(mcmc_list_fp, params = "eps4", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "eps5", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "eps9", ci = c(50, 89),  ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "eps10", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "eps11", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "gamma4", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "gamma5", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "gamma9", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "gamma10", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)
# MCMCvis::MCMCplot(mcmc_list_fp, params = "gamma11", ci = c(50, 89), ref = 0, ref_ovl = T, labels = specs, sz_labels = 0.5)

sp_int <- me_interaction_plots(dat_list = interaction_data, plot_type = "lines", selected_species = 1:51,
                               species_names = specs, show_significance = T, param_compile = param.compile, 
                               model_component = dyn_type,
                               y_label = dyn_lab, x_label = "Prop. HSF 1-10 yr",
                               palette = "viridis", climate_season = "mam", use_original_scale = TRUE)

me_interaction_plots(dat_list = interaction_data, plot_type = "lines", selected_species = 1)

me_interaction_plots(dat_list = interaction_data, plot_type = "heatmap", selected_species = 4)


param_sig <- param_compile |> 
  mutate(Sig = case_when(hici < 0 ~ "Negative",
                         lowci > 0 ~ "Positive",
                         TRUE ~ "Non")) |> 
  group_by(variable) |> 
  count(Sig)

ex.grid <- expand_grid(variable = unique(param_sig$variable), Sig = unique(param_sig$Sig))

facet_names <- c("Extinction" = "epsilon",
                 "Colonization" = "gamma",
                 "Occupancy" = "beta")

eff.sum <- param_sig |>
  ungroup() |> 
  full_join(ex.grid) |> 
  mutate(n = ifelse(is.na(n), 0, n)) |>
  mutate(Type = case_when(str_detect(variable, "eps") ~ "Extinction",
                          str_detect(variable, "gamma") ~ "Colonization",
                          str_detect(variable, "beta") ~ "Occupancy"
  )) |>
  mutate(variable = gsub("[[:alpha:]]", "", variable)) |> 
  mutate(PredPretty = case_when(Type == "Occupancy" & variable == 0 ~ "Intercept",
                                Type == "Occupancy" & variable == 1 ~ "Annual Tmax",
                                Type == "Occupancy" & variable == 2 ~ "Annual Tmax^2",
                                Type == "Occupancy" & variable == 3 ~ "Annual Ppt",
                                Type == "Occupancy" & variable == 4 ~ "Annual Ppt^2",
                                Type == "Occupancy" & variable == 5 ~ "Canopy Cover",
                                Type == "Occupancy" & variable == 6 ~ "Canopy Cover^2",
                                Type == "Colonization" & variable == 0 ~ "Intercept",
                                Type == "Colonization" & variable == 1 ~ "Tmax Anom",
                                Type == "Colonization" & variable == 2 ~ "Ppt Anom",
                                Type == "Colonization" & variable == 3 ~ "HSF",
                                Type == "Colonization" & variable == 4 ~ "Tmax Anom x HSF",
                                Type == "Colonization" & variable == 5 ~ "Ppt Anom x HSF",
                                Type == "Colonization" & variable == 6 ~ "Tmax Trend",
                                Type == "Colonization" & variable == 7 ~ "Tmin Trend",
                                Type == "Colonization" & variable == 8 ~ "Ppt Trend",
                                Type == "Colonization" & variable == 9 ~ "Tmax Trend x HSF",
                                Type == "Colonization" & variable == 10 ~ "Tmin Trend x HSF",
                                Type == "Colonization" & variable == 11 ~ "Ppt Trend x HSF",
                                Type == "Extinction" & variable == 0 ~ "Intercept",
                                Type == "Extinction" & variable == 1 ~ "Tmax Anom",
                                Type == "Extinction" & variable == 2 ~ "Ppt Anom",
                                Type == "Extinction" & variable == 3 ~ "HSF",
                                Type == "Extinction" & variable == 4 ~ "Tmax Anom x HSF",
                                Type == "Extinction" & variable == 5 ~ "Ppt Anom x HSF",
                                Type == "Extinction" & variable == 6 ~ "Tmax Trend",
                                Type == "Extinction" & variable == 7 ~ "Tmin Trend",
                                Type == "Extinction" & variable == 8 ~ "Ppt Trend",
                                Type == "Extinction" & variable == 9 ~ "Tmax Trend x HSF",
                                Type == "Extinction" & variable == 10 ~ "Tmin Trend x HSF",
                                Type == "Extinction" & variable == 11 ~ "Ppt Trend x HSF",
  )) |> 
  mutate(PredPretty = factor(PredPretty,
                             levels = rev(c("Intercept",
                                            "Tmax Anom",
                                            "Ppt Anom",
                                            "Tmax Trend",
                                            "Tmin Trend",
                                            "Ppt Trend",
                                            "HSF",
                                            "Tmax Anom x HSF",
                                            "Ppt Anom x HSF",
                                            "Tmax Trend x HSF",
                                            "Tmin Trend x HSF",
                                            "Ppt Trend x HSF",
                                            "Annual Tmax",
                                            "Annual Tmax^2",
                                            "Annual Ppt",
                                            "Annual Ppt^2",
                                            "Canopy Cover",
                                            "Canopy Cover^2")))) |> 
  filter(PredPretty != "Intercept") |> 
  mutate(SigPretty = factor(case_when(Sig == "Negative" ~ "-",
                                      Sig == "Positive" ~ "+",
                                      TRUE ~ "NS"),
                            levels = c("-", "NS", "+"))) |> 
  group_by(Type, PredPretty) |> 
  mutate(Pct = n/sum(n)) |> 
  ggplot(aes(x = SigPretty, y = PredPretty, fill = Pct)) +
  geom_tile(color = "white") + 
  scale_y_discrete(labels = function(x) {
    ifelse(grepl("\\^2", x), 
           parse(text = gsub(" ", "~", x)),  # Only parse labels with ^2
           x)  # Return original text for others
  }) +
  scale_fill_gradient(low = "#fff5f0", 
                      high = "#228B22",  #"#67000d", 
                      na.value = "grey90",
                      guide = guide_colorbar(barwidth = 15),
                      labels = scales::label_percent()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = unit(1, "lines"),
    aspect.ratio = 3.5,
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5)
  ) +
  ylab("") + 
  xlab("") +
  labs(fill = "% of species") + 
  facet_wrap(~Type, scale = "free_y", labeller = as_labeller(facet_names, default = label_parsed))

ggsave(eff.sum, filename = here("Figures/Effect_Summary_SpringClimHSF.png"),
       height = 8, width = 8, dpi = 600)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Additive Interactions
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ***********************************************************
##
## Section Notes: Logistic regression will estimate the interaction
## terms on the mulitplicative scale. We might want estimates
## for both multiplicative and additive scales
##
## Formula for RERI
## RERI = RR11 - R10 - RR01 + 1
## RR11 = exp(B1 + B2 + B3)
## RR10 = exp(B1)
## RR01 = exp(B2)
## ***********************************************************

glimpse(param_compile)



reri_fun <- function(b1, b2, b3) {
  RR11 <- exp(b1+b2+b3)
  RR10 <- exp(b1)
  RR01 <- exp(b2)
  RERI <- RR11 - RR10 - RR01 + 1
  return(RERI)
}

# Synergy Index (SI)
synergy_index <- function(b_fire, b_climate, b_interaction) {
  RR11 <- exp(b_fire + b_climate + b_interaction)
  RR10 <- exp(b_fire)
  RR01 <- exp(b_climate)
  
  # Ratio of observed to expected under additivity
  SI <- RR11 / (RR10 + RR01 - 1)
  return(SI)
}

# Calculate RERI using posterior samples
reri_from_samples <- samples_df |>
  filter(str_detect(variable, "^eps1\\[|^eps3\\[|^eps4\\[")) |>
  separate(col = variable, into = c("param", "sp.ind"), sep = "\\[") |>
  mutate(sp.ind = as.numeric(gsub("\\]", "", sp.ind))) |>
  pivot_wider(names_from = param, values_from = value) |>
  mutate(RERI = exp(eps1 + eps3 + eps4) - exp(eps1) - exp(eps3) + 1) |>
  group_by(sp.ind) |>
  summarise(
    mean = mean(RERI),
    median = median(RERI),
    lowci = quantile(RERI, 0.055),
    hici = quantile(RERI, 0.945),
    lower_80 = quantile(RERI, 0.1),
    upper_80 = quantile(RERI, 0.9),
    .groups = 'drop'
  ) |>
  mutate(variable = "RERI") |>
  select(variable, sp.ind, mean, median, lowci, hici, lower_80, upper_80) |> 
  left_join(sp.id |> mutate(sp.ind = as.numeric(sp.ind)))

si_from_samples <- samples_df |>
  filter(str_detect(variable, "^eps1\\[|^eps3\\[|^eps4\\[")) |>
  separate(col = variable, into = c("param", "sp.ind"), sep = "\\[") |>
  mutate(sp.ind = as.numeric(gsub("\\]", "", sp.ind))) |>
  pivot_wider(names_from = param, values_from = value) |>
  mutate(SI = exp(eps1 + eps3 + eps4) / (exp(eps1) + exp(eps3) - 1)) |>
  group_by(sp.ind) |>
  summarise(
    mean = mean(SI),
    median = median(SI),
    lowci = quantile(SI, 0.055),
    hici = quantile(SI, 0.945),
    lower_80 = quantile(SI, 0.1),
    upper_80 = quantile(SI, 0.9),
    .groups = 'drop'
  ) |>
  mutate(variable = "SI") |>
  select(variable, sp.ind, mean, median, lowci, hici, lower_80, upper_80) |> 
  left_join(sp.id |> mutate(sp.ind = as.numeric(sp.ind)))

ggplot(reri_from_samples) +
  geom_point(aes(x = species, y = mean)) +
  geom_errorbar(aes(ymax = hici, ymin = lowci, x = sp.ind)) +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggplot(si_from_samples) +
  geom_point(aes(x = species, y = mean)) +
  geom_errorbar(aes(ymax = hici, ymin = lowci, x = sp.ind)) +
  geom_hline(aes(yintercept = 1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Create comparison plot
multiplicative_vs_additive <- param_compile |>
  bind_rows(reri_from_samples |> select(-species)) |> 
  bind_rows(si_from_samples |> select(-species)) |> 
  left_join(sp.id |> mutate(sp.ind = as.numeric(sp.ind))) |> 
  filter(variable %in% c("eps4", "RERI", "SI")) |>
  select(variable, species, mean, lowci, hici) |>
  pivot_wider(names_from = variable, values_from = c(mean, lowci, hici))

# Plot both scales
ggplot(multiplicative_vs_additive, aes(x = mean_eps4, y = mean_SI, label = species)) +
  geom_point() +
  geom_errorbar(aes(ymin = lowci_SI, ymax = hici_SI), 
                width = 0, alpha = 0.5, color = "gray50") +
  geom_errorbar(aes(xmin = lowci_eps4, xmax = hici_eps4), 
                 height = 0, alpha = 0.5, color = "gray50") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel() +
  labs(x = "Multiplicative Interaction", 
       y = "Additive Interaction") +
  theme_bw()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Johnson-Neyman Figures
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test function for a single species
test_jnk_single <- function(mcmc_samples, species_id = 1, interaction_type = "tmax_anom", var_type = "temp", colext = "col",
                            species_key = specs) {
  
  cat("Testing JNK analysis for species", specs[species_id], "\n")
  
  # Define which parameters we need based on interaction type
  if(colext == "col"){
    if(interaction_type == "tmax_anom") {
      params_needed <- c("gamma3", "gamma1", "gamma4")  # fire, temp_trend, fire:temp
      col_names <- c("fire", "tmax_anom", "fire:tmax_anom")
    } else if (interaction_type == "precip_anom") {
      params_needed <- c("gamma3", "gamma2", "gamma5")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "precip_anom", "fire:precip_anom")
    } else if (interaction_type == "tmax_trend") {
      params_needed <- c("gamma3", "gamma6", "gamma9")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "tmax_trend", "fire:tmax_trend")
    } else if (interaction_type == "tmin_trend") {
      params_needed <- c("gamma3", "gamma7", "gamma10")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "tmin_trend", "fire:tmin_trend")
    } else if (interaction_type == "precip_trend") {
      params_needed <- c("gamma3", "gamma8", "gamma11")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "precip_trend", "fire:precip_trend")
    }
  } else {
    if(interaction_type == "tmax_anom") {
      params_needed <- c("eps3", "eps1", "eps4")  # fire, temp_trend, fire:temp
      col_names <- c("fire", "tmax_anom", "fire:tmax_anom")
    } else if (interaction_type == "precip_anom") {
      params_needed <- c("eps3", "eps2", "eps5")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "precip_anom", "fire:precip_anom")
    } else if (interaction_type == "tmax_trend") {
      params_needed <- c("eps3", "eps6", "eps9")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "tmax_trend", "fire:tmax_trend")
    } else if (interaction_type == "tmin_trend") {
      params_needed <- c("eps3", "eps7", "eps10")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "tmin_trend", "fire:tmin_trend")
    } else if (interaction_type == "precip_trend") {
      params_needed <- c("eps3", "eps8", "eps11")  # fire, precip_trend, fire:precip
      col_names <- c("fire", "precip_trend", "fire:precip_trend")
    }
  }
  
  cat("Looking for parameters:", paste(params_needed, collapse = ", "), "\n")
  
  # Extract samples for this species
  species_data <- mcmc_samples |>
    filter(str_detect(variable, paste0("\\[", species_id, "\\]$"))) |>
    filter(str_detect(variable, paste0("^(", paste(params_needed, collapse = "|"), ")"))) |>
    separate(variable, into = c("param", "sp"), sep = "\\[") |>
    mutate(sp = as.numeric(gsub("\\]", "", sp)))
  
  cat("Found", nrow(species_data), "rows of data\n")
  cat("Parameters found:", unique(species_data$param), "\n")
  
  # Check if we have all needed parameters
  params_found <- unique(species_data$param)
  missing_params <- setdiff(params_needed, params_found)
  if(length(missing_params) > 0) {
    stop("Missing parameters: ", paste(missing_params, collapse = ", "))
  }
  
  # Convert to matrix format
  species_matrix <- species_data |>
    group_by(param) |>
    mutate(iteration = row_number()) |>
    ungroup() |>
    select(param, iteration, value) |>
    pivot_wider(names_from = param, values_from = value) |>
    select(-iteration) |>
    as.matrix()
  
  # Reorder columns to match expected order
  species_matrix <- species_matrix[, params_needed]
  colnames(species_matrix) <- col_names
  
  cat("Matrix dimensions:", dim(species_matrix), "\n")
  cat("Column names:", colnames(species_matrix), "\n")
  cat("First few rows:\n")
  print(head(species_matrix, 3))
  
  # Define value ranges (conservative ranges first)
  fire_range <- seq(-2, 2, 1)
  climate_range <- seq(-2, 2, 1)
  
  cat("Testing JNK_bayes function...\n")
  
  # Run JNK_bayes
  if(var_type == "temp") {
    jnk_result <- JNK_bayes(
      x = species_matrix,
      theta_1 = "fire",
      theta_2 = interaction_type,
      theta_int_12 = paste0("fire:", interaction_type),
      theta_1_vals = fire_range,
      theta_2_vals = climate_range,
      thresholds = c(0.055, 0.945),
      save = FALSE
    )
  } else {
    jnk_result <- JNK_bayes(
      x = species_matrix,
      theta_1 = "fire", 
      theta_2 = "precip",
      theta_int_12 = "fire:precip",
      theta_1_vals = fire_range,
      theta_2_vals = climate_range,
      thresholds = c(0.055, 0.945),
      noTitle = paste("Species", species_id, "- Fire × Precipitation"),
      save = FALSE
    )
  }
  
  cat("JNK analysis completed successfully!\n")
  return(jnk_result)
}

# Test with temperature interaction
out <- vector(mode = "list", length = length(specs))
names(out) <- specs
for(i in 1:length(specs)){
  test_result <- test_jnk_single(mcmc_samples = samples_df, 
                                 species_id = i, 
                                 interaction_type = "tmax_trend",
                                 var_type = "temp",
                                 colext = "ext")
  out[[i]] <- test_result
}
# View the plot
test_result$fire_plot
test_result$fire_table
get_plot <- function(x){
  return(x$fire_plot)
}
get_table <- function(x){
  return(x$fire_table)
}

wrap_plots(lapply(out[1:9], get_plot), ncol = 3, nrow = 3)
wrap_plots(lapply(out[10:18], get_plot), ncol = 3, nrow = 3)
wrap_plots(lapply(out[19:27], get_plot), ncol = 3, nrow = 3)
wrap_plots(lapply(out[28:36], get_plot), ncol = 3, nrow = 3)
wrap_plots(lapply(out[37:45], get_plot), ncol = 3, nrow = 3)
wrap_plots(lapply(out[46:51], get_plot), ncol = 3, nrow = 3)

lapply(out, get_table)
