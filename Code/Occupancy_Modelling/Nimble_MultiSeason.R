## -------------------------------------------------------------
##
## Script name: Nimble Model for single species, multi-season colex
##
## Script purpose: Dynamic occ model
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
##
##
## -------------------------------------------------------------

## Defaults
options(scipen = 10, digits = 10)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(here)
library(nimble)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Model Code
##
## -------------------------------------------------------------

dynOcc <- nimbleCode({
  
  ## Priors
  ## Occupancy coefficients
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  #beta2 ~ dnorm(0, 0.01)
  
  ## Detection coefficients
  alpha0 ~ dnorm(0, 0.01)
  alpha1 ~ dnorm(0, 0.01)
  alpha2 ~ dnorm(0, 0.01)
  
  ## Persistence coefficients
  phi0 ~ dnorm(0, 0.01)
  phi1 ~ dnorm(0, 0.01)
  #phi2 ~ dnorm(0, 0.01)
  #phi3 ~ dnorm(0, 0.01)
  
  ## Colonization coefficients
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)
  #gamma2 ~ dnorm(0, 0.01)
  #gamma3 ~ dnorm(0, 0.01)
  
  ## Likelihood
  ## Ecological submodel
  for(i in 1:nsites){
    
    ## Predictors for intial occ
    logit(psi1[i]) <- beta0 + beta1 * clim[i, 1] #+ beta2 * fire[i, 1]
    
    ## Initial occ
    z[i, 1] ~ dbern(psi1[i])
    
    ## State process model (colex)
    for(t in 2:nyears){
      
      ## Colonization
      logit(gamma[i, t-1]) <- gamma0 + gamma1 * clim[i, t-1] #+ gamma2 * fire[i, t-1] + gamma3 * anom[i, t-1] * fire[i, t-1]
      
      ## Persistence
      logit(phi[i, t-1]) <- phi0 + phi1 * clim[i, t-1] #+ phi2 * fire[i, t-1] + phi3 * anom[i, t-1] * fire[i, t-1]
      
      ## Latent state
      z[i,t] ~ dbern(z[i, t-1] * phi[i, t-1] + (1-z[i,t-1]) * gamma[i, t-1])
    }
  }
  
  ## Observation submodel
  for(i in 1:nsites){
    for(j in 1:nreps){
      for(t in 1:nyears){
        logit(p[i,j,t]) <- alpha0 + alpha1 * eff.hrs[i,j,t] + alpha2 * eff.jday[i,j,t]
        y[i,j,t] ~ dbern(z[i,t] * p[i,j,t])
      }
    }
  }
  
  ## Derived parameters
  for(i in 1:nsites){
    for(t in 1:(nyears-1)){
      eps[i, t] <- 1 - phi[i, t]
    }
  }
  
  ## Mean occupancy
  psi.fs[1] <- mean(psi1[1:nsites])
  for(t in 2:nyears){
    psi.fs[t] <- mean(z[1:nsites, t])
  }
  
})

## -------------------------------------------------------------
##
## End Section: Model Code
##
## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Model data and inits 
##
## -------------------------------------------------------------
zst <- apply(win.data$y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() {
  list(z = zst,
       beta0 = rnorm(1),
       beta1 = rnorm(1),
       #beta2 = rnorm(1),
       alpha0 = rnorm(1),
       alpha1 = rnorm(1),
       alpha2 = rnorm(1),
       gamma0 = rnorm(1),
       gamma1 = rnorm(1),
       #gamma2 = rnorm(1),
       phi0 = rnorm(1),
       phi1 = rnorm(1))
} 

# Parameters monitored
params <- c("alpha0", 
            "alpha1",
            "alpha2",
            "beta0", 
            "beta1", 
            "gamma0",
            "gamma1",
            "phi0",
            "phi1",
            "eps",
            "psi.fs")

## Load data
bdata <- readRDS(here("Data/Occ_Data/LABU_Multi_Test.RDS"))  


## MCMC settings
nburnin = 2000
niter = 20000
nthin = 20

## Run the sampler
samples <- nimbleMCMC(
  code = dynOcc,
  constants = bdata,
  inits = inits,
  monitors = params,
  niter = niter,
  nburnin = nburnin,
  thin = nthin)
