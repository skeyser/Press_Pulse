## -------------------------------------------------------------
##
## Script name: Dynamic community model in NIMBLE
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

renv::load()

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: DCM Traits with continuous fire
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DCMTraitsContFire <- nimbleCode({
  
  ##----------------------
  ## Species-level priors
  ##----------------------
  for(k in 1:nspec){
    ## Occupancy coefficients
    beta0[k] ~ dnorm(mu.beta0, tau.beta0)
    beta1[k] ~ dnorm(mu.beta1, tau.beta1)
    beta2[k] ~ dnorm(mu.beta2, tau.beta2)
    beta3[k] ~ dnorm(mu.beta3, tau.beta3)
    beta4[k] ~ dnorm(mu.beta4, tau.beta4)
    beta5[k] ~ dnorm(mu.beta5, tau.beta5)
    beta6[k] ~ dnorm(mu.beta6, tau.beta6)
    
    ## Detection coefficients
    alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0)
    alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)
    alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2)
    alpha3[k] ~ dnorm(mu.alpha3, tau.alpha3)
    alpha4[k] ~ dnorm(mu.alpha4, tau.alpha4)
    
    ## Colonization coefficients
    gamma0[k] ~ dnorm(mu.gamma0[k], tau.gamma0)
    gamma1[k] ~ dnorm(mu.gamma1[k], tau.gamma1)
    gamma2[k] ~ dnorm(mu.gamma2[k], tau.gamma2)
    gamma3[k] ~ dnorm(mu.gamma3[k], tau.gamma3)
    gamma4[k] ~ dnorm(mu.gamma4[k], tau.gamma4)
    gamma5[k] ~ dnorm(mu.gamma5[k], tau.gamma5)
    gamma6[k] ~ dnorm(mu.gamma6[k], tau.gamma6)
    gamma7[k] ~ dnorm(mu.gamma7[k], tau.gamma7)
    gamma8[k] ~ dnorm(mu.gamma8[k], tau.gamma8)
    gamma9[k] ~ dnorm(mu.gamma9[k], tau.gamma9)
    gamma10[k] ~ dnorm(mu.gamma10[k], tau.gamma10)
    gamma11[k] ~ dnorm(mu.gamma11[k], tau.gamma11)
    
    ## Extinction coefficients
    eps0[k] ~ dnorm(mu.eps0[k], tau.eps0)
    eps1[k] ~ dnorm(mu.eps1[k], tau.eps1)
    eps2[k] ~ dnorm(mu.eps2[k], tau.eps2)
    eps3[k] ~ dnorm(mu.eps3[k], tau.eps3)
    eps4[k] ~ dnorm(mu.eps4[k], tau.eps4)
    eps5[k] ~ dnorm(mu.eps5[k], tau.eps5)
    eps6[k] ~ dnorm(mu.eps6[k], tau.eps6)
    eps7[k] ~ dnorm(mu.eps7[k], tau.eps7)
    eps8[k] ~ dnorm(mu.eps8[k], tau.eps8)
    eps9[k] ~ dnorm(mu.eps9[k], tau.eps9)
    eps10[k] ~ dnorm(mu.eps10[k], tau.eps10)
    eps11[k] ~ dnorm(mu.eps11[k], tau.eps11)
    
    ## Trait-based models of species effects
    ## Colonization
    mu.gamma0[k] <- delta0.gamma0 + 
      delta1.gamma0 * mass[k] + 
      delta2.gamma0 * hwi[k]
    
    mu.gamma1[k] <- delta0.gamma1 + 
      delta1.gamma1 * mass[k] + 
      delta2.gamma1 * hwi[k]
    
    mu.gamma2[k] <- delta0.gamma2 + 
      delta1.gamma2 * mass[k] + 
      delta2.gamma2 * hwi[k]
    
    mu.gamma3[k] <- delta0.gamma3 + 
      delta1.gamma3 * mass[k] + 
      delta2.gamma3 * hwi[k]
    
    mu.gamma4[k] <- delta0.gamma4 + 
      delta1.gamma4 * mass[k] + 
      delta2.gamma4 * hwi[k]
    
    mu.gamma5[k] <- delta0.gamma5 + 
      delta1.gamma5 * mass[k] + 
      delta2.gamma5 * hwi[k]
    
    mu.gamma6[k] <- delta0.gamma6 + 
      delta1.gamma6 * mass[k] + 
      delta2.gamma6 * hwi[k]
    
    mu.gamma7[k] <- delta0.gamma7 + 
      delta1.gamma7 * mass[k] + 
      delta2.gamma7 * hwi[k]
    
    mu.gamma8[k] <- delta0.gamma8 + 
      delta1.gamma8 * mass[k] + 
      delta2.gamma8 * hwi[k]
    
    mu.gamma9[k] <- delta0.gamma9 + 
      delta1.gamma9 * mass[k] + 
      delta2.gamma9 * hwi[k]
    
    mu.gamma10[k] <- delta0.gamma10 + 
      delta1.gamma10 * mass[k] + 
      delta2.gamma10 * hwi[k]
    
    mu.gamma11[k] <- delta0.gamma11 + 
      delta1.gamma11 * mass[k] + 
      delta2.gamma11 * hwi[k]
    
    ## Extinction
    mu.eps0[k] <- delta0.eps0 + 
      delta1.eps0 * mass[k] + 
      delta2.eps0 * hwi[k]
    
    mu.eps1[k] <- delta0.eps1 + 
      delta1.eps1 * mass[k] + 
      delta2.eps1 * hwi[k]
    
    mu.eps2[k] <- delta0.eps2 + 
      delta1.eps2 * mass[k] + 
      delta2.eps2 * hwi[k]
    
    mu.eps3[k] <- delta0.eps3 + 
      delta1.eps3 * mass[k] + 
      delta2.eps3 * hwi[k]
    
    mu.eps4[k] <- delta0.eps4 + 
      delta1.eps4 * mass[k] + 
      delta2.eps4 * hwi[k]
    
    mu.eps5[k] <- delta0.eps5 + 
      delta1.eps5 * mass[k] + 
      delta2.eps5 * hwi[k]
    
    mu.eps6[k] <- delta0.eps6 + 
      delta1.eps6 * mass[k] + 
      delta2.eps6 * hwi[k]
    
    mu.eps7[k] <- delta0.eps7 + 
      delta1.eps7 * mass[k] + 
      delta2.eps7 * hwi[k]
    
    mu.eps8[k] <- delta0.eps8 + 
      delta1.eps8 * mass[k] + 
      delta2.eps8 * hwi[k]
    
    mu.eps9[k] <- delta0.eps9 + 
      delta1.eps9 * mass[k] + 
      delta2.eps9 * hwi[k]
    
    mu.eps10[k] <- delta0.eps10 + 
      delta1.eps10 * mass[k] + 
      delta2.eps10 * hwi[k]
    
    mu.eps11[k] <- delta0.eps11 + 
      delta1.eps11 * mass[k] + 
      delta2.eps11 * hwi[k]
  }
  
  ## Cell-level raneff
  for(c in 1:n_cells){
    cell_det[c] ~ dnorm(0, tau.cell_det)
  }
  sd.cell_det ~ dunif(0, 2)
  tau.cell_det <- pow(sd.cell_det, -2)
  
  ## Year ranef
  # Year effect priors:
  for(k in 1:nspec) {
    for(t in 1:nyears) {
      alpha_year[k,t] ~ dnorm(0, tau.alpha_year)
    }
  }
  sd.alpha_year ~ dunif(0, 2)
  tau.alpha_year ~ pow(sd.alpha_year, 2)
  
  ## False-positive Priors
  ## Edit this for a beta distribution
  ## can estimate the parameters per species
  ## using validation data to try to make this
  ## a better representation of the data
  for(k in 1:nspec){
    #alphaFP[k] ~ dnorm(-2.944, 4)
    alphaFP[k] ~ dbeta(shape1 = 5, shape2 = 195)
  }
  
  ##----------------------
  ## Hyperpriors
  ##----------------------
  ## Occupancy hyperpriors
  mu.beta0 ~ dnorm(0, 0.01)
  sd.beta0 ~ dunif(0, 2)
  tau.beta0 <- pow(sd.beta0, -2)
  
  mu.beta1 ~ dnorm(0, 0.01)
  sd.beta1 ~ dunif(0, 2)
  tau.beta1 <- pow(sd.beta1, -2)
  
  mu.beta2 ~ dnorm(0, 0.01)
  sd.beta2 ~ dunif(0, 2)
  tau.beta2 <- pow(sd.beta2, -2)
  
  mu.beta3 ~ dnorm(0, 0.01)
  sd.beta3 ~ dunif(0, 2)
  tau.beta3 <- pow(sd.beta3, -2)
  
  mu.beta4 ~ dnorm(0, 0.01)
  sd.beta4 ~ dunif(0, 2)
  tau.beta4 <- pow(sd.beta4, -2)
  
  mu.beta5 ~ dnorm(0, 0.01)
  sd.beta5 ~ dunif(0, 2)
  tau.beta5 <- pow(sd.beta5, -2)
  
  mu.beta6 ~ dnorm(0, 0.01)
  sd.beta6 ~ dunif(0, 2)
  tau.beta6 <- pow(sd.beta6, -2)
  
  mu.beta7 ~ dnorm(0, 0.01)
  sd.beta7 ~ dunif(0, 2)
  tau.beta7 <- pow(sd.beta7, -2)
  
  ## Detection hyperpriors
  mu.alpha0 ~ dnorm(0, 0.01)
  sd.alpha0 ~ dunif(0, 2)
  tau.alpha0 <- pow(sd.alpha0, -2)
  
  mu.alpha1 ~ dnorm(0, 0.01)
  sd.alpha1 ~ dunif(0, 2)
  tau.alpha1 <- pow(sd.alpha1, -2)
  
  mu.alpha2 ~ dnorm(0, 0.01)
  sd.alpha2 ~ dunif(0, 2)
  tau.alpha2 <- pow(sd.alpha2, -2)
  
  mu.alpha3 ~ dnorm(0, 0.01)
  sd.alpha3 ~ dunif(0, 2)
  tau.alpha3 <- pow(sd.alpha3, -2)
  
  mu.alpha4 ~ dnorm(0, 0.01)
  sd.alpha4 ~ dunif(0, 2)
  tau.alpha4 <- pow(sd.alpha4, -2)
  
  ## Extinction hyperpriors
  #mu.eps0 ~ dnorm(0, 0.01)
  sd.eps0 ~ dunif(0, 2)
  tau.eps0 <- pow(sd.eps0, -2)
  
  #mu.eps1 ~ dnorm(0, 0.01)
  sd.eps1 ~ dunif(0, 2)
  tau.eps1 <- pow(sd.eps1, -2)
  
  #mu.eps2 ~ dnorm(0, 0.01)
  sd.eps2 ~ dunif(0, 2)
  tau.eps2 <- pow(sd.eps2, -2)
  
  #mu.eps3 ~ dnorm(0, 0.01)
  sd.eps3 ~ dunif(0, 2)
  tau.eps3 <- pow(sd.eps3, -2)
  
  #mu.eps4 ~ dnorm(0, 0.01)
  sd.eps4 ~ dunif(0, 2)
  tau.eps4 <- pow(sd.eps4, -2)
  
  #mu.eps5 ~ dnorm(0, 0.01)
  sd.eps5 ~ dunif(0, 2)
  tau.eps5 <- pow(sd.eps5, -2)
  
  #mu.eps6 ~ dnorm(0, 0.01)
  sd.eps6 ~ dunif(0, 2)
  tau.eps6 <- pow(sd.eps6, -2)
  
  #mu.eps7 ~ dnorm(0, 0.01)
  sd.eps7 ~ dunif(0, 2)
  tau.eps7 <- pow(sd.eps7, -2)
  
  #mu.eps8 ~ dnorm(0, 0.01)
  sd.eps8 ~ dunif(0, 2)
  tau.eps8 <- pow(sd.eps8, -2)
  
  #mu.eps9 ~ dnorm(0, 0.01)
  sd.eps9 ~ dunif(0, 2)
  tau.eps9 <- pow(sd.eps9, -2)
  
  #mu.eps10 ~ dnorm(0, 0.01)
  sd.eps10 ~ dunif(0, 2)
  tau.eps10 <- pow(sd.eps10, -2)
  
  #mu.eps11 ~ dnorm(0, 0.01)
  sd.eps11 ~ dunif(0, 2)
  tau.eps11 <- pow(sd.eps11, -2)
  
  ## Colonization hyperpriors
  #mu.gamma0 ~ dnorm(0, 0.01)
  sd.gamma0 ~ dunif(0, 2)
  tau.gamma0 <- pow(sd.gamma0, -2)
  
  #mu.gamma1 ~ dnorm(0, 0.01)
  sd.gamma1 ~ dunif(0, 2)
  tau.gamma1 <- pow(sd.gamma1, -2)
  
  #mu.gamma2 ~ dnorm(0, 0.01)
  sd.gamma2 ~ dunif(0, 2)
  tau.gamma2 <- pow(sd.gamma2, -2)
  
  ## Replaced below with mean trait response
  #mu.gamma3 ~ dnorm(0, 0.01)
  sd.gamma3 ~ dunif(0, 2)
  tau.gamma3 <- pow(sd.gamma3, -2)
  
  #mu.gamma4 ~ dnorm(0, 0.01)
  sd.gamma4 ~ dunif(0, 2)
  tau.gamma4 <- pow(sd.gamma4, -2)
  
  #mu.gamma5 ~ dnorm(0, 0.01)
  sd.gamma5 ~ dunif(0, 2)
  tau.gamma5 <- pow(sd.gamma5, -2)
  
  #mu.gamma6 ~ dnorm(0, 0.01)
  sd.gamma6 ~ dunif(0, 2)
  tau.gamma6 <- pow(sd.gamma6, -2)
  
  #mu.gamma7 ~ dnorm(0, 0.01)
  sd.gamma7 ~ dunif(0, 2)
  tau.gamma7 <- pow(sd.gamma7, -2)
  
  #mu.gamma8 ~ dnorm(0, 0.01)
  sd.gamma8 ~ dunif(0, 2)
  tau.gamma8 <- pow(sd.gamma8, -2)
  
  #mu.gamma9 ~ dnorm(0, 0.01)
  sd.gamma9 ~ dunif(0, 2)
  tau.gamma9 <- pow(sd.gamma9, -2)
  
  #mu.gamma10 ~ dnorm(0, 0.01)
  sd.gamma10 ~ dunif(0, 2)
  tau.gamma10 <- pow(sd.gamma10, -2)
  
  #mu.gamma11 ~ dnorm(0, 0.01)
  sd.gamma11 ~ dunif(0, 2)
  tau.gamma11 <- pow(sd.gamma11, -2)
  
  ## Trait-effect priors
  ## Colonization
  delta0.gamma0 ~ dnorm(0, 0.01)
  delta1.gamma0 ~ dnorm(0, 0.01)
  delta2.gamma0 ~ dnorm(0, 0.01)
  delta0.gamma1 ~ dnorm(0, 0.01)
  delta1.gamma1 ~ dnorm(0, 0.01)
  delta2.gamma1 ~ dnorm(0, 0.01)
  delta0.gamma2 ~ dnorm(0, 0.01)
  delta1.gamma2 ~ dnorm(0, 0.01)
  delta2.gamma2 ~ dnorm(0, 0.01)
  delta0.gamma3 ~ dnorm(0, 0.01)
  delta1.gamma3 ~ dnorm(0, 0.01)
  delta2.gamma3 ~ dnorm(0, 0.01)
  delta0.gamma4 ~ dnorm(0, 0.01)
  delta1.gamma4 ~ dnorm(0, 0.01)
  delta2.gamma4 ~ dnorm(0, 0.01)
  delta0.gamma5 ~ dnorm(0, 0.01)
  delta1.gamma5 ~ dnorm(0, 0.01)
  delta2.gamma5 ~ dnorm(0, 0.01)
  delta0.gamma6 ~ dnorm(0, 0.01)
  delta1.gamma6 ~ dnorm(0, 0.01)
  delta2.gamma6 ~ dnorm(0, 0.01)
  delta0.gamma7 ~ dnorm(0, 0.01)
  delta1.gamma7 ~ dnorm(0, 0.01)
  delta2.gamma7 ~ dnorm(0, 0.01)
  delta0.gamma8 ~ dnorm(0, 0.01)
  delta1.gamma8 ~ dnorm(0, 0.01)
  delta2.gamma8 ~ dnorm(0, 0.01)
  delta0.gamma9 ~ dnorm(0, 0.01)
  delta1.gamma9 ~ dnorm(0, 0.01)
  delta2.gamma9 ~ dnorm(0, 0.01)
  delta0.gamma10 ~ dnorm(0, 0.01)
  delta1.gamma10 ~ dnorm(0, 0.01)
  delta2.gamma10 ~ dnorm(0, 0.01)
  delta0.gamma11 ~ dnorm(0, 0.01)
  delta1.gamma11 ~ dnorm(0, 0.01)
  delta2.gamma11 ~ dnorm(0, 0.01)
  
  ## Extinction
  delta0.eps0 ~ dnorm(0, 0.01)
  delta1.eps0 ~ dnorm(0, 0.01)
  delta2.eps0 ~ dnorm(0, 0.01)
  delta0.eps1 ~ dnorm(0, 0.01)
  delta1.eps1 ~ dnorm(0, 0.01)
  delta2.eps1 ~ dnorm(0, 0.01)
  delta0.eps2 ~ dnorm(0, 0.01)
  delta1.eps2 ~ dnorm(0, 0.01)
  delta2.eps2 ~ dnorm(0, 0.01)
  delta0.eps3 ~ dnorm(0, 0.01)
  delta1.eps3 ~ dnorm(0, 0.01)
  delta2.eps3 ~ dnorm(0, 0.01)
  delta0.eps4 ~ dnorm(0, 0.01)
  delta1.eps4 ~ dnorm(0, 0.01)
  delta2.eps4 ~ dnorm(0, 0.01)
  delta0.eps5 ~ dnorm(0, 0.01)
  delta1.eps5 ~ dnorm(0, 0.01)
  delta2.eps5 ~ dnorm(0, 0.01)
  delta0.eps6 ~ dnorm(0, 0.01)
  delta1.eps6 ~ dnorm(0, 0.01)
  delta2.eps6 ~ dnorm(0, 0.01)
  delta0.eps7 ~ dnorm(0, 0.01)
  delta1.eps7 ~ dnorm(0, 0.01)
  delta2.eps7 ~ dnorm(0, 0.01)
  delta0.eps8 ~ dnorm(0, 0.01)
  delta1.eps8 ~ dnorm(0, 0.01)
  delta2.eps8 ~ dnorm(0, 0.01)
  delta0.eps9 ~ dnorm(0, 0.01)
  delta1.eps9 ~ dnorm(0, 0.01)
  delta2.eps9 ~ dnorm(0, 0.01)
  delta0.eps10 ~ dnorm(0, 0.01)
  delta1.eps10 ~ dnorm(0, 0.01)
  delta2.eps10 ~ dnorm(0, 0.01)
  delta0.eps11 ~ dnorm(0, 0.01)
  delta1.eps11 ~ dnorm(0, 0.01)
  delta2.eps11 ~ dnorm(0, 0.01)
  
  ##----------------------
  ## Ecological State Process
  ##----------------------
  for(k in 1:nspec){
    for(i in 1:nsites){
      
      ## Initial occupancy
      logit(psi1[i,k]) <- beta0[k] + 
        beta1[k] * btemp[i] +
        beta2[k] * bprec[i] +
        beta3[k] * cc[i] +
        beta4[k] * cc2[i] +
        beta5[k] * trendt[i] +
        beta6[k] * trendtmin[i] +
        beta7[k] * trendp[i]
      
      z[i,1,k] ~ dbern(psi1[i,k])
      psi[i,1,k] <- psi1[i,k]
      
      ## State process over years
      for(t in 2:nyears){
        ## Colonization
        logit(gamma[i,t-1,k]) <- gamma0[k] +
          gamma1[k] * tanom[i,t-1] +
          gamma2[k] * panom[i,t-1] +
          gamma3[k] * hsf_pland[i, t-1] +
          gamma4[k] * fire_pattern[i,t-1] +
          gamma5[k] * tanom[i,t-1] * hsf_pland[i,t-1] +
          gamma6[k] * tanom[i,t-1] * fire_pattern[i,t-1] +
          gamma7[k] * panom[i,t-1] * hsf_pland[i,t-1] +
          gamma8[k] * panom[i,t-1] * fire_pattern[i,t-1] +
          gamma9[k] * trendt[i] + 
          gamma10[k] * trendtmin[i] +
          gamma11[k] * trendp[i]
        
        ## Persistence
        logit(eps[i,t-1,k]) <- eps0[k] +
          eps1[k] * tanom[i,t-1] +
          eps2[k] * panom[i,t-1] +
          eps3[k] * hsf_pland[i,t-1] +
          eps4[k] * fire_pattern[i,t-1] +
          eps5[k] * tanom[i,t-1] * hsf_pland[i,t-1] +
          eps6[k] * tanom[i,t-1] * fire_pattern[i,t-1] +
          eps7[k] * panom[i,t-1] * hsf_pland[i,t-1] +
          eps8[k] * panom[i,t-1] * fire_pattern[i,t-1] +
          eps9[k] * trendt[i] +
          eps10[k] * trendtmin[i] +
          eps11[k] * trendp[i]
        
        ## Latent state
        z[i,t,k] ~ dbern(z[i,t-1,k] * (1-eps[i,t-1,k]) +
                           (1-z[i,t-1,k]) * gamma[i,t-1,k])
        
        ## Derived psi
        psi[i,t,k] <- psi[i,t-1,k] * (1-eps[i,t-1,k]) +
          (1-psi[i,t-1,k]) * gamma[i,t-1,k]
      }
    }
  }
  
  ## Obs sub-model ragged
  for(v in 1:nObs) {
    for(k in 1:nspec) {
      logit(p_obs[v,k]) <- alpha0[k] +
        alpha1[k] * eff.hrs[v] +
        alpha2[k] * eff.jday[v] + 
        alpha3[k] * eff.jday2[v] +
        alpha4[k] * ele[site_obs[v]] +
        alpha_year[k, year_obs[v]] + 
        cell_det[cell_id[site_obs[v]]]
      
      y[v,k] ~ dbern(z[site_obs[v], year_obs[v], k] * p_obs[v,k] +
                       (1 - z[site_obs[v], year_obs[v], k]) * alphaFP[k])
    }
  }
  
  ##----------------------
  ## Derived parameters
  ##----------------------
  for(k in 1:nspec){
    for(i in 1:nsites){
      for(t in 1:(nyears-1)){
        phi[i,t,k] <- 1 - eps[i,t,k]
      }
    }
  }
  
  ## Mean occupancy per species
  for(k in 1:nspec){
    psi.fs[1,k] <- mean(psi1[1:nsites,k])
    for(t in 2:nyears){
      psi.fs[t,k] <- mean(z[1:nsites,t,k])
    }
  }
  
  ## Mean richness
  for(i in 1:nsites){
    for(t in 1:nyears){
      richness[i,t] <- sum(z[i,t,1:nspec])
    }
  }
  
  ## FP Rate
  # for(k in 1:nspec){
  #   fp_rate[k] <- ilogit(alphaFP[k])
  # }
  
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: DCM Continuous Fire No traits
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: DCM Traits with continuous fire
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DCMContFire <- nimbleCode({
  
  ##----------------------
  ## Species-level priors
  ##----------------------
  for(k in 1:nspec){
    ## Occupancy coefficients
    beta0[k] ~ dnorm(mu.beta0, tau.beta0)
    beta1[k] ~ dnorm(mu.beta1, tau.beta1)
    beta2[k] ~ dnorm(mu.beta2, tau.beta2)
    beta3[k] ~ dnorm(mu.beta3, tau.beta3)
    beta4[k] ~ dnorm(mu.beta4, tau.beta4)
    beta5[k] ~ dnorm(mu.beta5, tau.beta5)
    beta6[k] ~ dnorm(mu.beta6, tau.beta6)
    
    ## Detection coefficients
    alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0)
    alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)
    alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2)
    alpha3[k] ~ dnorm(mu.alpha3, tau.alpha3)
    alpha4[k] ~ dnorm(mu.alpha4, tau.alpha4)
    
    ## Colonization coefficients
    gamma0[k] ~ dnorm(mu.gamma0, tau.gamma0)
    gamma1[k] ~ dnorm(mu.gamma1, tau.gamma1)
    gamma2[k] ~ dnorm(mu.gamma2, tau.gamma2)
    gamma3[k] ~ dnorm(mu.gamma3, tau.gamma3)
    gamma4[k] ~ dnorm(mu.gamma4, tau.gamma4)
    gamma5[k] ~ dnorm(mu.gamma5, tau.gamma5)
    gamma6[k] ~ dnorm(mu.gamma6, tau.gamma6)
    gamma7[k] ~ dnorm(mu.gamma7, tau.gamma7)
    gamma8[k] ~ dnorm(mu.gamma8, tau.gamma8)
    gamma9[k] ~ dnorm(mu.gamma9, tau.gamma9)
    gamma10[k] ~ dnorm(mu.gamma10, tau.gamma10)
    gamma11[k] ~ dnorm(mu.gamma11, tau.gamma11)

    ## Extinction coefficients
    eps0[k] ~ dnorm(mu.eps0, tau.eps0)
    eps1[k] ~ dnorm(mu.eps1, tau.eps1)
    eps2[k] ~ dnorm(mu.eps2, tau.eps2)
    eps3[k] ~ dnorm(mu.eps3, tau.eps3)
    eps4[k] ~ dnorm(mu.eps4, tau.eps4)
    eps5[k] ~ dnorm(mu.eps5, tau.eps5)
    eps6[k] ~ dnorm(mu.eps6, tau.eps6)
    eps7[k] ~ dnorm(mu.eps7, tau.eps7)
    eps8[k] ~ dnorm(mu.eps8, tau.eps8)
    eps9[k] ~ dnorm(mu.eps9, tau.eps9)
    eps10[k] ~ dnorm(mu.eps10, tau.eps10)
    eps11[k] ~ dnorm(mu.eps11, tau.eps11)
  }
  
  ## Cell-level raneff
  for(c in 1:n_cells){
    cell_det[c] ~ dnorm(0, tau.cell_det)
  }
  sd.cell_det ~ dunif(0, 2)
  tau.cell_det <- pow(sd.cell_det, -2)
  
  ## Year ranef
  # Year effect priors:
  for(k in 1:nspec) {
    for(t in 1:nyears) {
      alpha_year[k,t] ~ dnorm(0, tau.alpha_year)
    }
  }
  sd.alpha_year ~ dunif(0, 2)
  tau.alpha_year <- pow(sd.alpha_year, -2)
  
  ## False-positive Priors
  ## Edit this for a beta distribution
  ## can estimate the parameters per species
  ## using validation data to try to make this
  ## a better representation of the data
  for(k in 1:nspec){
    #alphaFP[k] ~ dnorm(-2.944, 4)
    alphaFP[k] ~ dbeta(shape1 = 5, shape2 = 195)
  }
  
  ##----------------------
  ## Hyperpriors
  ##----------------------
  ## Occupancy hyperpriors
  mu.beta0 ~ dnorm(0, 0.01)
  sd.beta0 ~ dunif(0, 2)
  tau.beta0 <- pow(sd.beta0, -2)
  
  mu.beta1 ~ dnorm(0, 0.01)
  sd.beta1 ~ dunif(0, 2)
  tau.beta1 <- pow(sd.beta1, -2)
  
  mu.beta2 ~ dnorm(0, 0.01)
  sd.beta2 ~ dunif(0, 2)
  tau.beta2 <- pow(sd.beta2, -2)
  
  mu.beta3 ~ dnorm(0, 0.01)
  sd.beta3 ~ dunif(0, 2)
  tau.beta3 <- pow(sd.beta3, -2)
  
  mu.beta4 ~ dnorm(0, 0.01)
  sd.beta4 ~ dunif(0, 2)
  tau.beta4 <- pow(sd.beta4, -2)
  
  mu.beta5 ~ dnorm(0, 0.01)
  sd.beta5 ~ dunif(0, 2)
  tau.beta5 <- pow(sd.beta5, -2)
  
  mu.beta6 ~ dnorm(0, 0.01)
  sd.beta6 ~ dunif(0, 2)
  tau.beta6 <- pow(sd.beta6, -2)
  
  ## Detection hyperpriors
  mu.alpha0 ~ dnorm(0, 0.01)
  sd.alpha0 ~ dunif(0, 2)
  tau.alpha0 <- pow(sd.alpha0, -2)
  
  mu.alpha1 ~ dnorm(0, 0.01)
  sd.alpha1 ~ dunif(0, 2)
  tau.alpha1 <- pow(sd.alpha1, -2)
  
  mu.alpha2 ~ dnorm(0, 0.01)
  sd.alpha2 ~ dunif(0, 2)
  tau.alpha2 <- pow(sd.alpha2, -2)
  
  mu.alpha3 ~ dnorm(0, 0.01)
  sd.alpha3 ~ dunif(0, 2)
  tau.alpha3 <- pow(sd.alpha3, -2)
  
  mu.alpha4 ~ dnorm(0, 0.01)
  sd.alpha4 ~ dunif(0, 2)
  tau.alpha4 <- pow(sd.alpha4, -2)
  
  ## Extinction hyperpriors
  mu.eps0 ~ dnorm(0, 0.01)
  sd.eps0 ~ dunif(0, 2)
  tau.eps0 <- pow(sd.eps0, -2)
  
  mu.eps1 ~ dnorm(0, 0.01)
  sd.eps1 ~ dunif(0, 2)
  tau.eps1 <- pow(sd.eps1, -2)
  
  mu.eps2 ~ dnorm(0, 0.01)
  sd.eps2 ~ dunif(0, 2)
  tau.eps2 <- pow(sd.eps2, -2)
  
  mu.eps3 ~ dnorm(0, 0.01)
  sd.eps3 ~ dunif(0, 2)
  tau.eps3 <- pow(sd.eps3, -2)
  
  mu.eps4 ~ dnorm(0, 0.01)
  sd.eps4 ~ dunif(0, 2)
  tau.eps4 <- pow(sd.eps4, -2)
  
  mu.eps5 ~ dnorm(0, 0.01)
  sd.eps5 ~ dunif(0, 2)
  tau.eps5 <- pow(sd.eps5, -2)
  
  mu.eps6 ~ dnorm(0, 0.01)
  sd.eps6 ~ dunif(0, 2)
  tau.eps6 <- pow(sd.eps6, -2)
  
  mu.eps7 ~ dnorm(0, 0.01)
  sd.eps7 ~ dunif(0, 2)
  tau.eps7 <- pow(sd.eps7, -2)
  
  mu.eps8 ~ dnorm(0, 0.01)
  sd.eps8 ~ dunif(0, 2)
  tau.eps8 <- pow(sd.eps8, -2)
  
  mu.eps9 ~ dnorm(0, 0.01)
  sd.eps9 ~ dunif(0, 2)
  tau.eps9 <- pow(sd.eps9, -2)

  mu.eps10 ~ dnorm(0, 0.01)
  sd.eps10 ~ dunif(0, 2)
  tau.eps10 <- pow(sd.eps10, -2)

  mu.eps11 ~ dnorm(0, 0.01)
  sd.eps11 ~ dunif(0, 2)
  tau.eps11 <- pow(sd.eps11, -2)
  
  ## Colonization hyperpriors
  mu.gamma0 ~ dnorm(0, 0.01)
  sd.gamma0 ~ dunif(0, 2)
  tau.gamma0 <- pow(sd.gamma0, -2)
  
  mu.gamma1 ~ dnorm(0, 0.01)
  sd.gamma1 ~ dunif(0, 2)
  tau.gamma1 <- pow(sd.gamma1, -2)
  
  mu.gamma2 ~ dnorm(0, 0.01)
  sd.gamma2 ~ dunif(0, 2)
  tau.gamma2 <- pow(sd.gamma2, -2)
  
  ## Replaced below with mean trait response
  mu.gamma3 ~ dnorm(0, 0.01)
  sd.gamma3 ~ dunif(0, 2)
  tau.gamma3 <- pow(sd.gamma3, -2)
  
  mu.gamma4 ~ dnorm(0, 0.01)
  sd.gamma4 ~ dunif(0, 2)
  tau.gamma4 <- pow(sd.gamma4, -2)
  
  mu.gamma5 ~ dnorm(0, 0.01)
  sd.gamma5 ~ dunif(0, 2)
  tau.gamma5 <- pow(sd.gamma5, -2)
  
  mu.gamma6 ~ dnorm(0, 0.01)
  sd.gamma6 ~ dunif(0, 2)
  tau.gamma6 <- pow(sd.gamma6, -2)
  
  mu.gamma7 ~ dnorm(0, 0.01)
  sd.gamma7 ~ dunif(0, 2)
  tau.gamma7 <- pow(sd.gamma7, -2)
  
  mu.gamma8 ~ dnorm(0, 0.01)
  sd.gamma8 ~ dunif(0, 2)
  tau.gamma8 <- pow(sd.gamma8, -2)
  
  mu.gamma9 ~ dnorm(0, 0.01)
  sd.gamma9 ~ dunif(0, 2)
  tau.gamma9 <- pow(sd.gamma9, -2)

  mu.gamma10 ~ dnorm(0, 0.01)
  sd.gamma10 ~ dunif(0, 2)
  tau.gamma10 <- pow(sd.gamma10, -2)

  mu.gamma11 ~ dnorm(0, 0.01)
  sd.gamma11 ~ dunif(0, 2)
  tau.gamma11 <- pow(sd.gamma11, -2)

  ##----------------------
  ## Ecological State Process
  ##----------------------
  for(k in 1:nspec){
    for(i in 1:nsites){
      
      ## Initial occupancy
      logit(psi1[i,k]) <- beta0[k] + 
        beta1[k] * btemp[i] +
        beta2[k] * btemp2[i] + 
        beta3[k] * bprec[i] +
        beta4[k] * bprec2[i] +
        beta5[k] * cc[i] +
        beta6[k] * cc2[i] 
        #beta7[k] * trendt[i] +
        #beta5[k] * trendtmin[i] +
        #beta6[k] * trendp[i]
      
      z[i,1,k] ~ dbern(psi1[i,k])
      psi[i,1,k] <- psi1[i,k]
      
      ## State process over years
      for(t in 2:nyears){
        ## Colonization
        logit(gamma[i,t-1,k]) <- gamma0[k] +
          gamma1[k] * tanom[i,t-1] + ## JJA tmax anom
          gamma2[k] * panom[i,t-1] + ## JJA sum p anom
          gamma3[k] * hsf_pland[i, t-1] + ## prop hsf (class 3)
          gamma4[k] * tanom[i,t-1] * hsf_pland[i,t-1] + #synergism
          gamma5[k] * panom[i,t-1] * hsf_pland[i,t-1] + #synergism
          gamma6[k] * trendt[i] + # trend tmax
          gamma7[k] * trendtmin[i] + # trend tmin
          gamma8[k] * trendp[i] + # trend precip
          ## Additions for long-term change interactions
          gamma9[k] * trendt[i] * hsf_pland[i, t-1] +
          gamma10[k] * trendtmin[i] * hsf_pland[i, t-1] +
          gamma11[k] * trendp[i] * hsf_pland[i, t-1]
        
        ## Persistence
        logit(eps[i,t-1,k]) <- eps0[k] +
          eps1[k] * tanom[i,t-1] +
          eps2[k] * panom[i,t-1] +
          eps3[k] * hsf_pland[i,t-1] +
          eps4[k] * tanom[i,t-1] * hsf_pland[i,t-1] +
          eps5[k] * panom[i,t-1] * hsf_pland[i,t-1] +
          eps6[k] * trendt[i] +
          eps7[k] * trendtmin[i] +
          eps8[k] * trendp[i] +
          eps9[k] * trendt[i] * hsf_pland[i,t-1] +
          eps10[k] * trendtmin[i] * hsf_pland[i,t-1] +
          eps11[k] * trendp[i] * hsf_pland[i,t-1]
        
        ## Latent state
        z[i,t,k] ~ dbern(z[i,t-1,k] * (1-eps[i,t-1,k]) +
                           (1-z[i,t-1,k]) * gamma[i,t-1,k])
        
        ## Derived psi
        psi[i,t,k] <- psi[i,t-1,k] * (1-eps[i,t-1,k]) +
          (1-psi[i,t-1,k]) * gamma[i,t-1,k]
      }
    }
  }
  
  ## Obs sub-model ragged
  for(v in 1:nObs) {
    for(k in 1:nspec) {
      logit(p_obs[v,k]) <- alpha0[k] +
        alpha1[k] * eff.hrs[v] +
        alpha2[k] * eff.jday[v] + 
        alpha3[k] * eff.jday2[v] +
        #alpha4[k] * ele[site_obs[v]] +
        alpha4[k] * cc[site_obs[v]] + ## per Kristin's Quail paper
        alpha_year[k, year_obs[v]] + 
        cell_det[cell_id[site_obs[v]]]
      
      y[v,k] ~ dbern(z[site_obs[v], year_obs[v], k] * p_obs[v,k] +
                       (1 - z[site_obs[v], year_obs[v], k]) * alphaFP[k])
    }
  }
  
  ##----------------------
  ## Derived parameters
  ##----------------------
  for(k in 1:nspec){
    for(i in 1:nsites){
      for(t in 1:(nyears-1)){
        phi[i,t,k] <- 1 - eps[i,t,k]
      }
    }
  }
  
  ## Mean occupancy per species
  for(k in 1:nspec){
    psi.fs[1,k] <- mean(psi1[1:nsites,k])
    for(t in 2:nyears){
      psi.fs[t,k] <- mean(z[1:nsites,t,k])
    }
  }
  
  ## Mean richness
  for(i in 1:nsites){
    for(t in 1:nyears){
      richness[i,t] <- sum(z[i,t,1:nspec])
    }
  }
  
  ## FP Rate
  # for(k in 1:nspec){
  #   fp_rate[k] <- ilogit(alphaFP[k])
  # }
  
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Data Loading and predictor scaling
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load data
#bdata <- readRDS(here("Data/Occ_Data/DCM_Ragged_Full.RDS"))  
bdata <- readRDS(here("Data/Occ_Data/DCM_Ragged_Full_2021_2025_FireCont_Thresh24.RDS"))  

str(bdata)

## For model validation in the way I have it set up I need some indicator variables
## n_visits = the number of site visits per year x site
## visit_ind = 
## sampled = wide format of whether or not the site has been sampled or not
sampled <- bdata$y_wide
sampled <- apply(sampled, c(1,3), function(x) ifelse(all(is.na(x)), 0, 1))
bdata$sampled <- sampled

## n_visits
nsites <- bdata$nsites
nyears <- bdata$nyears
max_reps <- bdata$nreps
n_visits <- matrix(0, nrow = nsites, ncol = nyears)
visit_ind <- array(NA, dim = c(nsites, nyears, max_reps))
for(i in 1:nsites){
  for(t in 1:nyears){
    v_ids <- which(bdata$site_obs == i & bdata$year_obs == t)
    n_visits[i,t] <- length(v_ids)
    if(n_visits[i,t] > 0){
      visit_ind[i,t,1:n_visits[i,t]] <- v_ids
    }
  }
}
visit_ind[is.na(visit_ind)] <- 1L

#print(n_visits)
#print(visit_ind)

## Replicate indicators
is_valid_visit <- array(0, dim = c(nsites, nyears, max_reps))
for(i in seq_len(nsites)) {
  for(t in seq_len(nyears)) {
    for(r in seq_len(max_reps)) {
      if(r <= n_visits[i,t]) {
        is_valid_visit[i,t,r] <- 1
      }
    }
  }
}


## Add as constants
bdata$n_visits <- n_visits
bdata$visit_ind <- visit_ind
bdata$max_reps <- max_reps
bdata$is_valid_visit <- is_valid_visit

## ***********************************************************
##
## Section Notes: Fixed! Check DCM_Data_Prep.R...misaligned
## from not removing the invalid ARUs at the effort level data
## delete block below after next successful run
##
## ***********************************************************
## Check predictors
## Jday should be > 0
## These should be NA if they match missing values for the actual data
## Something is off here...we should have all zero.hrs = NA for the occ data
# zero.hrs <- which(bdata$eff.hrs == 0)
# zero.day <- which(bdata$eff.jday == 0)
# str(bdata$eff.hrs)
# str(bdata$y)
# bdata$y[zero.hrs, 2]
# bdata$y[zero.day, 2]

## Fix predictor scaling
#bdata$eff.jday2 <- bdata$eff.jday^2
#bdata$cc2 <- bdata$cc^2

mean.ele <- mean(bdata$ele)
sd.ele <- sd(bdata$ele)
bdata$ele <- (bdata$ele - mean.ele)/sd.ele

mean.cc <- mean(bdata$cc)
sd.cc <- sd(bdata$cc)
#bdata$cc <- (bdata$cc - mean.cc)/sd.cc

cc.poly <- poly(bdata$cc, degree = 2)
bdata$cc <- cc.poly[,1]
bdata$cc2 <- cc.poly[,2]

plot(bdata$cc, bdata$cc2)

mean.cc2 <- mean(bdata$cc2)
sd.cc2 <- sd(bdata$cc2)
#bdata$cc2 <- (bdata$cc2 - mean.cc2)/sd.cc2

mean.btemp <- mean(bdata$btmax)
sd.btemp <- sd(bdata$btmax)

btmax.poly <- poly(bdata$btmax, degree = 2)
bdata$btemp <- btmax.poly[,1]
bdata$btemp2 <- btmax.poly[,2]
#bdata$btemp <- (bdata$btmax - mean.btemp)/sd.btemp
plot(bdata$btmax, bdata$btmax2)

mean.bprec <- mean(bdata$bprec)
sd.bprec <- sd(bdata$bprec)

bprec.poly <- poly(bdata$bprec, degree = 2)
bdata$bprec <- bprec.poly[,1]
bdata$bprec2 <- bprec.poly[,2]
#bdata$bprec <- (bdata$bprec - mean.bprec)/sd.bprec
plot(bdata$bprec, bdata$bprec2)

mean.lat <- mean(bdata$Lat)
sd.lat <- sd(bdata$Lat)
bdata$Lat <- (bdata$Lat - mean.lat)/sd.lat

# mean.jday <- mean(bdata$eff.jday)
# sd.jday <- sd(bdata$eff.jday)
# bdata$eff.jday <- (bdata$eff.jday - mean.jday)/sd.jday
# 
# mean.jday2 <- mean(bdata$eff.jday2)
# sd.jday2 <- sd(bdata$eff.jday2)
# bdata$eff.jday2 <- (bdata$eff.jday2 - mean.jday2)/sd.jday2

jday.poly <- poly(bdata$eff.jday, degree = 2)
bdata$eff.jday <- jday.poly[,1]
bdata$eff.jday2 <- jday.poly[,2]

plot(bdata$eff.jday, bdata$eff.jday2)

mean.hrs <- mean(bdata$eff.hrs)
sd.hrs <- sd(bdata$eff.hrs)
bdata$eff.hrs <- (bdata$eff.hrs - mean.hrs)/sd.hrs

sd.panom <- sd(bdata$panom)
bdata$panom <- bdata$panom/sd.panom

sd.tanom <- sd(bdata$tanom)
bdata$tanom <- bdata$tanom/sd.tanom

sd.trendp <- sd(bdata$trendp)
bdata$trendp <- bdata$trendp/sd.trendp

sd.trendt <- sd(bdata$trendt)
bdata$trendt <- bdata$trendt/sd.trendt

sd.trendtmin <- sd(bdata$trendtmin)
bdata$trendtmin <- bdata$trendtmin/sd.trendtmin

## Move from % to proportion
#bdata$hsf_pland15 <- bdata$hsf_pland15/100
#bdata$hsf_pland610 <- bdata$hsf_pland610/100
bdata$hsf_pland <- bdata$hsf_pland/100

## Standardize RMUTINF
mean_fp <- mean(bdata$fire_pattern)
sd_fp <- sd(bdata$fire_pattern)
bdata$fire_pattern <- (bdata$fire_pattern - mean_fp)/sd_fp

#str(bdata)

## Traits
trait <- read.csv(here("Data/Trait/Sierra_HWI_BM_DCM.csv"))
trait <- trait[,-1]
specs <- unlist(dimnames(bdata$y_wide)[4], recursive = T, use.names = F)

## Put the traits together with the combined species
merged_sp <- trait  |> 
  filter(grepl("Sapsucker", Com_Name) | grepl("Cassin's Vireo|Plumbeous Vireo", Com_Name))  |> 
  mutate(Group = ifelse(grepl("Sapsucker", Com_Name), "Sphyrapicus spp.", "Vireo spp.")) |> 
  group_by(Group) |> 
  summarise(
    Family = first(Family),
    Order = first(Order),
    HWI = mean(HWI),
    BM = mean(BM),
    Count = n(),
    .groups = 'drop'
  ) |> 
  rename(Scientific = Group) |>
  mutate(Com_Name = Scientific) |> 
  select(colnames(trait))

## Final trait df
trait <- trait |>
  bind_rows(merged_sp) |> 
  mutate(Com_Name = ifelse(Com_Name == "Western Flycatcher", "Pacific-slope Flycatcher", Com_Name)) |>
  filter(Com_Name %in% specs) |>
  arrange(match(Com_Name, specs)) 

## Extract what we need
bdata$hwi <- as.vector(scale(trait$HWI))
bdata$mass <- as.vector(scale(log(trait$BM)))

str(bdata)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Check collinearity
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# preds <- c("trendt", "trendtmin", "trendp", "tanom", "panom", "hsf_pland", "fire_pattern", "cc", "cc2", "btemp", "btemp2", "bprec", "bprec2")
# 
# cor.pred <- bdata[preds]
# 
# ## Occupancy 
# occ.pred <- data.frame(btemp = cor.pred$btemp,
#                        btemp2 = cor.pred$btemp2,
#                        bprec = cor.pred$bprec,
#                        bprec2 = cor.pred$bprec2,
#                        cc = cor.pred$cc,
#                        cc2 = cor.pred$cc2
# )
# 
# occ.cor <- cor(occ.pred)
# corrplot::corrplot(occ.cor, method = "number")
# 
# dyn.pred <- data.frame(tanom = cor.pred$tanom,
#                        panom = cor.pred$panom,
#                        hsf = cor.pred$hsf_pland,
#                        trendt = cor.pred$trendt,
#                        trendtmin = cor.pred$trendtmin,
#                        trendp = cor.pred$trendp)
# 
# dyn.cor <- cor(dyn.pred)
# corrplot::corrplot(dyn.cor, method = "number")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Inits
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Filling in the gaps
zst <- apply(bdata$y_wide, c(1,3,4), max, na.rm = T)
zst[is.infinite(zst)] <- NA
zst[is.na(zst)] <- 0

z_init <- zst

## Force potential mismatches between zst and y
for (v in seq_len(bdata$nObs)) {
  i <- bdata$site_obs[v]
  t <- bdata$year_obs[v]
  for (k in seq_len(bdata$nspec)) {
    if (!is.na(bdata$y[v, k]) && bdata$y[v, k] == 1) {
      z_init[i, t, k] <- 1
    }
  }
}

str(z_init)
any(is.na(z_init))

inits <- function() {
  
  ncells <- bdata$n_cells
  nspec <- bdata$nspec
  nyears <- bdata$nyears
  
  list(z = z_init,
       
       ## Species-level
       beta0 = rnorm(nspec, 0, 1),
       beta1 = rnorm(nspec, 0, 1),
       beta2 = rnorm(nspec, 0, 1),
       beta3 = rnorm(nspec, 0, 1),
       beta4 = rnorm(nspec, 0, 1),
       beta5 = rnorm(nspec, 0, 1),
       beta6 = rnorm(nspec, 0, 1),
       #beta7 = rnorm(nspec, 0, 1),
       
       alpha0 = rnorm(nspec, 0, 1),
       alpha1 = rnorm(nspec, 0, 1),
       alpha2 = rnorm(nspec, 0, 1),
       alpha3 = rnorm(nspec, 0, 1),
       alpha4 = rnorm(nspec, 0, 1),
       
       gamma0 = rnorm(nspec, 0, 1),
       gamma1 = rnorm(nspec, 0, 1),
       gamma2 = rnorm(nspec, 0, 1),
       gamma3 = rnorm(nspec, 0, 1),
       gamma4 = rnorm(nspec, 0, 1),
       gamma5 = rnorm(nspec, 0, 1),
       gamma6 = rnorm(nspec, 0, 1),
       gamma7 = rnorm(nspec, 0, 1),
       gamma8 = rnorm(nspec, 0, 1),
       gamma9 = rnorm(nspec, 0, 1),
       gamma10 = rnorm(nspec, 0, 1),
       gamma11 = rnorm(nspec, 0, 1),
       
       eps0 = rnorm(nspec, 0, 1),
       eps1 = rnorm(nspec, 0, 1),
       eps2 = rnorm(nspec, 0, 1),
       eps3 = rnorm(nspec, 0, 1),
       eps4 = rnorm(nspec, 0, 1),
       eps5 = rnorm(nspec, 0, 1),
       eps6 = rnorm(nspec, 0, 1),
       eps7 = rnorm(nspec, 0, 1),
       eps8 = rnorm(nspec, 0, 1),
       eps9 = rnorm(nspec, 0, 1),
       eps10 = rnorm(nspec, 0, 1),
       eps11 = rnorm(nspec, 0, 1),
       
       ## Cell raneff
       cell_det = rnorm(ncells, 0, 0.1),
       sd.cell_det = runif(1, 0, 1),
       
       ## Detect year
       alpha_year = array(rnorm(nspec * nyears, 0, 0.1), dim = c(nspec, nyears)),
       sd.alpha_year = runif(1, 0, 1),
       
       ## FP
       #alphaFP <- rnorm(nspec, mean = -2.944, sd = 0.1),
       alphaFP <- rbeta(nspec, shape1 = 5, shape2 = 195),
       
       ## Community level
       mu.beta0 = 0, sd.beta0 = 1,
       mu.beta1 = 0, sd.beta1 = 1,
       mu.beta2 = 0, sd.beta2 = 1,
       mu.beta3 = 0, sd.beta3 = 1,
       mu.beta4 = 0, sd.beta4 = 1,
       mu.beta5 = 0, sd.beta5 = 1,
       mu.beta6 = 0, sd.beta6 = 1,
       #mu.beta7 = 0, sd.beta7 = 1,
       
       mu.alpha0 = 0, sd.alpha0 = 1,
       mu.alpha1 = 0, sd.alpha1 = 1,
       mu.alpha2 = 0, sd.alpha2 = 1,
       mu.alpha3 = 0, sd.alpha3 = 1,
       mu.alpha4 = 0, sd.alpha4 = 1,
       
       mu.eps0 = 0, sd.eps0 = 1,
       mu.eps1 = 0, sd.eps1 = 1,
       mu.eps2 = 0, sd.eps2 = 1,
       mu.eps3 = 0, sd.eps3 = 1,
       mu.eps4 = 0, sd.eps4 = 1,
       mu.eps5 = 0, sd.eps5 = 1,
       mu.eps6 = 0, sd.eps6 = 1,
       mu.eps7 = 0, sd.eps7 = 1,
       mu.eps8 = 0, sd.eps8 = 1,
       #mu.eps9 = 0, sd.eps9 = 1,
       #mu.eps10 = 0, sd.eps10 = 1,
       #mu.eps11 = 0,  sd.eps11 = 1,
       
       mu.gamma0 = 0, sd.gamma0 = 1,
       mu.gamma1 = 0, sd.gamma1 = 1,
       mu.gamma2 = 0, sd.gamma2 = 1,
       mu.gamma3 = 0, sd.gamma3 = 1,
       mu.gamma4 = 0, sd.gamma4 = 1,
       mu.gamma5 = 0, sd.gamma5 = 1,
       mu.gamma6 = 0, sd.gamma6 = 1,
       mu.gamma7 = 0, sd.gamma7 = 1,
       mu.gamma8 = 0, sd.gamma8 = 1)#,
       #mu.gamma9 = 0, sd.gamma9 = 1)#,
       #mu.gamma10 = 0, sd.gamma10 = 1,
       #mu.gamma11 = 0, sd.gamma11 = 1)
}

## Inits for trait model
# inits <- function() {
#   
#   ncells <- bdata$n_cells
#   nspec <- bdata$nspec
#   nyears <- bdata$nyears
#   
#   list(z = z_init,
#        
#        ## Species-level
#        beta0 = rnorm(nspec, 0, 1),
#        beta1 = rnorm(nspec, 0, 1),
#        beta2 = rnorm(nspec, 0, 1),
#        beta3 = rnorm(nspec, 0, 1),
#        beta4 = rnorm(nspec, 0, 1),
#        beta5 = rnorm(nspec, 0, 1),
#        beta6 = rnorm(nspec, 0, 1),
#        
#        alpha0 = rnorm(nspec, 0, 1),
#        alpha1 = rnorm(nspec, 0, 1),
#        alpha2 = rnorm(nspec, 0, 1),
#        alpha3 = rnorm(nspec, 0, 1),
#        alpha4 = rnorm(nspec, 0, 1),
#        
#        gamma0 = rnorm(nspec, 0, 1),
#        gamma1 = rnorm(nspec, 0, 1),
#        gamma2 = rnorm(nspec, 0, 1),
#        gamma3 = rnorm(nspec, 0, 1),
#        gamma4 = rnorm(nspec, 0, 1),
#        gamma5 = rnorm(nspec, 0, 1),
#        gamma6 = rnorm(nspec, 0, 1),
#        gamma7 = rnorm(nspec, 0, 1),
#        gamma8 = rnorm(nspec, 0, 1),
#        gamma9 = rnorm(nspec, 0, 1),
#        gamma10 = rnorm(nspec, 0, 1),
#        gamma11 = rnorm(nspec, 0, 1),
#        
#        eps0 = rnorm(nspec, 0, 1),
#        eps1 = rnorm(nspec, 0, 1),
#        eps2 = rnorm(nspec, 0, 1),
#        eps3 = rnorm(nspec, 0, 1),
#        eps4 = rnorm(nspec, 0, 1),
#        eps5 = rnorm(nspec, 0, 1),
#        eps6 = rnorm(nspec, 0, 1),
#        eps7 = rnorm(nspec, 0, 1),
#        eps8 = rnorm(nspec, 0, 1),
#        eps9 = rnorm(nspec, 0, 1),
#        eps10 = rnorm(nspec, 0, 1),
#        eps11 = rnorm(nspec, 0, 1),
#        
#        ## Cell raneff
#        cell_det = rnorm(ncells, 0, 0.1),
#        sd.cell_det = runif(1, 0, 1),
#        
#        ## Detect year
#        alpha_year = rnorm(nspec * nyears, 0, 0.1),
#        sd.alpha_year = runif(1, 0, 1),
#        
#        ## FP
#        #alphaFP <- rnorm(nspec, mean = -2.944, sd = 0.1),
#        alphaFP <- rbeta(nspec, shape1 = 5, shape2 = 195),
#        
#        ## Community level
#        mu.beta0 = 0,  sd.beta0 = 1,
#        mu.beta1 = 0,  sd.beta1 = 1,
#        mu.beta2 = 0,  sd.beta2 = 1,
#        mu.beta3 = 0,  sd.beta3 = 1,
#        mu.beta4 = 0,  sd.beta4 = 1,
#        mu.beta5 = 0,  sd.beta5 = 1,
#        mu.beta6 = 0,  sd.beta6 = 1,
#        
#        mu.alpha0 = 0, sd.alpha0 = 1,
#        mu.alpha1 = 0, sd.alpha1 = 1,
#        mu.alpha2 = 0, sd.alpha2 = 1,
#        mu.alpha3 = 0, sd.alpha3 = 1,
#        mu.alpha4 = 0, sd.alpha4 = 1,
#        
#        mu.eps0 = rep(0, nspec),  sd.eps0 = 1,
#        mu.eps1 = rep(0, nspec),  sd.eps1 = 1,
#        mu.eps2 = rep(0, nspec),  sd.eps2 = 1,
#        mu.eps3 = rep(0, nspec),  sd.eps3 = 1,
#        mu.eps4 = rep(0, nspec),  sd.eps4 = 1,
#        mu.eps5 = rep(0, nspec),  sd.eps5 = 1,
#        mu.eps6 = rep(0, nspec),  sd.eps6 = 1,
#        mu.eps7 = rep(0, nspec),  sd.eps7 = 1,
#        mu.eps8 = rep(0, nspec),  sd.eps8 = 1,
#        mu.eps9 = rep(0, nspec),  sd.eps9 = 1,
#        mu.eps10 = rep(0, nspec),  sd.eps10 = 1,
#        mu.eps11 = rep(0, nspec),  sd.eps11 = 1,
#        
#        mu.gamma0 = rep(0, nspec), sd.gamma0 = 1,
#        mu.gamma1 = rep(0, nspec), sd.gamma1 = 1,
#        mu.gamma2 = rep(0, nspec), sd.gamma2 = 1,
#        mu.gamma3 = rep(0, nspec), sd.gamma3 = 1,
#        mu.gamma4 = rep(0, nspec), sd.gamma4 = 1,
#        mu.gamma5 = rep(0, nspec), sd.gamma5 = 1,
#        mu.gamma6 = rep(0, nspec), sd.gamma6 = 1,
#        mu.gamma7 = rep(0, nspec), sd.gamma7 = 1,
#        mu.gamma8 = rep(0, nspec), sd.gamma8 = 1,
#        mu.gamma9 = rep(0, nspec), sd.gamma9 = 1,
#        mu.gamma10 = rep(0, nspec), sd.gamma10 = 1,
#        mu.gamma11 = rep(0, nspec), sd.gamma11 = 1,
#        
#        ## Trait model inits
#        delta0.gamma0  = rnorm(1, 0, 1),
#        delta1.gamma0  = rnorm(1, 0, 1),
#        delta2.gamma0  = rnorm(1, 0, 1),
#        delta0.gamma1  = rnorm(1, 0, 1),
#        delta1.gamma1  = rnorm(1, 0, 1),
#        delta2.gamma1  = rnorm(1, 0, 1),
#        delta0.gamma2  = rnorm(1, 0, 1),
#        delta1.gamma2  = rnorm(1, 0, 1),
#        delta2.gamma2  = rnorm(1, 0, 1),
#        delta0.gamma3  = rnorm(1, 0, 1),
#        delta1.gamma3  = rnorm(1, 0, 1),
#        delta2.gamma3  = rnorm(1, 0, 1),
#        delta0.gamma4  = rnorm(1, 0, 1),
#        delta1.gamma4  = rnorm(1, 0, 1),
#        delta2.gamma4  = rnorm(1, 0, 1),
#        delta0.gamma5  = rnorm(1, 0, 1),
#        delta1.gamma5  = rnorm(1, 0, 1),
#        delta2.gamma5  = rnorm(1, 0, 1),
#        delta0.gamma6  = rnorm(1, 0, 1),
#        delta1.gamma6  = rnorm(1, 0, 1),
#        delta2.gamma6  = rnorm(1, 0, 1),
#        delta0.gamma7  = rnorm(1, 0, 1),
#        delta1.gamma7  = rnorm(1, 0, 1),
#        delta2.gamma7  = rnorm(1, 0, 1),
#        delta0.gamma8  = rnorm(1, 0, 1),
#        delta1.gamma8  = rnorm(1, 0, 1),
#        delta2.gamma8  = rnorm(1, 0, 1),
#        delta0.gamma9  = rnorm(1, 0, 1),
#        delta1.gamma9  = rnorm(1, 0, 1),
#        delta2.gamma9  = rnorm(1, 0, 1),
#        delta0.gamma10 = rnorm(1, 0, 1),
#        delta1.gamma10 = rnorm(1, 0, 1),
#        delta2.gamma10 = rnorm(1, 0, 1),
#        delta0.gamma11 = rnorm(1, 0, 1),
#        delta1.gamma11 = rnorm(1, 0, 1),
#        delta2.gamma11 = rnorm(1, 0, 1),
# 
#        delta0.eps0  = rnorm(1, 0, 1),
#        delta1.eps0  = rnorm(1, 0, 1),
#        delta2.eps0  = rnorm(1, 0, 1),
#        delta0.eps1  = rnorm(1, 0, 1),
#        delta1.eps1  = rnorm(1, 0, 1),
#        delta2.eps1  = rnorm(1, 0, 1),
#        delta0.eps2  = rnorm(1, 0, 1),
#        delta1.eps2  = rnorm(1, 0, 1),
#        delta2.eps2  = rnorm(1, 0, 1),
#        delta0.eps3  = rnorm(1, 0, 1),
#        delta1.eps3  = rnorm(1, 0, 1),
#        delta2.eps3  = rnorm(1, 0, 1),
#        delta0.eps4  = rnorm(1, 0, 1),
#        delta1.eps4  = rnorm(1, 0, 1),
#        delta2.eps4  = rnorm(1, 0, 1),
#        delta0.eps5  = rnorm(1, 0, 1),
#        delta1.eps5  = rnorm(1, 0, 1),
#        delta2.eps5  = rnorm(1, 0, 1),
#        delta0.eps6  = rnorm(1, 0, 1),
#        delta1.eps6  = rnorm(1, 0, 1),
#        delta2.eps6  = rnorm(1, 0, 1),
#        delta0.eps7  = rnorm(1, 0, 1),
#        delta1.eps7  = rnorm(1, 0, 1),
#        delta2.eps7  = rnorm(1, 0, 1),
#        delta0.eps8  = rnorm(1, 0, 1),
#        delta1.eps8  = rnorm(1, 0, 1),
#        delta2.eps8  = rnorm(1, 0, 1),
#        delta0.eps9  = rnorm(1, 0, 1),
#        delta1.eps9  = rnorm(1, 0, 1),
#        delta2.eps9  = rnorm(1, 0, 1),
#        delta0.eps10 = rnorm(1, 0, 1),
#        delta1.eps10 = rnorm(1, 0, 1),
#        delta2.eps10 = rnorm(1, 0, 1),
#        delta0.eps11 = rnorm(1, 0, 1),
#        delta1.eps11 = rnorm(1, 0, 1),
#        delta2.eps11 = rnorm(1, 0, 1))
# }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Monitors
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Parameters monitored
nspec <- bdata$nspec

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: MCMC Settings
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nburnin = 10000
niter = 40000
nthin = 20
nchain = 3

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Model Compilation, MCMC Config, and Building
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DCMmodel <- nimbleModel(
  code = DCMContFire,
  data = list(y = bdata$y),   ## match name in nimbleCode
  constants = bdata[names(bdata) != "y"],
  inits = inits()
)

## Check the node names
node_names <- DCMmodel$getNodeNames()
var_names <- DCMmodel$getVarNames()

## Grep the names we want from available nodes
nodes_by_pattern <- function(pattern) {
  grep(pattern, node_names, value = TRUE)
}

## Adding different monitors for two flavors of runs
## Effects monitors
mon_beta  <- nodes_by_pattern("^beta[0-9]+\\[")
mon_alpha <- nodes_by_pattern("^alpha[0-9]+\\[")
mon_gamma <- nodes_by_pattern("^gamma[0-9]+\\[")
mon_eps   <- nodes_by_pattern("^eps[0-9]+\\[")
mon_mu <- nodes_by_pattern("^mu.")
mon_cell <- nodes_by_pattern("cell_det\\[")
mon_yeardet <- nodes_by_pattern("alpha_year\\[")
#mon_delta <- nodes_by_pattern("delta")

## Single parameters
singles <- c("sd.cell_det", 
             "richness", 
             "alphaFP",
             "sd.alpha_year"
             #"fp_rate"
)

## Put all together for params
params_eff <- c(mon_beta, 
                mon_alpha, 
                mon_gamma, 
                mon_eps, 
                mon_mu, 
                mon_cell,
                mon_yeardet,
                singles)

## Monitors for PPCs
# mon_z <- nodes_by_pattern("z\\[")
# mon_col <- nodes_by_pattern("gamma\\[")
# mon_ext <- nodes_by_pattern("eps\\[")
# 
# params_ppc <- c(mon_z, 
#                 mon_col, 
#                 mon_ext,
#                 "alphaFP",
#                 "p_obs")
# 
# 
# length(params) * (niter / nthin) * nchain
## Compile
Cmodel <- compileNimble(DCMmodel)

## Configure model
conf <- configureMCMC(DCMmodel, monitors = params_eff)

#conf$enableWAIC(TRUE)
print(conf$getMonitors())

## Build
Rmcmc <- buildMCMC(conf)

## Compile 2x
Cmcmc <- compileNimble(Rmcmc, project = DCMmodel)

# # ## Run for a dry 
# Cmcmc$run(niter = 1000, nburnin = 500)
# 
# samples_dry <- as.matrix(Cmcmc$mvSamples)
# head(samples_dry)

samples <- runMCMC(Cmcmc,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchain,
                   samplesAsCodaMCMC = FALSE,
                   summary = FALSE)

## Write to file
saveRDS(samples, file = "D:/DCM_Samples/DCMmodel_mcmc_output_FP_21_25_Thresh24_NoTraits_V3.rds")