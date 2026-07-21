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
library(stringr)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: Model Code
##
## -------------------------------------------------------------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: DCM Continuous Fire No traits
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
  ## Note: Strong prior on low FP rate
  ## to try to avoid multi-modality and sign switching
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
        beta1[k] * btmax_annual_sc[i] +
        beta2[k] * btmax_annual_sc2[i] + 
        beta3[k] * bprec_annual_sc[i] +
        beta4[k] * bprec_annual_sc2[i] +
        beta5[k] * cc_sc[i] +
        beta6[k] * cc_sc2[i]
      
      z[i,1,k] ~ dbern(psi1[i,k])
      psi[i,1,k] <- psi1[i,k]
      
      ## State process over years
      for(t in 2:nyears){
        ## Colonization
        logit(gamma[i,t-1,k]) <- gamma0[k] +
          gamma1[k] * tmax_anom_mam_sc[i, t-1] +
          gamma2[k] * p_anom_mam_sc[i, t-1] +
          gamma3[k] * hsf_pland1_10[i,t-1] +
          gamma4[k] * tmax_anom_mam_sc[i, t-1] * hsf_pland1_10[i,t-1] + #synergism
          gamma5[k] * p_anom_mam_sc[i, t-1] * hsf_pland1_10[i,t-1] + #synergism
          gamma6[k] * tmax_trend_mam_sc[i] + 
          gamma7[k] * tmin_trend_mam_sc[i] + 
          gamma8[k] * prec_trend_mam_sc[i] + 
          ## Additions for long-term change interactions
          gamma9[k] * tmax_trend_mam_sc[i] * hsf_pland1_10[i, t-1] +
          gamma10[k] * tmin_trend_mam_sc[i] * hsf_pland1_10[i, t-1] +
          gamma11[k] * prec_trend_mam_sc[i] * hsf_pland1_10[i, t-1]
        
        ## Persistence
        logit(eps[i,t-1,k]) <-  eps0[k] +
          eps1[k] * tmax_anom_mam_sc[i, t-1] +
          eps2[k] * p_anom_mam_sc[i, t-1] +
          eps3[k] * hsf_pland1_10[i,t-1] +
          eps4[k] * tmax_anom_mam_sc[i, t-1] * hsf_pland1_10[i,t-1] +
          eps5[k] * p_anom_mam_sc[i, t-1] * hsf_pland1_10[i,t-1] +
          eps6[k] * tmax_trend_mam_sc[i] + 
          eps7[k] * tmin_trend_mam_sc[i] + 
          eps8[k] * prec_trend_mam_sc[i] + 
          ## Additions for long-term change interactions
          eps9[k] * tmax_trend_mam_sc[i] * hsf_pland1_10[i, t-1] +
          eps10[k] * tmin_trend_mam_sc[i] * hsf_pland1_10[i, t-1] +
          eps11[k] * prec_trend_mam_sc[i] * hsf_pland1_10[i, t-1]
        
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
        alpha1[k] * eff.hrs_sc[v] +
        alpha2[k] * eff.jday_sc[v] + 
        alpha3[k] * eff.jday_sc2[v] +
        alpha4[k] * cc_sc[site_obs[v]] +
        alpha_year[k, year_obs[v]] + 
        cell_det[cell_id[site_obs[v]]]
      
      ## Observed data model w/ FPs built in
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
  
})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: DCM Climate only with rel temp
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DCMClimCh <- nimbleCode({
  
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
  ## Note: Strong prior on low FP rate
  ## to try to avoid multi-modality and sign switching
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
        beta1[k] * btmax_annual_sc[i] +
        beta2[k] * btmax_annual_sc2[i] + 
        beta3[k] * bprec_annual_sc[i] +
        beta4[k] * bprec_annual_sc2[i] +
        beta5[k] * cc_sc[i] +
        beta6[k] * cc_sc2[i]
      
      z[i,1,k] ~ dbern(psi1[i,k])
      psi[i,1,k] <- psi1[i,k]
      
      ## State process over years
      for(t in 2:nyears){
        ## Colonization
        logit(gamma[i,t-1,k]) <- gamma0[k] +
          gamma1[k] * tmax_anom_mam_sc[i, t-1] +
          gamma2[k] * p_anom_mam_sc[i, t-1] +
          gamma3[k] * rt[i] +
          gamma4[k] * tmax_anom_mam_sc[i, t-1] * rt[i] + #synergism
          gamma5[k] * p_anom_mam_sc[i, t-1] * rp[i] + #synergism
          gamma6[k] * tmax_trend_mam_sc[i] + 
          gamma7[k] * tmin_trend_mam_sc[i] + 
          gamma8[k] * prec_trend_mam_sc[i] + 
          ## Additions for long-term change interactions
          gamma9[k] * tmax_trend_mam_sc[i] * rt[i] +
          gamma10[k] * tmin_trend_mam_sc[i] * rt[i] +
          gamma11[k] * prec_trend_mam_sc[i] * rp[i]
        
        ## Persistence
        logit(eps[i,t-1,k]) <-  eps0[k] +
          eps1[k] * tmax_anom_mam_sc[i, t-1] +
          eps2[k] * p_anom_mam_sc[i, t-1] +
          eps3[k] * rt[i] +
          eps4[k] * tmax_anom_mam_sc[i, t-1] * rt[i] +
          eps5[k] * p_anom_mam_sc[i, t-1] * rp[i] +
          eps6[k] * tmax_trend_mam_sc[i] + 
          eps7[k] * tmin_trend_mam_sc[i] + 
          eps8[k] * prec_trend_mam_sc[i] + 
          ## Additions for long-term change interactions
          eps9[k] * tmax_trend_mam_sc[i] * rt[i] +
          eps10[k] * tmin_trend_mam_sc[i] * rt[i] +
          eps11[k] * prec_trend_mam_sc[i] * rp[i]
        
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
        alpha1[k] * eff.hrs_sc[v] +
        alpha2[k] * eff.jday_sc[v] + 
        alpha3[k] * eff.jday_sc2[v] +
        alpha4[k] * cc_sc[site_obs[v]] +
        alpha_year[k, year_obs[v]] + 
        cell_det[cell_id[site_obs[v]]]
      
      ## Observed data model w/ FPs built in
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
#bdata <- readRDS(here("Data/Occ_Data/DCM_Ragged_Full_2021_2024_RevisedCovs_Tester4yrOnly.RDS"))  
bdata <- readRDS(here("Data/Occ_Data/DCM_Ragged_Complete_Survey_Only.RDS"))  

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

## ***********************************************************
##
## Section Notes: Fixed! Check DCM_Data_Prep.R...misaligned
## from not removing the invalid ARUs at the effort level data
## delete block below after next successful run
##
## ***********************************************************
## Predictor scaling

## Canopy Cover
bdata$cc_sc <- scale(bdata$cc.nldc)[,1]
bdata$cc_sc2 <- bdata$cc_sc^2
plot(bdata$cc_sc, bdata$cc_sc2)

## Baseline Tmax
bdata$btmax_annual_sc <- scale(bdata$btmax_annual)[,1] 
bdata$btmax_annual_sc2 <- bdata$btmax_annual_sc^2
plot(bdata$btmax_annual_sc, bdata$btmax_annual_sc2)

## Baseline Prec
bdata$bprec_annual_sc <- scale(bdata$bprec_annual)[,1]
bdata$bprec_annual_sc2 <- bdata$bprec_annual_sc^2
plot(bdata$bprec_annual_sc, bdata$bprec_annual_sc2)

## Scale the trends, but don't center
bdata$tmin_trend_jja_sc <- scale(bdata$tmin_trend_jja, center = FALSE)[,1]
bdata$tmax_trend_jja_sc <- scale(bdata$tmax_trend_jja, center = FALSE)[,1]
bdata$prec_trend_jja_sc <- scale(bdata$prec_trend_jja, center = FALSE)[,1]

# hist(bdata$tmin_trend_jja_sc)
# hist(bdata$tmax_trend_jja_sc)
# hist(bdata$prec_trend_jja_sc)

bdata$tmin_trend_mam_sc <- scale(bdata$tmin_trend_mam, center = FALSE)[,1]
bdata$tmax_trend_mam_sc <- scale(bdata$tmax_trend_mam, center = FALSE)[,1]
bdata$prec_trend_mam_sc <- scale(bdata$prec_trend_mam, center = FALSE)[,1]

# hist(bdata$tmin_trend_mam_sc)
# hist(bdata$tmax_trend_mam_sc)
# hist(bdata$prec_trend_mam_sc)

## Detection covariates
## Jday
bdata$eff.jday_sc <- scale(bdata$eff.jday)[,1]
bdata$eff.jday_sc2 <- bdata$eff.jday_sc^2
plot(bdata$eff.jday_sc, bdata$eff.jday_sc2)

## Hours
bdata$eff.hrs_sc <- scale(bdata$eff.hrs)[,1]

## Dynamic covariates
## Scale panomalies and trends by sd but leave intercept intact
bdata$p_anom_jja_sc <- scale(bdata$p_anom_jja, center = FALSE)[,2:nyears]
bdata$tmax_anom_jja_sc <- scale(bdata$tmax_anom_jja, center = FALSE)[,2:nyears]

bdata$p_anom_mam_sc <- scale(bdata$p_anom_mam, center = FALSE)[,2:nyears]
bdata$tmax_anom_mam_sc <- scale(bdata$tmax_anom_mam, center = FALSE)[,2:nyears]

## Fire data
## Move from % to proportion
#bdata$hsf_pland15 <- bdata$hsf_pland15/100
#bdata$hsf_pland610 <- bdata$hsf_pland610/100
bdata$hsf_pland1_10 <- as.matrix(bdata$hsf_pland1_10[,2:nyears])
mean_hsf <- mean(bdata$hsf_pland1_10, na.rm = TRUE)
sd_hsf <- sd(bdata$hsf_pland1_10, na.rm = TRUE)
bdata$hsf_pland1_10 <- (bdata$hsf_pland1_10 - mean_hsf) / sd_hsf
hist(bdata$hsf_pland1_10)

## TSF
bdata$tsf_hsf <- as.matrix(bdata$tsf_hsf)
bdata$tsf_hsf_sc <- scale(bdata$tsf_hsf)
hist(bdata$tsf_hsf_sc)

str(bdata)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Check collinearity
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor.pred <- bdata[str_detect(names(bdata), pattern = "_sc|trend")]
str(cor.pred)

## Occupancy
occ.pred <- data.frame(btemp = cor.pred$btmax_annual_sc,
                       btemp2 = cor.pred$btmax_annual_sc2,
                       bprec = cor.pred$bprec_annual_sc,
                       bprec2 = cor.pred$bprec_annual_sc2,
                       cc = cor.pred$cc_sc,
                       cc2 = cor.pred$cc_sc2
)

occ.cor <- cor(occ.pred)
corrplot::corrplot(occ.cor, method = "number")

dyn.pred <- data.frame(tanom = cor.pred$tmax_anom_mam_sc,
                       panom = cor.pred$p_anom_mam_sc,
                       tsf = cor.pred$tsf_hsf_sc,
                       trendt = cor.pred$tmax_trend_mam,
                       trendtmin = cor.pred$tmin_trend_mam,
                       trendp = cor.pred$prec_trend_mam)

dyn.cor <- cor(dyn.pred)
corrplot::corrplot(dyn.cor, method = "number")

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Intialization
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function() {
  
  ncells <- bdata$n_cells
  nspec <- bdata$nspec
  nyears <- bdata$nyears
  
  list(z = z_init,
       
       ## Species-level
       ## Occupancy coeffs
       beta0 = rnorm(nspec, 0, 1),
       beta1 = rnorm(nspec, 0, 1),
       beta2 = rnorm(nspec, 0, 1),
       beta3 = rnorm(nspec, 0, 1),
       beta4 = rnorm(nspec, 0, 1),
       beta5 = rnorm(nspec, 0, 1),
       beta6 = rnorm(nspec, 0, 1),
       
       ## Detection coeffs
       alpha0 = rnorm(nspec, 0, 1),
       alpha1 = rnorm(nspec, 0, 1),
       alpha2 = rnorm(nspec, 0, 1),
       alpha3 = rnorm(nspec, 0, 1),
       alpha4 = rnorm(nspec, 0, 1),
       
       ## Colonization coeffs
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
       
       ## Extinction coeffs
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
       mu.eps9 = 0, sd.eps9 = 1,
       mu.eps10 = 0, sd.eps10 = 1,
       mu.eps11 = 0,  sd.eps11 = 1,
       
       mu.gamma0 = 0, sd.gamma0 = 1,
       mu.gamma1 = 0, sd.gamma1 = 1,
       mu.gamma2 = 0, sd.gamma2 = 1,
       mu.gamma3 = 0, sd.gamma3 = 1,
       mu.gamma4 = 0, sd.gamma4 = 1,
       mu.gamma5 = 0, sd.gamma5 = 1,
       mu.gamma6 = 0, sd.gamma6 = 1,
       mu.gamma7 = 0, sd.gamma7 = 1,
       mu.gamma8 = 0, sd.gamma8 = 1,
       mu.gamma9 = 0, sd.gamma9 = 1,
       mu.gamma10 = 0, sd.gamma10 = 1,
       mu.gamma11 = 0, sd.gamma11 = 1)
}

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
niter = 30000
nthin = 20
nchain = 3

## Total posterior samples
((niter - nburnin) / nthin) * nchain

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

## Compile (allow Nimble to assign default samplers)
Cmodel <- compileNimble(DCMmodel)

## Configure model
conf <- configureMCMC(DCMmodel, monitors = params_eff)

#conf$enableWAIC(TRUE)
print(conf$getMonitors())

## Build
Rmcmc <- buildMCMC(conf)

## Compile 2x
Cmcmc <- compileNimble(Rmcmc, project = DCMmodel)

# ## Run for a test 
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
saveRDS(samples, file = "D:/DCM_Samples/DCMmodel_mcmc_output_FP_21_25_Spring_Climate_CompSurv.rds")

range(bdata$tmax_trend_mam, na.rm = TRUE)
range(bdata$hsf_pland1_10, na.rm = TRUE)
hist(as.matrix(bdata$tmax_trend_mam * bdata$hsf_pland1_10[,-1]))  ## interaction values
