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

# dynOcc <- nimbleCode({
#   
#   ## Priors
#   ## Occupancy coefficients
#   beta0 ~ dnorm(0, 0.01)
#   beta1 ~ dnorm(0, 0.01)
#   beta2 ~ dnorm(0, 0.01)
#   
#   ## Detection coefficients
#   alpha0 ~ dnorm(0, 0.01)
#   alpha1 ~ dnorm(0, 0.01)
#   alpha2 ~ dnorm(0, 0.01)
#   
#   ## Persistence coefficients
#   eps0 ~ dnorm(0, 0.01)
#   eps1 ~ dnorm(0, 0.01)
#   eps2 ~ dnorm(0, 0.01)
#   eps3 ~ dnorm(0, 0.01)
#   
#   ## Colonization coefficients
#   gamma0 ~ dnorm(0, 0.01)
#   gamma1 ~ dnorm(0, 0.01)
#   gamma2 ~ dnorm(0, 0.01)
#   gamma3 ~ dnorm(0, 0.01)
#   
#   ## Likelihood
#   ## Ecological submodel
#   for(i in 1:nsites){
#     
#     ## Predictors for intial occ
#     logit(psi1[i]) <- beta0 + beta1 * bclim[i, 1] + beta2 * cc[i,1]
#     
#     ## Initial occ
#     z[i, 1] ~ dbern(psi1[i])
#     
#     ## State process model (colex)
#     for(t in 2:nyears){
#       
#       ## Colonization
#       logit(gamma[i, t-1]) <- gamma0 + gamma1 * anom[i, t-1] + gamma2 * fire[i, t-1] + gamma3 * anom[i, t-1] * fire[i, t-1]
#       
#       ## Persistence
#       logit(eps[i, t-1]) <- eps0 + eps1 * anom[i, t-1] + eps2 * fire[i, t-1] + eps3 * anom[i, t-1] * fire[i, t-1]
#       
#       ## Latent state
#       z[i,t] ~ dbern(z[i, t-1] * (1 - eps[i, t-1]) + (1-z[i,t-1]) * gamma[i, t-1])
#     }
#   }
#   
#   ## Observation submodel
#   for(i in 1:nsites){
#     for(j in 1:nreps){
#       for(t in 1:nyears){
#         logit(p[i,j,t]) <- alpha0 + alpha1 * eff.hrs[i,j,t] + alpha2 * eff.jday[i,j,t]
#         y[i,j,t] ~ dbern(z[i,t] * p[i,j,t])
#       }
#     }
#   }
#   
#   ## Derived parameters
#   for(i in 1:nsites){
#     for(t in 1:(nyears-1)){
#       phi[i, t] <- 1 - eps[i, t]
#     }
#   }
#   
#   ## Mean occupancy
#   psi.fs[1] <- mean(psi1[1:nsites])
#   for(t in 2:nyears){
#     psi.fs[t] <- mean(z[1:nsites, t])
#   }
#   
# })
# 
# dynOccSite <- nimbleCode({
#   
#   ## Priors
#   ## Occupancy coefficients
#   beta0 ~ dnorm(0, 0.01)
#   beta1 ~ dnorm(0, 0.01)
#   beta2 ~ dnorm(0, 0.01)
#   
#   ## Raneff for Site
#   sigma.site ~ dunif(0, 10)
#   tau.site <- 1/(sigma.site * sigma.site)
#   
#   for(i in 1:nsite){
#     eps.site[i] ~ dnorm(0, tau.site)
#   }
#   
#   ## Detection coefficients
#   alpha0 ~ dnorm(0, 0.01)
#   alpha1 ~ dnorm(0, 0.01)
#   alpha2 ~ dnorm(0, 0.01)
#   
#   ## Persistence coefficients
#   eps0 ~ dnorm(0, 0.01)
#   eps1 ~ dnorm(0, 0.01)
#   eps2 ~ dnorm(0, 0.01)
#   eps3 ~ dnorm(0, 0.01)
#   
#   ## Colonization coefficients
#   gamma0 ~ dnorm(0, 0.01)
#   gamma1 ~ dnorm(0, 0.01)
#   gamma2 ~ dnorm(0, 0.01)
#   gamma3 ~ dnorm(0, 0.01)
#   
#   ## Likelihood
#   ## Ecological submodel
#   for(i in 1:nsites){
#     
#     ## Predictors for intial occ
#     logit(psi1[i]) <- beta0 + beta1 * bclim[i, 1] + beta2 * cc[i,1] + eps.site[i]
#     
#     ## Initial occ
#     z[i, 1] ~ dbern(psi1[i])
#     
#     ## State process model (colex)
#     for(t in 2:nyears){
#       
#       ## Colonization
#       logit(gamma[i, t-1]) <- gamma0 + gamma1 * anom[i, t-1] + gamma2 * fire[i, t-1] + gamma3 * anom[i, t-1] * fire[i, t-1] + eps.site[i]
#       
#       ## Persistence
#       logit(eps[i, t-1]) <- eps0 + eps1 * anom[i, t-1] + eps2 * fire[i, t-1] + eps3 * anom[i, t-1] * fire[i, t-1] + eps.site[i]
#       
#       ## Latent state
#       z[i,t] ~ dbern(z[i, t-1] * (1 - eps[i, t-1]) + (1-z[i,t-1]) * gamma[i, t-1])
#     }
#   }
#   
#   ## Observation submodel
#   for(i in 1:nsites){
#     for(j in 1:nreps){
#       for(t in 1:nyears){
#         logit(p[i,j,t]) <- alpha0 + alpha1 * eff.hrs[i,j,t] + alpha2 * eff.jday[i,j,t]
#         y[i,j,t] ~ dbern(z[i,t] * p[i,j,t])
#       }
#     }
#   }
#   
#   ## Derived parameters
#   for(i in 1:nsites){
#     for(t in 1:(nyears-1)){
#       phi[i, t] <- 1 - eps[i, t]
#     }
#   }
#   
#   ## Mean occupancy
#   psi.fs[1] <- mean(psi1[1:nsites])
#   for(t in 2:nyears){
#     psi.fs[t] <- mean(z[1:nsites, t])
#   }
#   
# })
# 
# dynOccVert <- nimbleCode({
#   
#   ## Priors
#   ## Occupancy coefficients
#   beta0 ~ dnorm(0, 0.01)
#   beta1 ~ dnorm(0, 0.01)
#   beta2 ~ dnorm(0, 0.01)
#   
#   ## Detection coefficients
#   alpha0 ~ dnorm(0, 0.01)
#   alpha1 ~ dnorm(0, 0.01)
#   alpha2 ~ dnorm(0, 0.01)
#   
#   ## Persistence coefficients
#   eps0 ~ dnorm(0, 0.01)
#   eps1 ~ dnorm(0, 0.01)
#   eps2 ~ dnorm(0, 0.01)
#   eps3 ~ dnorm(0, 0.01)
#   
#   ## Colonization coefficients
#   gamma0 ~ dnorm(0, 0.01)
#   gamma1 ~ dnorm(0, 0.01)
#   gamma2 ~ dnorm(0, 0.01)
#   gamma3 ~ dnorm(0, 0.01)
#   
#   ## Likelihood
#   ## Ecological submodel
#   for(i in 1:nsites){
#     
#     ## Predictors for initial occ
#     logit(psi1[i]) <- beta0 + beta1 * bclim[i] + beta2 * cc[i]
#     
#     ## Initial occ
#     z[i, 1] ~ dbern(psi1[i])
#     
#     ## State process model (colex)
#     for(t in 2:nyears){
#       
#       ## Colonization
#       logit(gamma[i, t-1]) <- gamma0 + gamma1 * temp_anomaly[i, t-1] + 
#         gamma2 * fire_severity[i, t-1] + 
#         gamma3 * temp_anomaly[i, t-1] * fire_severity[i, t-1]
#       
#       ## Persistence
#       logit(eps[i, t-1]) <- eps0 + eps1 * temp_anomaly[i, t-1] + 
#         eps2 * fire_severity[i, t-1] + 
#         eps3 * temp_anomaly[i, t-1] * fire_severity[i, t-1]
#       
#       ## Latent state
#       z[i,t] ~ dbern(z[i, t-1] * (1 - eps[i, t-1]) + (1-z[i,t-1]) * gamma[i, t-1])
#     }
#   }
#   
#   ## Observation submodel
#   for(i in 1:nobs){
#     logit(p[i]) <- alpha0 + alpha1 * hours[i] + alpha2 * date[i]
#     y[i] ~ dbern(z[site[i], year[i]] * p[i])
#   }
#   
#   ## Derived parameters
#   for(i in 1:nsites){
#     for(t in 1:(nyears-1)){
#       phi[i, t] <- 1 - eps[i, t]
#     }
#   }
#   
#   ## Mean occupancy
#   psi.fs[1] <- mean(psi1[1:nsites])
#   for(t in 2:nyears){
#     psi.fs[t] <- mean(z[1:nsites, t])
#   }
#   
#   
#   # (2) GoF computation part of code
#   # (based on posterior predictive distributions)
#   # --------------------------------------------
#   # Draw a replicate data set under the fitted model
#   for (i in 1:obs){
#     yrep[i] ~ dbern(z[site[i], year[i]] * p[i])
#   }
#   
#   # (2a) Computations for the GoF of the open part of the model
#   # (based on number of state transitions)
#   # ----------------------------------------------------------
#   # Compute observed z matrix for observed and replicated data
#   for (i in 1:nsites){
#     for (t in 1:nyear){
#       zobs[i,t] <- max(y[site == i & year == t])
#       zrep[i,t] <- max(yrep[site == i & year == t]) 
#     }
#     
#     # Identify extinctions, persistence, colonization and non-colonizations
#     for (t in 2:nyears){
#       # ... for observed data
#       ext[i,(t-1)] <- (zobs[i,t-1] == 1) * (zobs[i,t] == 0)
#       nonext[i,(t-1)] <- (zobs[i,t-1] == 1) * (zobs[i,t] == 1)
#       colo[i,(t-1)] <- (zobs[i,t-1] == 0) * (zobs[i,t] == 1)
#       noncolo[i,(t-1)] <- (zobs[i,t-1] == 0) * (zobs[i,t] == 0)
#       # ... for replicated data
#       extrep[i,(t-1)] <- (zobsrep[i,t-1] == 1) * (zobsrep[i,t] == 0)
#       nonextrep[i,(t-1)] <- (zobsrep[i,t-1] == 1) * (zobsrep[i,t] == 1)
#       colorep[i,(t-1)] <- (zobsrep[i,t-1] == 0) * (zobsrep[i,t] == 1)
#       noncolorep[i,(t-1)] <- (zobsrep[i,t-1] == 0) * (zobsrep[i,t] == 0)
#     }
#   }
#   
#   # Tally up number of transitions and put into a matrix for each year
#   for(t in 1:(nyears-1)){
#     # ... for observed data
#     tm[1,1,t] <- sum(noncolo[1:nsites,t]) # transition mat for obs. data
#     tm[1,2,t] <- sum(colo[1:nsites,t])
#     tm[2,1,t] <- sum(ext[1:nsites,t])
#     tm[2,2,t] <- sum(nonext[1:nsites,t])
#     # ... for replicated data
#     tmrep[1,1,t] <- sum(noncolorep[1:nsites,t]) # transition mat for rep. data
#     tmrep[1,2,t] <- sum(colorep[1:nsites,t])
#     tmrep[2,1,t] <- sum(extrep[1:nsites,t])
#     tmrep[2,2,t] <- sum(nonextrep[1:nsites,t])
#   }
#   
#   # Compute expected numbers of transitions under the model
#   # Probability of each individual transition
#   for(i in 1:nsites){
#     for(t in 1:(nyear-1)){
#       noncolo.exp[i,t] <- (1-eps[i,t]) * (1-gamma[i,t])
#       colo.exp[i,t] <- (1-eps[i,t]) * gamma[i,t]
#       ext.exp[i,t] <- psi[i,t] * (1-phi[i,t])
#       nonext.exp[i,t] <- psi[i,t] * phi[i,t]
#     }
#   }
#   # Sum up over sites to obtain the expected number of those transitions
#   for(t in 1:(nyear-1)){
#     Etm[1,1,t] <- sum(noncolo.exp[1:nsite,t])
#     Etm[1,2,t] <- sum(colo.exp[1:nsite,t])
#     Etm[2,1,t] <- sum(ext.exp[1:nsite,t])
#     Etm[2,2,t] <- sum(nonext.exp[1:nsite,t])
#   }
#   
#   # Compute Chi-square discrepancy ~~~ see Errata 2021-10-09
#   for(t in 1:(nyear-1)){
#     # ... for observed data
#     x2Open[1,1,t] <- pow((tm[1,1,t] - Etm[1,1,t]), 2) / (Etm[1,1,t]+e)
#     x2Open[1,2,t] <- pow((tm[1,2,t] - Etm[1,2,t]), 2) / (Etm[1,2,t]+e)
#     x2Open[2,1,t] <- pow((tm[2,1,t] - Etm[2,1,t]), 2) / (Etm[2,1,t]+e)
#     x2Open[2,2,t] <- pow((tm[2,2,t] - Etm[2,2,t]), 2) / (Etm[2,2,t]+e)
#     # ... for replicated data
#     x2repOpen[1,1,t] <- pow((tmrep[1,1,t]-Etm[1,1,t]),2)/(Etm[1,1,t]+e)
#     x2repOpen[1,2,t] <- pow((tmrep[1,2,t]-Etm[1,2,t]),2)/(Etm[1,2,t]+e)
#     x2repOpen[2,1,t] <- pow((tmrep[2,1,t]-Etm[2,1,t]),2)/(Etm[2,1,t]+e)
#     x2repOpen[2,2,t] <- pow((tmrep[2,2,t]-Etm[2,2,t]),2)/(Etm[2,2,t]+e)
#   }
#   
#   # Add up overall test statistic and compute fit stat ratio (open part)
#   for(i in 1:2){
#     for(j in 1:2){
#       Chi2OpenSum[i,j] <- sum(x2Open[i, j, 1:(nyear-1)])       # Chisq. statistic for observed data
#       Chi2repOpenSum[i,j] <- sum(x2repOpen[i, j, 1:(nyear-1)]) # Chisq. statistic for replicated data  
#     }
#   }
#   Chi2Open <- sum(Chi2OpenSum[1:2, 1:2])
#   Chi2repOpen <- sum(Chi2repOpenSum[1:2, 1:2])
#   Chi2ratioOpen <- Chi2Open / Chi2repOpen
#   
#   
#   # (2b) Computations for the GoF of the closed part of the model
#   # (based on the number of times detected, i.e., detection freqiencies)
#   # --------------------------------------------------------------------
#   # Compute detection frequencies for observed and replicated data
#   for (i in 1:nsite){
#     for (t in 1:nyear){
#       # Det. frequencies for observed and replicated data
#       detfreq[i,t] <- sum(y[i,1:nsurvey,t])
#       detfreqrep[i,t] <- sum(yrep[i,1:nsurvey,t])
#       
#       # Expected detection frequencies under the model
#       for (j in 1:nsurvey){
#         tmp[i,j,t] <- z[i, t] * p[i,j,t]
#       }
#       E[i,t] <- sum(tmp[i,1:nsurvey,t])     # Expected number of detections
#       
#       # Chi-square and Freeman-Tukey discrepancy measures
#       # ..... for actual data set
#       x2Closed[i,t] <- pow((detfreq[i,t] - E[i,t]),2) / (E[i,t]+e)
#       ftClosed[i,t] <- pow((sqrt(detfreq[i,t]) - sqrt(E[i,t])),2)
#       
#       # ..... for replicated data set
#       x2repClosed[i,t] <- pow((detfreqrep[i,t] - E[i,t]),2) / (E[i,t]+e)
#       ftrepClosed[i,t] <- pow((sqrt(detfreqrep[i,t]) - sqrt(E[i,t])),2)
#     }
#   }
#   
#   # Add up Chi-square and FT discrepancies and compute fit stat ratio (closed part)
#   Chi2Closed <- sum(x2Closed[1:nsite, 1:nyear])
#   FTClosed <- sum(ftClosed[1:nsite, 1:nyear])
#   Chi2repClosed <- sum(x2repClosed[1:nsite, 1:nyear])
#   FTrepClosed <- sum(ftrepClosed[1:nsite, 1:nyear])
#   Chi2ratioClosed <- Chi2Closed / Chi2repClosed
#   FTratioClosed <- FTClosed / FTrepClosed
# })
# 
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##
# ## Subsection: Multi-species models
# ##
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# dynOccMSOM <- nimbleCode({
#   
#   ## i = sites; j = reps; t = years; k = species (4-d array)
#   ## We want to accommodate a community-level hyperprior
#   
#   ## Priors
#   ## Occupancy coefficients
#   for(k in 1:nspec){
#     
#     beta0[k] ~ dnorm(0, 0.01)
#     beta1[k] ~ dnorm(0, 0.01)
#     beta2[k] ~ dnorm(0, 0.01)
#     
#     ## Detection coefficients
#     alpha0[k] ~ dnorm(0, 0.01)
#     alpha1[k] ~ dnorm(0, 0.01)
#     alpha2[k] ~ dnorm(0, 0.01)
#     
#     ## Persistence coefficients
#     eps0[k] ~ dnorm(0, 0.01)
#     eps1[k] ~ dnorm(0, 0.01)
#     eps2[k] ~ dnorm(0, 0.01)
#     eps3[k] ~ dnorm(0, 0.01)
#     
#     ## Colonization coefficients
#     gamma0[k] ~ dnorm(0, 0.01)
#     gamma1[k] ~ dnorm(0, 0.01)
#     gamma2[k] ~ dnorm(0, 0.01)
#     gamma3[k] ~ dnorm(0, 0.01)
#   }
#   
#   ## Hyperpriors
#   mu.lpsi1 ~ dnorm(0, 0.01)
#   tau.lpsi ~ pow(sd.lpsi, -2)
#   sd.lpsi ~ dunif(0, 2)
#   mu.lp ~ dnorm(0, 0.1)
#   mu.eps ~ dnorm(0, 0.01)
#   mu.gamma ~ dnorm(0, 0.01)
#   
#   
#   ## Likelihood
#   ## Ecological submodel
#   for(i in 1:nsites){
#     
#     ## Predictors for intial occ
#     logit(psi1[i]) <- lpsi1[k] + beta1 * bclim[i, 1] + beta2 * cc[i,1]
#     
#     ## Initial occ
#     z[i, 1] ~ dbern(psi1[i])
#     
#     ## State process model (colex)
#     for(t in 2:nyears){
#       
#       ## Colonization
#       logit(gamma[i, t-1]) <- gamma0 + gamma1 * anom[i, t-1] + gamma2 * fire[i, t-1] + gamma3 * anom[i, t-1] * fire[i, t-1]
#       
#       ## Persistence
#       logit(eps[i, t-1]) <- eps0 + eps1 * anom[i, t-1] + eps2 * fire[i, t-1] + eps3 * anom[i, t-1] * fire[i, t-1]
#       
#       ## Latent state
#       z[i,t] ~ dbern(z[i, t-1] * (1 - eps[i, t-1]) + (1-z[i,t-1]) * gamma[i, t-1])
#     }
#   }
#   
#   ## Observation submodel
#   for(i in 1:nsites){
#     for(j in 1:nreps){
#       for(t in 1:nyears){
#         logit(p[i,j,t]) <- alpha0 + alpha1 * eff.hrs[i,j,t] + alpha2 * eff.jday[i,j,t]
#         y[i,j,t] ~ dbern(z[i,t] * p[i,j,t])
#       }
#     }
#   }
#   
#   ## Derived parameters
#   for(i in 1:nsites){
#     for(t in 1:(nyears-1)){
#       phi[i, t] <- 1 - eps[i, t]
#     }
#   }
#   
#   ## Mean occupancy
#   psi.fs[1] <- mean(psi1[1:nsites])
#   for(t in 2:nyears){
#     psi.fs[t] <- mean(z[1:nsites, t])
#   }
#   
# })
# 
# 
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##
# ## Subsection: dynOcc MSOM
# ##
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# dynOccMSOM <- nimbleCode({
#   
#   ## These will need to be looped
#   ## for each species so we can 
#   ## accommodate species specific effects
#   
#   ## The community-level effect should be 
#   ## a common prior for the overall effect
#   ## of the community
#   for(k in 1:nspec){
#     ## Priors (mu and tau governed by hyperpriors)
#     ## Occupancy coefficients
#     beta0[k] ~ dnorm(mu.beta0, tau.beta0)
#     beta1[k] ~ dnorm(mu.beta1, tau.beta1)
#     beta2[k] ~ dnorm(mu.beta2, tau.beta2)
#     
#     ## Detection coefficients
#     ## Add some more restrictive priors
#     ## after some PPCs
#     alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0)
#     alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)
#     alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2)
#     
#     ## Persistence coefficients
#     eps0[k] ~ dnorm(mu.eps0, tau.eps0)
#     eps1[k] ~ dnorm(mu.eps1, tau.eps1)
#     eps2[k] ~ dnorm(mu.eps2, tau.eps2)
#     eps3[k] ~ dnorm(mu.eps3, tau.eps3)
#     
#     ## Colonization coefficients
#     gamma0[k] ~ dnorm(mu.gamma0, tau.gamma0)
#     gamma1[k] ~ dnorm(mu.gamma1, tau.gamma1)
#     gamma2[k] ~ dnorm(mu.gamma2, tau.gamma2)
#     gamma3[k] ~ dnorm(mu.gamma3, tau.gamma3)
#   }
#   
#   ## Hyperpriors
#   
#   ## Occupancy
#   mu.beta0 ~ dnorm(0, 0.01)
#   tau.beta0 <- pow(sd.beta0, -2)
#   sd.beta0 ~ dunif(0, 2)
#   
#   mu.beta1 ~ dnorm(0, 0.01)
#   tau.beta1 <- pow(sd.beta1, -2)
#   sd.beta1 ~ dunif(0, 2)
#   
#   mu.beta2 ~ dnorm(0, 0.01)
#   tau.beta2 <- pow(sd.beta2, -2)
#   sd.beta2 ~ dunif(0, 2)
#   
#   ## Detection
#   mu.alpha0 ~ dnorm(0, 0.01)
#   tau.alpha0 <- pow(sd.alpha0, -2)
#   sd.alpha0 ~ dunif(0, 2)
#   
#   mu.alpha1 ~ dnorm(0, 0.01)
#   tau.alpha1 <- pow(sd.alpha1, -2)
#   sd.alpha1 ~ dunif(0, 2)
#   
#   mu.alpha2 ~ dnorm(0, 0.01)
#   tau.alpha2 <- pow(sd.alpha2, -2)
#   sd.alpha2 ~ dunif(0, 2)
#   
#   ## Extinction
#   mu.eps0 ~ dnorm(0, 0.01)
#   tau.eps0 <- pow(sd.eps0, -2)
#   sd.eps0 ~ dunif(0, 2)
#   
#   mu.eps1 ~ dnorm(0, 0.01)
#   tau.eps1 <- pow(sd.eps1, -2)
#   sd.eps1 ~ dunif(0, 2)
#   
#   mu.eps2 ~ dnorm(0, 0.01)
#   tau.eps2 <- pow(sd.eps2, -2)
#   sd.eps2 ~ dunif(0, 2)
#   
#   mu.eps3 ~ dnorm(0, 0.01)
#   tau.eps3 <- pow(sd.eps3, -2)
#   sd.eps3 ~ dunif(0, 2)
#   
#   ## Colonization
#   mu.gamma0 ~ dnorm(0, 0.01)
#   tau.gamma0 <- pow(sd.gamma0, -2)
#   sd.gamma0 ~ dunif(0, 2)
#   
#   mu.gamma1 ~ dnorm(0, 0.01)
#   tau.gamma1 <- pow(sd.gamma1, -2)
#   sd.gamma1 ~ dunif(0, 2)
#   
#   mu.gamma2 ~ dnorm(0, 0.01)
#   tau.gamma2 <- pow(sd.gamma2, -2)
#   sd.gamma2 ~ dunif(0, 2)
#   
#   mu.gamma3 ~ dnorm(0, 0.01)
#   tau.gamma3 <- pow(sd.gamma3, -2)
#   sd.gamma3 ~ dunif(0, 2)
#   
#   
#   ## Likelihood
#   ## Ecological submodel
#   for(k in 1:nspec){
#     for(i in 1:nsites){
#       
#       ## Predictors for intial occ
#       logit(psi1[i,k]) <- beta0[k] + beta1[k] * bclim[i, 1] + beta2[k] * cc[i,1]
#       
#       ## Initial occ
#       z[i,1,k] ~ dbern(psi1[i,k])
#       
#       psi[i,1,k] <- psi1[i,k]
#       
#       ## State process model (colex)
#       for(t in 2:nyears){
#         
#         ## Colonization
#         logit(gamma[i,t-1,k]) <- gamma0[k] + gamma1[k] * anom[i, t-1] + gamma2[k] * fire[i, t-1] + gamma3[k] * anom[i, t-1] * fire[i, t-1]
#         
#         ## Persistence
#         logit(eps[i,t-1,k]) <- eps0[k] + eps1[k] * anom[i, t-1] + eps2[k] * fire[i, t-1] + eps3[k] * anom[i, t-1] * fire[i, t-1]
#         
#         ## Latent state
#         z[i,t,k] ~ dbern(z[i,t-1,k] * (1 - eps[i,t-1,k]) + (1-z[i,t-1,k]) * gamma[i,t-1,k])
#         
#         ## Derived estimate of psi
#         psi[i,t,k] <- psi[i,t-1,k]*(1-eps[i,t-1,k])+(1-psi[i,t-1,k])*gamma[i,t-1,k]
#       }
#     }
#   }
#   ## Observation submodel
#   for(k in 1:nspec){
#     for(i in 1:nsites){
#       for(j in 1:nreps){
#         for(t in 1:nyears){
#           logit(p[i,j,t,k]) <- alpha0[k] + alpha1[k] * eff.hrs[i,j,t] + alpha2[k] * eff.jday[i,j,t]
#           y[i,j,t,k] ~ dbern(z[i,t,k] * p[i,j,t,k])
#         }
#       }
#     }
#   }
#   
#   ## Derived parameters
#   for(k in 1:nspec){
#     for(i in 1:nsites){
#       for(t in 1:(nyears-1)){
#         phi[i,t,k] <- 1 - eps[i,t,k]
#       }
#     }
#   }
#   
#   ## Mean occupancy per species
#   for(k in 1:nspec){
#     psi.fs[1,k] <- mean(psi1[1:nsites,k])
#     for(t in 2:nyears){
#       psi.fs[t,k] <- mean(z[1:nsites,t,k])
#     }
#   }
#   
#   ## Mean richness across each site per year
#   for(i in 1:nsites){
#     for(t in 1:nyears){
#       richness[i,t] <- sum(z[i,t,1:nspec]) 
#     }
#   }
#   
# })
# 
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##
# ## Subsection: DCM Categorical
# ##
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# dynOccMSOMCat <- nimbleCode({
#   
#   ## These will need to be looped
#   ## for each species so we can 
#   ## accommodate species specific effects
#   
#   ## The community-level effect should be 
#   ## a common prior for the overall effect
#   ## of the community
#   for(k in 1:nspec){
#     ## Priors (mu and tau governed by hyperpriors)
#     ## Occupancy coefficients
#     beta0[k] ~ dnorm(mu.beta0, tau.beta0)
#     beta1[k] ~ dnorm(mu.beta1, tau.beta1)
#     beta2[k] ~ dnorm(mu.beta2, tau.beta2)
#     beta3[k] ~ dnorm(mu.beta3, tau.beta3)
#     beta4[k] ~ dnorm(mu.beta4, tau.beta4)
#     
#     ## Detection coefficients
#     ## Add some more restrictive priors
#     ## after some PPCs
#     alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0)
#     alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)
#     alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2)
#     
#     ## Persistence coefficients
#     eps0[k] ~ dnorm(mu.eps0, tau.eps0)
#     eps1[k] ~ dnorm(mu.eps1, tau.eps1)
#     eps2[k] ~ dnorm(mu.eps2, tau.eps2)
#     eps3[k] ~ dnorm(mu.eps3, tau.eps3)
#     eps4[k] ~ dnorm(mu.eps4, tau.eps4)
#     eps5[k] ~ dnorm(mu.eps5, tau.eps5)
#     
#     ## Colonization coefficients
#     gamma0[k] ~ dnorm(mu.gamma0, tau.gamma0)
#     gamma1[k] ~ dnorm(mu.gamma1, tau.gamma1)
#     gamma2[k] ~ dnorm(mu.gamma2, tau.gamma2)
#     gamma3[k] ~ dnorm(mu.gamma3, tau.gamma3)
#     gamma4[k] ~ dnorm(mu.gamma4, tau.gamma4)
#     gamma5[k] ~ dnorm(mu.gamma5, tau.gamma5)
#   }
#   
#   ## Cell-level raneff
#   for(c in 1:n_cells){
#     cell_det[c] ~ dnorm(0, tau.cell_det)
#   }
#   sd.cell_det ~ dunif(0, 2)
#   tau.cell_det <- pow(sd.cell_det, -2)
#   
#   ## Hyperpriors
#   
#   ## Occupancy
#   mu.beta0 ~ dnorm(0, 0.01)
#   tau.beta0 <- pow(sd.beta0, -2)
#   sd.beta0 ~ dunif(0, 2)
#   
#   mu.beta1 ~ dnorm(0, 0.01)
#   tau.beta1 <- pow(sd.beta1, -2)
#   sd.beta1 ~ dunif(0, 2)
#   
#   mu.beta2 ~ dnorm(0, 0.01)
#   tau.beta2 <- pow(sd.beta2, -2)
#   sd.beta2 ~ dunif(0, 2)
#   
#   mu.beta3 ~ dnorm(0, 0.01)
#   tau.beta3 <- pow(sd.beta2, -2)
#   sd.beta3 ~ dunif(0, 2)
#   
#   mu.beta4 ~ dnorm(0, 0.01)
#   tau.beta4 <- pow(sd.beta2, -2)
#   sd.beta4 ~ dunif(0, 2)
#   
#   ## Detection
#   mu.alpha0 ~ dnorm(0, 0.01)
#   tau.alpha0 <- pow(sd.alpha0, -2)
#   sd.alpha0 ~ dunif(0, 2)
#   
#   mu.alpha1 ~ dnorm(0, 0.01)
#   tau.alpha1 <- pow(sd.alpha1, -2)
#   sd.alpha1 ~ dunif(0, 2)
#   
#   mu.alpha2 ~ dnorm(0, 0.01)
#   tau.alpha2 <- pow(sd.alpha2, -2)
#   sd.alpha2 ~ dunif(0, 2)
#   
#   ## Extinction
#   mu.eps0 ~ dnorm(0, 0.01)
#   tau.eps0 <- pow(sd.eps0, -2)
#   sd.eps0 ~ dunif(0, 2)
#   
#   mu.eps1 ~ dnorm(0, 0.01)
#   tau.eps1 <- pow(sd.eps1, -2)
#   sd.eps1 ~ dunif(0, 2)
#   
#   mu.eps2 ~ dnorm(0, 0.01)
#   tau.eps2 <- pow(sd.eps2, -2)
#   sd.eps2 ~ dunif(0, 2)
#   
#   mu.eps3 ~ dnorm(0, 0.01)
#   tau.eps3 <- pow(sd.eps3, -2)
#   sd.eps3 ~ dunif(0, 2)
#   
#   mu.eps4 ~ dnorm(0, 0.01)
#   tau.eps4 <- pow(sd.eps4, -2)
#   sd.eps4 ~ dunif(0, 2)
#   
#   mu.eps5 ~ dnorm(0, 0.01)
#   tau.eps5 <- pow(sd.eps5, -2)
#   sd.eps5 ~ dunif(0, 2)
#   
#   ## Colonization
#   mu.gamma0 ~ dnorm(0, 0.01)
#   tau.gamma0 <- pow(sd.gamma0, -2)
#   sd.gamma0 ~ dunif(0, 2)
#   
#   mu.gamma1 ~ dnorm(0, 0.01)
#   tau.gamma1 <- pow(sd.gamma1, -2)
#   sd.gamma1 ~ dunif(0, 2)
#   
#   mu.gamma2 ~ dnorm(0, 0.01)
#   tau.gamma2 <- pow(sd.gamma2, -2)
#   sd.gamma2 ~ dunif(0, 2)
#   
#   mu.gamma3 ~ dnorm(0, 0.01)
#   tau.gamma3 <- pow(sd.gamma3, -2)
#   sd.gamma3 ~ dunif(0, 2)
#   
#   mu.gamma4 ~ dnorm(0, 0.01)
#   tau.gamma4 <- pow(sd.gamma4, -2)
#   sd.gamma4 ~ dunif(0, 2)
#   
#   mu.gamma5 ~ dnorm(0, 0.01)
#   tau.gamma5 <- pow(sd.gamma5, -2)
#   sd.gamma5 ~ dunif(0, 2)
#   
#   
#   ## Likelihood
#   ## Ecological submodel
#   for(k in 1:nspec){
#     for(i in 1:nsites){
#       
#       ## Predictors for intial occ
#       logit(psi1[i,k]) <- beta0[k] + 
#         beta1[k] * bclim[i, 1] + 
#         beta2[k] * cc[i,1] + 
#         beta3[k] * fire_init[i,2] + 
#         beta4[k] * fire_init[i,3]
#       
#       ## Initial occ
#       z[i,1,k] ~ dbern(psi1[i,k])
#       
#       psi[i,1,k] <- psi1[i,k]
#       
#       ## State process model (colex)
#       for(t in 2:nyears){
#         
#         ## Colonization
#         logit(gamma[i,t-1,k]) <- gamma0[k] + 
#           gamma1[k] * anom[i, t-1] + 
#           gamma2[k] * fire[i, t-1, 1] + 
#           gamma3[k] * fire[i,t-1,2] + 
#           gamma4[k] * anom[i,t-1] * fire[i,t-1,1] + 
#           gamma5[k] * anom[i,t-1] * fire[i,t-1,2]
#         
#         ## Persistence
#         logit(eps[i,t-1,k]) <- eps0[k] + 
#           eps1[k] * anom[i, t-1] + 
#           eps2[k] * fire[i, t-1, 1] + 
#           eps3[k] * fire[i,t-1,2] + 
#           eps4[k] * anom[i,t-1] * fire[i,t-1,1] + 
#           eps5[k] * anom[i,t-1] * fire[i,t-1,2]
#         
#         ## NEED to DOUBLE CHECK LIKELIHOOD FOR COLONIZATION AND EXTINCTION
#         ## I think below has been corrected
#         ## Latent state
#         z[i,t,k] ~ dbern(z[i,t-1,k] * (1-eps[i,t-1,k]) + 
#                            (1-z[i,t-1,k]) * gamma[i,t-1,k])
#         
#         ## Derived estimate of psi
#         psi[i,t,k] <- psi[i,t-1,k] * (1-eps[i,t-1,k]) + 
#           (1-psi[i,t-1,k]) * gamma[i,t-1,k]
#       }
#     }
#   }
#   ## Observation submodel
#   for(k in 1:nspec){
#     for(i in 1:nsites){
#       for(j in 1:nreps){
#         for(t in 1:nyears){
#           
#           logit(p[i,j,t,k]) <- alpha0[k] + 
#             alpha1[k] * eff.hrs[i,j,t] + 
#             alpha2[k] * eff.jday[i,j,t]
#           
#           y[i,j,t,k] ~ dbern(z[i,t,k] * p[i,j,t,k])
#         }
#       }
#     }
#   }
#   
#   ## Derived parameters
#   for(k in 1:nspec){
#     for(i in 1:nsites){
#       for(t in 1:(nyears-1)){
#         phi[i,t,k] <- 1 - eps[i,t,k]
#       }
#     }
#   }
#   
#   ## Mean occupancy per species
#   for(k in 1:nspec){
#     psi.fs[1,k] <- mean(psi1[1:nsites,k])
#     for(t in 2:nyears){
#       psi.fs[t,k] <- mean(z[1:nsites,t,k])
#     }
#   }
#   
#   ## Mean richness across each site per year
#   for(i in 1:nsites){
#     for(t in 1:nyears){
#       richness[i,t] <- sum(z[i,t,1:nspec]) 
#     }
#   }
#   
# })
# 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: DCM Categorical Vertical structure
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DCMCatVert <- nimbleCode({
  
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
    
    ## Persistence coefficients
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
  }
  
  ## Cell-level raneff
  for(c in 1:n_cells){
    cell_det[c] ~ dnorm(0, tau.cell_det)
  }
  sd.cell_det ~ dunif(0, 2)
  tau.cell_det <- pow(sd.cell_det, -2)
  
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
  
  ## Persistence hyperpriors
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
        beta2[k] * bprec[i] +
        beta3[k] * cc[i] +
        beta4[k] * trendt[i] +
        beta5[k] * trendtmin[i] +
        beta6[k] * trendp[i]
      
      z[i,1,k] ~ dbern(psi1[i,k])
      psi[i,1,k] <- psi1[i,k]
      
      ## State process over years
      for(t in 2:nyears){
        ## Colonization
        logit(gamma[i,t-1,k]) <- gamma0[k] +
          gamma1[k] * tanom[i, t-1] +
          gamma2[k] * panom[i,t-1] +
          gamma3[k] * fire[i, t-1, 1] +
          gamma4[k] * fire[i,t-1,2] +
          gamma5[k] * tanom[i,t-1] * fire[i,t-1,1] +
          gamma6[k] * tanom[i,t-1] * fire[i,t-1,2] +
          gamma7[k] * panom[i,t-1] * fire[i,t-1,1] +
          gamma8[k] * panom[i,t-1] * fire[i,t-1,2] +
          gamma9[k] * trendt[i] + 
          gamma10[k] * trendtmin[i] +
          gamma11[k] * trendp[i]
        
        ## Persistence
        logit(eps[i,t-1,k]) <- eps0[k] +
          eps1[k] * tanom[i,t-1] +
          eps2[k] * panom[i,t-1] +
          eps3[k] * fire[i,t-1,1] +
          eps4[k] * fire[i,t-1,2] +
          eps5[k] * tanom[i,t-1] * fire[i,t-1,1] +
          eps6[k] * tanom[i,t-1] * fire[i,t-1,2] +
          eps7[k] * panom[i,t-1] * fire[i,t-1,1] +
          eps8[k] * panom[i,t-1] * fire[i,t-1,2] +
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
        cell_det[cell_id[site_obs[v]]]
      
      y[v,k] ~ dbern(z[site_obs[v], year_obs[v], k] * p_obs[v,k])
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


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: DCM Categorical Vertical structure
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DCMCatVertNP <- nimbleCode({
  
  ##----------------------
  ## Species-level priors
  ##----------------------
  for(k in 1:nspec){
    ## Occupancy coefficients
    beta0[k] ~ dnorm(0, 0.01)
    beta1[k] ~ dnorm(0, 0.01)
    beta2[k] ~ dnorm(0, 0.01)
    beta3[k] ~ dnorm(0, 0.01)
    beta4[k] ~ dnorm(0, 0.01)
    beta5[k] ~ dnorm(0, 0.01)
    beta6[k] ~ dnorm(0, 0.01)
    
    ## Detection coefficients
    alpha0[k] ~ dnorm(0, 0.01)
    alpha1[k] ~ dnorm(0, 0.01)
    alpha2[k] ~ dnorm(0, 0.01)
    alpha3[k] ~ dnorm(0, 0.01)
    alpha4[k] ~ dnorm(0, 0.01)
    
    ## Persistence coefficients
    eps0[k] ~ dnorm(0, 0.01)
    eps1[k] ~ dnorm(0, 0.01)
    eps2[k] ~ dnorm(0, 0.01)
    eps3[k] ~ dnorm(0, 0.01)
    eps4[k] ~ dnorm(0, 0.01)
    eps5[k] ~ dnorm(0, 0.01)
    eps6[k] ~ dnorm(0, 0.01)
    eps7[k] ~ dnorm(0, 0.01)
    eps8[k] ~ dnorm(0, 0.01)
    eps9[k] ~ dnorm(0, 0.01)
    eps10[k] ~ dnorm(0, 0.01)
    eps11[k] ~ dnorm(0, 0.01)
    
    ## Colonization coefficients
    gamma0[k] ~ dnorm(0, 0.01)
    gamma1[k] ~ dnorm(0, 0.01)
    gamma2[k] ~ dnorm(0, 0.01)
    gamma3[k] ~ dnorm(0, 0.01)
    gamma4[k] ~ dnorm(0, 0.01)
    gamma5[k] ~ dnorm(0, 0.01)
    gamma6[k] ~ dnorm(0, 0.01)
    gamma7[k] ~ dnorm(0, 0.01)
    gamma8[k] ~ dnorm(0, 0.01)
    gamma9[k] ~ dnorm(0, 0.01)
    gamma10[k] ~ dnorm(0, 0.01)
    gamma11[k] ~ dnorm(0, 0.01)
  }
  
  ## Cell-level raneff
  for(c in 1:n_cells){
    cell_det[c] ~ dnorm(0, tau.cell_det)
  }
  sd.cell_det ~ dunif(0, 2)
  tau.cell_det <- pow(sd.cell_det, -2)
  
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
        beta4[k] * trendt[i] +
        beta5[k] * trendtmin[i] +
        beta6[k] * trendp[i]
      
      z[i,1,k] ~ dbern(psi1[i,k])
      psi[i,1,k] <- psi1[i,k]
      
      ## State process over years
      for(t in 2:nyears){
        ## Colonization
        logit(gamma[i,t-1,k]) <- gamma0[k] +
          gamma1[k] * tanom[i, t-1] +
          gamma2[k] * panom[i,t-1] +
          gamma3[k] * fire[i, t-1, 1] +
          gamma4[k] * fire[i,t-1,2] +
          gamma5[k] * tanom[i,t-1] * fire[i,t-1,1] +
          gamma6[k] * tanom[i,t-1] * fire[i,t-1,2] +
          gamma7[k] * panom[i,t-1] * fire[i,t-1,1] +
          gamma8[k] * panom[i,t-1] * fire[i,t-1,2] +
          gamma9[k] * trendt[i] + 
          gamma10[k] * trendtmin[i] +
          gamma11[k] * trendp[i]
        
        ## Persistence
        logit(eps[i,t-1,k]) <- eps0[k] +
          eps1[k] * tanom[i,t-1] +
          eps2[k] * panom[i,t-1] +
          eps3[k] * fire[i,t-1,1] +
          eps4[k] * fire[i,t-1,2] +
          eps5[k] * tanom[i,t-1] * fire[i,t-1,1] +
          eps6[k] * tanom[i,t-1] * fire[i,t-1,2] +
          eps7[k] * panom[i,t-1] * fire[i,t-1,1] +
          eps8[k] * panom[i,t-1] * fire[i,t-1,2] +
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
        cell_det[cell_id[site_obs[v]]]
      
      y[v,k] ~ dbern(z[site_obs[v], year_obs[v], k] * p_obs[v,k])
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
  
  nspec <- bdata$nspec
  
  list(z = z_init,
       
       ## Species-level
       beta0 = rnorm(nspec, 0, 1),
       beta1 = rnorm(nspec, 0, 1),
       beta2 = rnorm(nspec, 0, 1),
       beta3 = rnorm(nspec, 0, 1),
       beta4 = rnorm(nspec, 0, 1),
       beta5 = rnorm(nspec, 0, 1),
       alpha0 = rnorm(nspec, 0, 1),
       alpha1 = rnorm(nspec, 0, 1),
       alpha2 = rnorm(nspec, 0, 1),
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
       eps10 = rnorm(nspec, 0, 1))
  
  # ## Community level
  # mu.beta0 = 0,  sd.beta0 = 1,
  # mu.beta1 = 0,  sd.beta1 = 1,
  # mu.beta2 = 0,  sd.beta2 = 1,
  # mu.beta3 = 0,  sd.beta3 = 1,
  # #mu.beta4 = 0,  sd.beta4 = 1,
  # 
  # mu.alpha0 = 0, sd.alpha0 = 1,
  # mu.alpha1 = 0, sd.alpha1 = 1,
  # mu.alpha2 = 0, sd.alpha2 = 1,
  # 
  # mu.eps0 = 0,  sd.eps0 = 1,
  # mu.eps1 = 0,  sd.eps1 = 1,
  # mu.eps2 = 0,  sd.eps2 = 1,
  # mu.eps3 = 0,  sd.eps3 = 1,
  # mu.eps4 = 0,  sd.eps4 = 1,
  # mu.eps5 = 0,  sd.eps5 = 1,
  # 
  # mu.gamma0 = 0, sd.gamma0 = 1,
  # mu.gamma1 = 0, sd.gamma1 = 1,
  # mu.gamma2 = 0, sd.gamma2 = 1,
  # mu.gamma3 = 0, sd.gamma3 = 1,
  # mu.gamma4 = 0, sd.gamma4 = 1,
  # mu.gamma5 = 0, sd.gamma5 = 1)
}

# Parameters monitored
params <- c(paste0("alpha",
                  0:4),
            paste0("beta",
                  0:6),
            paste0("gamma",
                  0:11),
            paste0("eps",
                  0:11),
            ## Community
            paste0("mu.alpha",
                   0:4),
            paste0("mu.beta",
                   0:6),
            paste0("mu.gamma",
                   0:11),
            paste0("mu.eps",
                  0:11),
            ## Derived params
            "richness",
            "phi",
            "psi"
            )


## Load data
bdata <- readRDS(here("Data/Occ_Data/LABU_Multi_Test.RDS"))  


## MCMC settings
nburnin = 2000
niter = 20000
nthin = 20
nchain = 3

## Run the sampler
samples <- nimbleMCMC(
  code = dynOcc,
  constants = bdata,
  inits = inits,
  monitors = params,
  niter = niter,
  nburnin = nburnin,
  thin = nthin,
  nchains = nchain)
