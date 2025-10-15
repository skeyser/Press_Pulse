## -------------------------------------------------------------
##
## Script name: Prior predictive checks for multi-season occupancy
## model
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2025-09-26
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

## -------------------------------------------------------------

## Values of interest
# Prior predictive check for dynamic occupancy model
prior_pred_check <- function(nsites = 100, nyears = 4, nsurveys = 3,
                             nsim = 1000) {
  
  ## Number of simulations
  nsims <- nsim
  
  ## Storage for simulated values
  occ_prob <- matrix(NA, nsims, nyears)  # Occupancy probability
  det_prob <- numeric(nsims)             # Detection probability
  col_prob <- numeric(nsims)             # Colonization probability
  ext_prob <- numeric(nsims)             # Extinction probability
  
  ## Simulate from priors
  for(sim in 1:nsims) {
    ## Draw parameters from priors (sqrt)
    beta0 <- rnorm(1, 0, 10)     # Occupancy intercept
    beta1 <- rnorm(1, 0, 10)     # Occupancy climate effect
    beta2 <- rnorm(1, 0, 10)     # Occupancy canopy effect
    
    alpha0 <- rnorm(1, 0, 3.16)    # Detection intercept
    alpha1 <- rnorm(1, 0, 3.16)    # Detection hours effect
    alpha2 <- rnorm(1, 0, 3.16)    # Detection date effect
    alpha3 <- rnorm(1, 0, 3.16)    # Detection date^2
    
    gamma0 <- rnorm(1, 0, 10)    # Colonization intercept
    gamma1 <- rnorm(1, 0, 10)    # Colonization temp effect
    gamma2 <- rnorm(1, 0, 10)    # Colonization fire effect
    gamma3 <- rnorm(1, 0, 10)    # Colonization interaction
    
    eps0 <- rnorm(1, 0, 10)      # Extinction intercept
    eps1 <- rnorm(1, 0, 10)      # Extinction temp effect
    eps2 <- rnorm(1, 0, 10)      # Extinction fire effect
    eps3 <- rnorm(1, 0, 10)      # Extinction interaction
    
    # Simulate states and observations
    z <- matrix(NA, nsites, nyears)
    y <- array(NA, dim = c(nsites, nsurveys, nyears))
    
    # Initial occupancy
    psi1 <- plogis(beta0 + beta1 * rnorm(nsites) + beta2 * rnorm(nsites))
    z[,1] <- rbinom(nsites, 1, psi1)
    occ_prob[sim, 1] <- mean(psi1)
    
    ## Subsequent years
    for(t in 2:nyears) {
      ## Simulate covariates
      temp <- rnorm(nsites)
      fire <- rnorm(nsites)
      
      ## Calculate transition probabilities
      gamma <- plogis(gamma0 + gamma1*temp + gamma2*fire + gamma3*temp*fire)
      eps <- plogis(eps0 + eps1*temp + eps2*fire + eps3*temp*fire)
      
      # State process
      for(i in 1:nsites) {
        if(z[i,t-1] == 1) {
          z[i,t] <- rbinom(1, 1, 1-eps[i])
        } else {
          z[i,t] <- rbinom(1, 1, gamma[i])
        }
      }
      
      occ_prob[sim, t] <- mean(z[,t])
    }
    
    # Store average transition probabilities
    col_prob[sim] <- mean(gamma)
    ext_prob[sim] <- mean(eps)
    
    # Detection process
    p <- plogis(alpha0 + alpha1*rnorm(1) + alpha2*rnorm(1) + alpha3*(rnorm(1)^2))
    det_prob[sim] <- p
  }
  
  # Create plots
  par(mfrow=c(2,2))
  
  # Plot 1: Occupancy trajectories
  matplot(1:nyears, t(occ_prob), type='l', col=rgb(0,0,0,0.1),
          xlab="Year", ylab="Occupancy probability",
          main="Simulated occupancy trajectories")
  lines(1:nyears, colMeans(occ_prob), col="red", lwd=2)
  
  # Plot 2: Detection probability
  hist(det_prob, breaks=30, main="Detection probability",
       xlab="Probability", col="lightblue")
  abline(v=mean(det_prob), col="red", lwd=2)
  
  # Plot 3: Colonization probability
  hist(col_prob, breaks=30, main="Colonization probability",
       xlab="Probability", col="lightgreen")
  abline(v=mean(col_prob), col="red", lwd=2)
  
  # Plot 4: Extinction probability
  hist(ext_prob, breaks=30, main="Extinction probability",
       xlab="Probability", col="pink")
  abline(v=mean(ext_prob), col="red", lwd=2)
  
  # Return summary statistics
  list(
    occ_prob_summary = apply(occ_prob, 2, quantile, probs=c(0.025, 0.5, 0.975)),
    det_prob_summary = quantile(det_prob, probs=c(0.025, 0.5, 0.975)),
    col_prob_summary = quantile(col_prob, probs=c(0.025, 0.5, 0.975)),
    ext_prob_summary = quantile(ext_prob, probs=c(0.025, 0.5, 0.975))
  )
}

# Run the prior predictive check
set.seed(123)
results <- prior_pred_check(nsites = 2095, nsurveys = 5, nyears = 4)

# Print summaries
print("Prior predictive summaries:")
print("Occupancy probability by year:")
print(results$occ_prob_summary)
print("\nDetection probability:")
print(results$det_prob_summary)
print("\nColonization probability:")
print(results$col_prob_summary)
print("\nExtinction probability:")
print(results$ext_prob_summary)


## MSOM flavor

## Values of interest
# Prior predictive check for dynamic occupancy model
prior_pred_check_nimble <- function(nsites = 100, nyears = 4, nsurveys = 5,
                                    nspec = 2, nsim = 1000) {
  
  ## Storage for simulated values
  occ_prob <- array(NA, dim = c(nsim, nyears, nspec))
  det_prob <- matrix(NA, nsim, nspec)
  col_prob <- matrix(NA, nsim, nspec)
  ext_prob <- matrix(NA, nsim, nspec)
  richness <- array(NA, dim = c(nsim, nyears))
  
  ## Simulate predictors
  bclim <- rnorm(nsites)  # baseline climate
  cc <- rnorm(nsites)     # canopy cover
  
  # Time-varying predictors
  anom <- array(rnorm(nsites * nyears), dim = c(nsites, nyears))  # temperature anomaly
  fire <- array(rnorm(nsites * nyears), dim = c(nsites, nyears))  # fire
  
  # Detection covariates (assuming these are fixed across years)
  eff.hrs <- array(rnorm(nsites * nsurveys * nyears), 
                   dim = c(nsites, nsurveys, nyears))
  eff.jday <- array(rnorm(nsites * nsurveys * nyears), 
                    dim = c(nsites, nsurveys, nyears))
  
  for(sim in 1:nsim) {
    # Generate community hyperparameters
    # Occupancy hyperparameters
    mu.beta0 <- rnorm(1, 0, sqrt(1/0.01))
    sd.beta0 <- runif(1, 0, 2)
    mu.beta1 <- rnorm(1, 0, sqrt(1/0.01))
    sd.beta1 <- runif(1, 0, 2)
    mu.beta2 <- rnorm(1, 0, sqrt(1/0.01))
    sd.beta2 <- runif(1, 0, 2)
    
    # Detection hyperparameters
    mu.alpha0 <- rnorm(1, 0, sqrt(1/0.01))
    sd.alpha0 <- runif(1, 0, 2)
    mu.alpha1 <- rnorm(1, 0, sqrt(1/0.01))
    sd.alpha1 <- runif(1, 0, 2)
    mu.alpha2 <- rnorm(1, 0, sqrt(1/0.01))
    sd.alpha2 <- runif(1, 0, 2)
    
    # Colonization hyperparameters
    mu.gamma0 <- rnorm(1, 0, sqrt(1/0.01))
    sd.gamma0 <- runif(1, 0, 2)
    mu.gamma1 <- rnorm(1, 0, sqrt(1/0.01))
    sd.gamma1 <- runif(1, 0, 2)
    mu.gamma2 <- rnorm(1, 0, sqrt(1/0.01))
    sd.gamma2 <- runif(1, 0, 2)
    mu.gamma3 <- rnorm(1, 0, sqrt(1/0.01))
    sd.gamma3 <- runif(1, 0, 2)
    
    # Extinction hyperparameters
    mu.eps0 <- rnorm(1, 0, sqrt(1/0.01))
    sd.eps0 <- runif(1, 0, 2)
    mu.eps1 <- rnorm(1, 0, sqrt(1/0.01))
    sd.eps1 <- runif(1, 0, 2)
    ## Strong influence on prior
    mu.eps2 <- rnorm(1, 0, sqrt(1/0.01))
    sd.eps2 <- runif(1, 0, 2)
    mu.eps3 <- rnorm(1, 0, sqrt(1/0.01))
    sd.eps3 <- runif(1, 0, 2)
    
    # Create z array for all species
    z <- array(NA, dim = c(nsites, nyears, nspec))
    
    for(spec in 1:nspec) {
      # Species-specific parameters
      beta0 <- rnorm(1, mu.beta0, sd.beta0)
      beta1 <- rnorm(1, mu.beta1, sd.beta1)
      beta2 <- rnorm(1, mu.beta2, sd.beta2)
      
      alpha0 <- rnorm(1, mu.alpha0, sd.alpha0)
      alpha1 <- rnorm(1, mu.alpha1, sd.alpha1)
      alpha2 <- rnorm(1, mu.alpha2, sd.alpha2)
      
      gamma0 <- rnorm(1, mu.gamma0, sd.gamma0)
      gamma1 <- rnorm(1, mu.gamma1, sd.gamma1)
      gamma2 <- rnorm(1, mu.gamma2, sd.gamma2)
      gamma3 <- rnorm(1, mu.gamma3, sd.gamma3)
      
      eps0 <- rnorm(1, mu.eps0, sd.eps0)
      eps1 <- rnorm(1, mu.eps1, sd.eps1)
      eps2 <- rnorm(1, mu.eps2, sd.eps2)
      eps3 <- rnorm(1, mu.eps3, sd.eps3)
      
      # Initial occupancy
      psi1 <- plogis(beta0 + beta1 * bclim + beta2 * cc)
      z[,1,spec] <- rbinom(nsites, 1, psi1)
      occ_prob[sim, 1, spec] <- mean(psi1)
      
      # Subsequent years
      for(t in 2:nyears) {
        gamma <- plogis(gamma0 + gamma1*anom[,t] + gamma2*fire[,t] + 
                          gamma3*anom[,t]*fire[,t])
        eps <- plogis(eps0 + eps1*anom[,t] + eps2*fire[,t] + 
                        eps3*anom[,t]*fire[,t])
        
        for(i in 1:nsites) {
          prob <- z[i,t-1,spec] * (1-eps[i]) + (1-z[i,t-1,spec]) * gamma[i]
          z[i,t,spec] <- rbinom(1, 1, prob)
        }
        
        occ_prob[sim, t, spec] <- mean(z[,t,spec])
      }
      
      # Store average transition probabilities
      col_prob[sim,spec] <- mean(gamma)
      ext_prob[sim,spec] <- mean(eps)
      
      # Detection probability (using fixed covariates)
      # Calculate mean detection probability across all surveys and years
      p <- array(NA, dim = c(nsites, nsurveys, nyears))
      for(t in 1:nyears) {
        for(j in 1:nsurveys) {
          p[,j,t] <- plogis(alpha0 + alpha1*eff.hrs[,j,t] + alpha2*eff.jday[,j,t])
        }
      }
      det_prob[sim,spec] <- mean(p)
    }
    
    # Calculate richness for each year
    for(t in 1:nyears) {
      richness[sim,t] <- mean(rowSums(z[,t,]))
    }
  }
  
  # Create plots
  par(mfrow=c(3,3))
  
  # 1. Community-level occupancy trajectories
  matplot(1:nyears, t(apply(occ_prob, c(1,2), mean)), type='l', 
          col=rgb(0,0,0,0.1), xlab="Year", ylab="Occupancy probability",
          main="Community occupancy trajectories")
  lines(1:nyears, apply(apply(occ_prob, c(1,2), mean), 2, mean), 
        col="red", lwd=2)
  
  # 2. Species-specific occupancy trajectories
  for(k in 1:nspec) {
    if(k == 1) {
      plot(1:nyears, colMeans(occ_prob[,,k]), type="l", col=k,
           ylim=range(occ_prob), xlab="Year", ylab="Occupancy probability",
           main="Species-specific occupancy")
    } else {
      lines(1:nyears, colMeans(occ_prob[,,k]), col=k)
    }
  }
  legend("topright", legend=paste("Species", 1:nspec), col=1:nspec, lwd=1)
  
  # 3. Expected richness
  plot(1:nyears, colMeans(richness), type="l", col="blue", lwd=2,
       xlab="Year", ylab="Species richness",
       main="Expected species richness")
  rich_ci <- apply(richness, 2, quantile, probs=c(0.025, 0.975))
  polygon(c(1:nyears, rev(1:nyears)), 
          c(rich_ci[1,], rev(rich_ci[2,])),
          col=rgb(0,0,1,0.2), border=NA)
  
  # 4. Detection probability by species
  boxplot(det_prob, main="Detection probability by species",
          xlab="Species", ylab="Probability", col="lightblue")
  
  # 5. Colonization probability by species
  boxplot(col_prob, main="Colonization probability by species",
          xlab="Species", ylab="Probability", col="lightgreen")
  
  # 6. Extinction probability by species
  boxplot(ext_prob, main="Extinction probability by species",
          xlab="Species", ylab="Probability", col="pink")
  
  # Return summary statistics
  results <- list(
    # Community-level summaries
    community_occ = apply(apply(occ_prob, c(1,2), mean), 2, 
                          quantile, probs=c(0.025, 0.5, 0.975)),
    community_det = apply(det_prob, 1, mean) %>% 
      quantile(probs=c(0.025, 0.5, 0.975)),
    community_col = apply(col_prob, 1, mean) %>% 
      quantile(probs=c(0.025, 0.5, 0.975)),
    community_ext = apply(ext_prob, 1, mean) %>% 
      quantile(probs=c(0.025, 0.5, 0.975)),
    richness_summary = apply(richness, 2, 
                             quantile, probs=c(0.025, 0.5, 0.975)),
    
    # Species-specific summaries
    species_occ = lapply(1:nspec, function(k) 
      apply(occ_prob[,,k], 2, quantile, probs=c(0.025, 0.5, 0.975))),
    species_det = lapply(1:nspec, function(k)
      quantile(det_prob[,k], probs=c(0.025, 0.5, 0.975))),
    species_col = lapply(1:nspec, function(k)
      quantile(col_prob[,k], probs=c(0.025, 0.5, 0.975))),
    species_ext = lapply(1:nspec, function(k)
      quantile(ext_prob[,k], probs=c(0.025, 0.5, 0.975)))
  )
  
  return(results)
}

# Run the prior predictive check
set.seed(123)
results <- prior_pred_check_nimble(nsites = 100, nsurveys = 5, nyears = 4, nspec = 2)

# Print summaries
print("Community-level summaries:")
print("Mean occupancy probability by year:")
print(results$community_occ)
print("\nMean detection probability:")
print(results$community_det)
print("\nMean colonization probability:")
print(results$community_col)
print("\nMean extinction probability:")
print(results$community_ext)
print("\nExpected richness by year:")
print(results$richness_summary)

print("\nSpecies-specific summaries:")
for(k in 1:nspec) {
  cat(sprintf("\nSpecies %d:\n", k))
  cat("Occupancy probability by year:\n")
  print(results$species_occ[[k]])
  cat("Detection probability:\n")
  print(results$species_det[[k]])
  cat("Colonization probability:\n")
  print(results$species_col[[k]])
  cat("Extinction probability:\n")
  print(results$species_ext[[k]])
}

