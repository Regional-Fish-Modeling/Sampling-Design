model{
  
  # Abundance model  
  
  for(i in 1:nSites){
    for(j in 1:nYears){
      N[i,j] ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- mu + trend*(j-1) + #mu + (trend + slope.ran[i])*(j-1) +
        site.ran[i] + year.ran[j] + eps[i,j] + log(siteLength[i,j]) + 
        b[1,i]*fallPrcp[i,j] + b[2,i]*fallTmean[i,j] + 
        b[3,i]*winterPrcp[i,j] + b[4,i]*winterTmean[i,j] +
        b[5,i]*springPrcp[i,j] + b[6,i]*springTmean[i,j]
    }
  }
  
  ## priors
  mu ~ dunif(0, 5)      # overall intercept
#  mu ~ dnorm(0, 0.01)      # overall intercept
  trend ~ dnorm(0, 0.1)   # linear trend
#  for(i in 1:nSites){
#    slope.ran[i] ~ dnorm(0, tau.slope)  # random slope effects (variation in trend among sites)
#  }
#  sd.slope ~ dunif(0,5)
#  tau.slope <- pow(sd.slope, -2)
  
  for(i in 1:nSites){
    site.ran[i] ~ dnorm(0,tau.site)     # random site effects
  }
  tau.site <- pow(sd.site, -2)
  sd.site ~ dunif(0,2)
  sd2.site <- pow(sd.site, 2)
  
  for (j in 1:nYears){
    year.ran[j] ~ dnorm(0, tau.year)  # Random year effect
  }
  tau.year <- pow(sd.year, -2) 
  sd.year ~ dunif(0,2)
  sd2.year <- pow(sd.year, 2)
  
  for(i in 1:nSites){
    for(j in 1:nYears){
      eps[i,j] ~ dnorm(0, tau)  # Over-dispersion
    }
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 2)
  sigma2 <- pow(sigma, 2)

  ## weather effect
  for(h in 1:nCovs){
    for(i in 1:nSites){
      b[h,i] ~ dnorm(mu.b[h], tau.b[h])
    }  
    mu.b[h] ~ dnorm(0, 0.01)  
    tau.b[h] <- pow(sigma.b[h], -2)
    sigma.b[h] ~ dunif(0, 2)
    sigma2.b[h] <- pow(sigma.b[h], 2)
  }
  
  # Detection model
  
  for(i in 1:nSites){
    for(j in 1:nYears){
      y[i,j,1] ~ dbin(p[i,j], N[i,j])
      y[i,j,2] ~ dbin(p[i,j], N[i,j]-y[i,j,1])
      y[i,j,3] ~ dbin(p[i,j], N[i,j]-y[i,j,1]-y[i,j,2])
      
      p[i,j] <- 1/(1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      lp[i,j] <- p.mu + p.b*sampday[i,j] #+ p.eps[i]
      # removed "prcp7day" for its small effect
    }
  }
  
  ## priors  
  p.mean ~ dunif(0.2,0.8)
  p.mu <- log(p.mean/(1-p.mean))
  p.b ~ dnorm(0, 0.37)
  
#  for(i in 1:nSites){
#    p.eps[i] ~ dnorm(0, p.tau)  # site-level variation
#  }
#  p.tau <- pow(p.sigma,-2)
#  p.sigma ~ dunif(0,3)
#  p.sigma2 <- pow(p.sigma, 2)
}