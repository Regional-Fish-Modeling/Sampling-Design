model{
  
  # Abundance model  
  
  for(i in 1:nSites){
    for(j in 1:nYears){
      N[i,j] ~ dpois(lambda[i,j])
      log(lambda[i,j]) <- mu + (trend + slope.ran[i])*(j-1) + site.ran[i] + year.ran[j] + eps[i,j]
    }
  }
  
  ## priors
  mu ~ dnorm(0, 0.01)      # overall intercept
  trend ~ dnorm(0, 0.01)   # linear trend
  for(i in 1:nSites){
    slope.ran[i] ~ dnorm(0, tau.slope)  # random slope effects
  }
  sd.slope ~ dunif(0,5)
  tau.slope <- pow(sd.slope, -2)
  
  ######### Random walk on spatial variance
  # the first year estimates of spatial variability
  # sd.site[1] ~ dunif(0,5)
  # tau.site[1] <- 1/(sd.site[1]^2)
  # wn[1] <- 0  # white noise
  
  # estimating spatail variance as AR1 process
  # for(j in 2:nYears){
  #  tau.site[j] <- tau.site[j-1] + wn[j]
  #  wn[j] ~ dnorm(0,tau.wn)
  #}
  #tau.wn <- 1/(sd.wn^2)
  #sd.wn ~ dunif(0,5) # first order random walk model variance parameter 
  
  for(i in 1:nSites){
    site.ran[i] ~ dnorm(0,tau.site)     # random site effects
  }
  tau.site <- pow(sd.site, -2)
  sd.site ~ dunif(0,5)
  sd2.site <- pow(sd.site, 2)
  
  
  for (j in 1:nYears){
    year.ran[j] ~ dnorm(0, tau.year)  # Random year effect
  }
  tau.year <- pow(sd.year, -2) 
  sd.year ~ dunif(0,5)
  sd2.year <- pow(sd.year, 2)
  
  for(i in 1:nSites){
    for(j in 1:nYears){
      eps[i,j] ~ dnorm(0, tau)  # Over-dispersion
    }
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 5)
  sigma2 <- pow(sigma, 2)
  
  # Detection model
  
  for(i in 1:nSites){
    for(j in 1:nYears){
      y[i,j,1] ~ dbin(p[i,j], N[i,j])
      y[i,j,2] ~ dbin(p[i,j], N[i,j]-y[i,j,1])
      y[i,j,3] ~ dbin(p[i,j], N[i,j]-y[i,j,1]-y[i,j,2])
      
      p[i,j] <- 1/(1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      lp[i,j] <- p.mu + p.b*sampday[i,j] + p.eps[i]
      # removed "prcp7day" for its small effect
    }
  }
  
  ## priors  
  p.mean ~ dunif(0,1)
  p.mu <- log(p.mean/(1-p.mean))
  p.b ~ dnorm(0, 0.37)
  
  for(i in 1:nSites){
    p.eps[i] ~ dnorm(0, p.tau)  # site-level variation
  }
  p.tau <- pow(p.sigma,-2)
  p.sigma ~ dunif(0,3)
  p.sigma2 <- pow(p.sigma, 2)
}  