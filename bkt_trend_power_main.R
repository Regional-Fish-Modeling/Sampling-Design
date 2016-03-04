## working directory & libraries

# setwd("D:/Kanno/bkt_sample_design/sim1")
outDir<-"simTest"
getwd()
library(reshape2); library(jagsUI); library(plyr); library(ggplot2)
library(knitr); library(arm); library(boot)

## read in model output based on analysis of adult fish in southern range
load("bkt trend power out.rdata")


#########################
## Running Simulations ##
#########################   

## Simulation settings
startTime<-Sys.time()
nSims <- 1
nPasses <- 3

# parameters to save
pars.to.save <- c("mu","trend","sd.site","sd.slope","sd.year","sigma",
                  "p.mean","p.b","p.sigma")

# model name
model = "bkt trend power model.r"

# some MCMC settings
n.chains = 3  		# number of chains
n.adapt = 1000		# number of sampler tuning iterations
n.burnin = 1 # number of iterations to discard
n.iter = 25000	  # total iterations
thin = 1				  # number to thin by

#loop through settings
for(nSites in c(50,100,200)){   # number of sites
  for(nYears in c(5,10,20)){    # number of years
    for(r in c(-0.01,-0.025,-0.05)){   # percent annual decline (1%, 2.5%, 5%)
      
      # create array to capture results
      res = array(NA, dim=c(nSims, length(pars.to.save), 9), 
                  dimnames=list(1:nSims, pars.to.save, 
                                c('Mean','2.5%','7.5%','12.5%','Median','87.5%','92.5%','97.5%','SD')))
      
      # create list to save posteriors
      models = vector('list', nSims)
      
      # create array to check convergence  
      gelmanR = array(NA, dim=c(nSims, length(pars.to.save)),
                      dimnames=list(1:nSims, pars.to.save))
      
      #create array to save true parameter values
      truePars = array(NA,dim=c(nSims,length(pars.to.save)),
                       dimnames=list(1:nSims,pars.to.save))
      
      #loop through simulations
      for (s in 1:nSims){   # number of simulations
        ## Data generation
        N <- lambda <- p <- array(NA, dim=c(nSites, nYears),
                                  dimnames=list(paste("site",1:nSites), 
                                                paste("year",1:nYears)))
        
        y <- array(NA, dim=c(nSites, nYears, nPasses),
                   dimnames=list(paste("site",1:nSites),
                                 paste("year",1:nYears),
                                 paste("pass",1:nPasses)))      
        
        ## Parameters for population abundance
        ## Grab a set of params from MCMC samples
        out.new = as.matrix(out)
        chainLength = dim(out.new)[1]
        i = sample(1:nrow(out.new), 1)
        
        mu = out.new[i,"mu"]             # overall mean abundance at a site on a log scale
        sd.site = out.new[i,"sd.site"]   # variation among sites
        sd.year = out.new[i,"sd.year"]   # variation among years 
        sd.slope = out.new[i,"sd.slope"] # variation in trend among sites
        sigma = out.new[i,"sigma"]        # over-dispersion
        #r = -0.05         # annual rate of population decrease (e.g. 5% decrease) 
        trend = log(1+r)   # convert to log scale for linear model: Dauwalter et al. (2010)
        
        ## Parameters for detection
        p.mean = out.new[i,"p.mean"]    # mean detection prob
        p.mu = log(p.mean/(1-p.mean))   # convert to logit scale
        p.b = out.new[i,"p.b"]          # effect size of a det cov (day of year)
        p.sigma = out.new[i,"p.sigma"]  # variation among sites
        
        ### save true values - grabbed using names in pars.to.save
        truePars[s,] <- sapply(pars.to.save,get)
        
        #--------------------------------------------------------------------
        ## JAGS set up
        
        init.vals <- function() list(mu=runif(1,0,10), 
                                     N=array(1000, dim=c(nSites, nYears)),
                                     sd.slope=runif(1,0,5), sd.site=runif(1,0,5),
                                     sd.year=runif(1,0,5), sigma=runif(1,0,5),
                                     p.mean=runif(1,0,1), p.b=rnorm(1))
        
        #--------------------------------------------------------------------
        ## run the simulations
        set.seed(123)
        
        ## stochastic factors affecting abundance
        site.ran = rnorm(nSites, 0, sd.site)    # variation among sites
        year.ran = rnorm(nYears, 0, sd.year)    # variation among years
        slope.ran = rnorm(nSites, 0, sd.slope)  # variation in trend among sites
        eps = array(rnorm(nSites*nYears, 0, sigma), dim=c(nSites,nYears))  # over-dispersion
        
        ## stochastic factors affecting detection
        sampday = array(rnorm(nSites*nYears, 0, 1), dim=c(nSites,nYears)) # standardized cov
        p.eps = array(rnorm(nSites*nYears, 0, p.sigma), dim=c(nSites,nYears))  # over-d
        
        ## Simulated data in a site-by-year format     
        for(i in 1:nSites){
          for(j in 1:nYears){
            # abundance
            lambda[i,j] <- exp(mu + (trend + slope.ran[i])*(j-1) + 
                                 site.ran[i] + year.ran[j] + eps[i,j])
            N[i,j] <- rpois(1,lambda[i,j])
            # observed count  
            p[i,j] <- plogis(p.mu + p.b*sampday[i,j] + p.eps[i,j])
            y[i,j,1] <- rbinom(1, N[i,j], p[i,j])
            y[i,j,2] <- rbinom(1, N[i,j]-y[i,j,1], p[i,j])
            y[i,j,3] <- rbinom(1, N[i,j]-y[i,j,1]-y[i,j,2], p[i,j])
          }
        }
        
        #--------------------------------------------------------------------
        ## Bundle data
        jags.data <- list(nSites=nSites, nYears=nYears, sampday=sampday, y=y)
        
        cat(paste('nSites=',nSites,'nYears=',nYears,'r=',r,'Simulation',s," started ",Sys.time(),'\n'))
        
        #--------------------------------------------------------------------
        ## Run using jags::jagsUI in parallel
        dm.mcmc=jags(data=jags.data,inits=init.vals,parameters.to.save=pars.to.save,model.file=model,
                     n.iter=n.iter,n.thin=thin,n.chains=n.chains,n.adapt=n.adapt,n.burnin=n.burnin,
                     parallel=T)
        
          ## save posterior samples
        models[[s]][[1]] = jags.data
        models[[s]][[2]] = dm.mcmc$samples
        
        summ = summary(dm.mcmc$samples)$statistics
        quan = summary(dm.mcmc$samples, quantile=c(0.025,0.075,0.125,0.5,0.875,0.925,0.975))$quantile
        
        res[s,,1] = summ[1:length(pars.to.save),1]
        res[s,,2] = quan[1:length(pars.to.save),'2.5%']
        res[s,,3] = quan[1:length(pars.to.save),'7.5%']
        res[s,,4] = quan[1:length(pars.to.save),'12.5%']
        res[s,,5] = quan[1:length(pars.to.save),'50%']
        res[s,,6] = quan[1:length(pars.to.save),'87.5%']
        res[s,,7] = quan[1:length(pars.to.save),'92.5%']
        res[s,,8] = quan[1:length(pars.to.save),'97.5%']
        res[s,,9] = summ[1:length(pars.to.save),2]
        
        gelmanR[s,] = dm.mcmc$summary[1:length(pars.to.save),
                                      which(dimnames(dm.mcmc$summary)[[2]]=="Rhat")]
      }
        saveRDS(res, file=paste0(outDir,'/res','nSites',nSites,'nYears',nYears,'r',r,'.rds'))
        saveRDS(models, file=paste0(outDir,'/models','nSites',nSites,'nYears',nYears,'r',r,'_MCMC.rds'))
        saveRDS(gelmanR, file=paste0(outDir,'/gelmanR','nSites',nSites,'nYears',nYears,'r',r,'.rds'))
        saveRDS(truePars, file=paste0(outDir,'/truePars','nSites',nSites,'nYears',nYears,'r',r,'.rds'))
    }
  }
}
endTime<-Sys.time()
cat("simulations took ",round(difftime(endTime,startTime,units="h"),2)," hours")
