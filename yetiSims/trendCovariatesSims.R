#This runs a set of simulations based on a control structure that uses simNum to choose setttings

## working directory & libraries
#using home directory instead of setting explicitly ("~/")
library(jagsUI); library(dplyr); library(data.table)

#########################
## Running Simulations ##
#########################   

#grab the simulation number from the environment passed from the slurm call
simNum<- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

simControl<-readRDS("simsToDo.rds") %>%
  .[simNum==whichSim]

cat("\n",simNum,"\n")
print(simControl)

## Simulation settings
nPasses <- 3


# some MCMC settings
n.chains = 3  		# number of chains
n.adapt = 5000		# number of sampler tuning iterations
n.burnin = 60000 # number of iterations to discard
n.iter = 65000	  # total iterations

#n.chains = 3  		# number of chains
#n.adapt = 5		# number of sampler tuning iterations
#n.burnin = 10 # number of iterations to discard
#n.iter = 20	  # total iterations



thin = 1				  # number to thin by
results<-NULL#create empty results to bind real results into
for(s in 1:nrow(simControl)){
  
  #include environmental covariates in simulations and models?
  includeCovs<-simControl$covariates[s]
  
  # model name
  model = ifelse(includeCovs,"trendModelCov.r","trendModelNoCov.r")
  
  #load data specific to stage and covariate inclusion
  load(paste0(simControl$stage[s],
              ifelse(includeCovs,"Cov","NoCov"),
              ".RData"))
  
  #get simulation settings from control structure
  nYears<-simControl$nYears[s]
  nSites<-simControl$nSites[s]
  r<-exp(simControl$trend[s])-1
  
  #set random seed based on simNum because all iterations were returning identical results
  set.seed(simNum*s)


      #give an update of progress
      cat("starting model ",s," of ",nrow(simControl))
      quantsToSave<-c(0.025,0.05,0.075,0.1,0.125,0.5,0.875,0.9,0.925,0.95,0.975)


# parameters to save
      if(includeCovs){
        #with environmental covariates
		pars.to.save <- c("mu","trend","sd.site","sd.year","sigma",
                 "p.mean","p.b","b1","b2","b3","b4","b5","b6")
	} else {
	  #without environmental covariates
		pars.to.save <- c("mu","trend","sd.site","sd.year","sigma",
                 "p.mean","p.b")
	}
      #set up results structures
      resultCols<-c('parameter','Mean',paste0("q",quantsToSave*100),'SD','rHat',
                    'trueValue')
      
      res = array(NA,dim=c(length(pars.to.save),length(resultCols))) %>%
        data.frame()
      names(res)<-resultCols
      res$parameter<-pars.to.save
    
      ## Data generation  
      #generate environment
      if(includeCovs){
        fallPrcp<-rnorm(nYears)
        fallTmean<-rnorm(nYears)
        winterPrcp<-rnorm(nYears)
        winterTmean<-rnorm(nYears)
        springPrcp<-rnorm(nYears)
        springTmean<-rnorm(nYears)
      }
      
      #set up data structures
      N <- lambda <- p <- array(NA_real_, dim=c(nSites, nYears),
                                dimnames=list(paste("site",1:nSites), 
                                              paste("year",1:nYears)))
      
      y <- array(NA_integer_, dim=c(nSites, nYears, nPasses),
                 dimnames=list(paste("site",1:nSites),
                               paste("year",1:nYears),
                               paste("pass",1:nPasses)))      
      
      ## Parameters for population abundance
      ## Grab a set of params from MCMC samples
      out.new = out2$sims.list
      chainLength = length(out.new$mu)
      i = sample(1:chainLength, 1)
      
      mu = out.new$mu[i]             # overall mean abundance at a site on a log scale
      sd.site = out.new$sd.site[i]   # variation among sites
      sd.year = out.new$sd.year[i]   # variation among years 
      sigma = out.new$sigma[i]       # over-dispersion
      #r = -0.05         # annual rate of population decrease (e.g. 5% decrease) 
      trend = log(1+r)   # convert to log scale for linear model: Dauwalter et al. (2010)
      
      #change format of environmental covariate betas to facilitate output storage
      if(simControl$covariates[s]==TRUE){
        for(bee in 1:6){
          assign(paste0("b",bee),out.new$b[i,bee])
        }
        b = out.new$b[i,] #environmental covariate betas
      } else {
        for(bee in 1:6){
          assign(paste0("b",bee),NA)
        }
      }
      
      ## Parameters for detection
      p.mean = out.new$p.mean[i]    # mean detection prob
      p.mu = log(p.mean/(1-p.mean))   # convert to logit scale
      p.b = out.new$p.b[i]          # effect size of a det cov (day of year)
      #p.sigma = out.new$p.sigma[i]  # variation among sites
      
      ### save true values - grabbed using names in pars.to.save
      truePars <- sapply(pars.to.save,get)
      
      #--------------------------------------------------------------------
      ## JAGS set up
      nInit<-ifelse(simControl$stage[s]=="yoy",5000,1200)
      init.vals <- function() list(N=array(nInit, dim=c(nSites, nYears)))
      
      #--------------------------------------------------------------------
      ## run the simulations
      
      ## stochastic factors affecting abundance
      site.ran = rnorm(nSites, 0, sd.site)    # variation among sites
      
      if(includeCovs){
        year.ran = b[1]*fallPrcp + b[2]*fallTmean + #environmental effects
                   b[3]*winterPrcp + b[4]*winterTmean +
                   b[5]*springPrcp + b[6]*springTmean +
                   rnorm(nYears, 0, sd.year)    # random variation among years
      }else{
       year.ran = rnorm(nYears,0,sd.year) 
      }
      eps = array(rnorm(nSites*nYears, 0, sigma), dim=c(nSites,nYears))  # over-dispersion
      
      ## stochastic factors affecting detection
      sampday = array(rnorm(nSites*nYears, 0, 1), dim=c(nSites,nYears)) # standardized cov
      #p.eps = array(rnorm(nSites*nYears, 0, p.sigma), dim=c(nSites,nYears))  # over-d
      
      ## Simulated data in a site-by-year format     
      for(i in 1:nSites){
        for(j in 1:nYears){
          # abundance
          lambda[i,j] <- exp(mu + (trend)*(j-1) + 
                               site.ran[i] + year.ran[j] + eps[i,j])
          N[i,j] <- rpois(1,lambda[i,j])
          # observed count  
          p[i,j] <- plogis(p.mu + p.b*sampday[i,j]) #+ p.eps[i,j])
          y[i,j,1] <- rbinom(1, N[i,j], p[i,j])
          y[i,j,2] <- rbinom(1, N[i,j]-y[i,j,1], p[i,j])
          y[i,j,3] <- rbinom(1, N[i,j]-y[i,j,1]-y[i,j,2], p[i,j])
        }
      }
      
      #--------------------------------------------------------------------
      ## Bundle data
      if(includeCovs){
        jags.data <- list(nSites=nSites, nYears=nYears, sampday=sampday, y=y,
        fallPrcp=fallPrcp, fallTmean=fallTmean, winterPrcp=winterPrcp,
        winterTmean=winterTmean, springPrcp=springPrcp, springTmean=springTmean)
      } else{
          jags.data <- list(nSites=nSites, nYears=nYears, sampday=sampday, y=y)              
      }

      #--------------------------------------------------------------------
      
      ## Run using jags::jagsUI in parallel
      
      #trying to balance convergence with high enough N inits to not throw initialization errors
      # but an error could happen, so using try() to avoid nullifying the whole iteration
      jagsWorked<-FALSE
      try(expr={
        dm.mcmc=jags(data=jags.data,inits=init.vals,parameters.to.save=pars.to.save,model.file=model,
                     n.iter=n.iter,n.thin=thin,n.chains=n.chains,n.adapt=n.adapt,n.burnin=n.burnin,
                     parallel=T)
        jagsWorked<-TRUE
      })
      
      if(!jagsWorked){
        #if jags throws an error (could happen if initial N values aren't high enough, but should be very rare),
        #then just save NAs as the result to indicate which models failed
        res$f<-as.numeric(NA)
        res$nYears<-nYears
        res$nSites<-nSites
        res$stage<-simControl$stage[s]
        res$simNum<-paste0("2.",simNum,".",s-1)
        results<-rbind(results,res)
        next
      }
      
      #save the results summary and simulation info
      summ = summary(dm.mcmc$samples)$statistics
      quan = summary(dm.mcmc$samples, quantile=quantsToSave)$quantile
      res[,paste0("q",quantsToSave*100)]<-quan[dimnames(quan)[[1]]!="deviance",]
      res[,"rHat"]<-dm.mcmc$summary[1:length(pars.to.save),
                                    which(dimnames(dm.mcmc$summary)[[2]]=="Rhat")]
      res[,c("Mean","SD")]<-summ[dimnames(quan)[[1]]!="deviance",c("Mean","SD")]
      res[,"f"]<-dm.mcmc$summary[1:length(pars.to.save),"f"]
      res[,"trueValue"]<-truePars
      res$nYears<-nYears
      res$nSites<-nSites
      res$stage<-simControl$stage[s]
      res$covariates<-includeCovs
      res$simNum<-paste0("2.",simNum,".",s-1)
      
      #bind with previous simulations
      results<-rbind(results,res)

      #save to output folder
      saveRDS(results, file=paste0("~/output/sim",simNum,'.rds'))
      
    }
#   }
# }
#save to output folder
saveRDS(results, file=paste0("~/output/sim",simNum,'.rds'))
cat("simNum: ",simNum) 