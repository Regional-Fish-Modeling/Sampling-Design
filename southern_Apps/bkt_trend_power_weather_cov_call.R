# 4/25/2016
# Power analysis for detecting a temporal trend in regional brook trout abundance
# Fit a model to southern Apps data compiled by Than Hitt
# Adult count is the response variable
# Put Poisson over-dispersion term back

# libraries
library(reshape2); library(rjags); library(plyr); library(ggplot2)
library(knitr); library(arm); library(boot); library(dplyr); library(jagsUI)
load.module("glm")

# Read in data from Than's data
## Trout count data
load("Data_FishCountAr.rData")
## Imports pass-specific count data for each site and year

## Detection covriates data
load("Data_DetectionCovsStd.RData")
## Imports sampling day-of-year and precip in prior 7 days from DAYMET data (standardized)

## Seasonal weather data 
load("Data_SeasonalClimateStd.RData")

#Site width for offset
#read data, make area, and grab columns of interest
siteLength<-read.csv("BKT_COVS_SAMPLES.csv",stringsAsFactors=F) %>%
  mutate(siteArea100m2=SiteWidth_m*SiteWidth_m/100) %>%
  mutate(siteLength100m=SiteLength_m/100) %>%
  dplyr::select(SiteID,siteArea100m2,siteLength100m,SiteWidth_m,Year) %>%
  data.frame()

siteLength[siteLength$Year==2105,"Year"]<-2015 #fix typo on year

siteLength<-siteLength %>%
  right_join(data.frame(SiteID=rep(unique(siteLength$SiteID),each=length(unique(siteLength$Year))),
                        Year=rep(unique(siteLength$Year),length(unique(siteLength$SiteID))),
                        stringsAsFactors=F),
             by=c("SiteID","Year")) %>%
  dplyr::group_by(SiteID) %>%
  dplyr::mutate(meanSiteLength100m=mean(siteLength100m,na.rm=T)) %>%
  ungroup()

#assign the mean by site for unsampled sites (and one sampled, but unmeasured site, with no interannual variability in site length)
siteLength[is.na(siteLength$siteLength100m),"siteLength100m"]<-
  siteLength[is.na(siteLength$siteLength100m),"meanSiteLength100m"]
siteLength<-dplyr::select(siteLength,-meanSiteLength100m)

#make an array, can change b/w area and length here by switching 3rd arugment of select
siteLength<-siteLength %>%
  dplyr::select(SiteID,Year,siteLength100m) %>%
  melt(id.vars=c("SiteID","Year")) %>%
  acast(SiteID~Year)
#subset to sites that are in ADUFish
siteLength<-siteLength[match(dimnames(ADUFish)[[1]],dimnames(siteLength)[[1]]),]


## Calculate "regional" mean by year
RegFallPrcpStd <- apply(FallPrcpStd, 2, mean)
RegFallTmeanStd <- apply(FallTmeanStd, 2, mean)
RegWinterPrcpStd <- apply(WinterPrcpStd, 2, mean)
RegWinterTmeanStd <- apply(WinterTmeanStd, 2, mean)
RegSpringPrcpStd <- apply(SpringPrcpStd, 2, mean)
RegSpringTmeanStd <- apply(SpringTmeanStd, 2, mean)


# prep for JAGS
## data structure
nSites = dim(ADUFish)[1]
nYears = dim(ADUFish)[2]
nCovs = 6 # number of seasonal weather covariates

## bundle data
dat <- list(nSites=nSites, nYears=nYears, y=ADUFish, nCovs=nCovs, 
            prcp7day=prcp7day.std, sampday=sampday.std,
            siteLength=siteLength,
            fallPrcp=RegFallPrcpStd, fallTmean=RegFallTmeanStd,
            winterPrcp=RegWinterPrcpStd, winterTmean=RegWinterTmeanStd,
            springPrcp=RegSpringPrcpStd, springTmean=RegSpringTmeanStd)

# set initial values
init <- function() list(mu=runif(1,0,5), 
                        N=array(1000, dim=c(nSites, nYears)),
                        sd.site=runif(1,0,2), 
                        sd.year=runif(1,0,2), sigma=runif(1,0,2),
                        p.mean=runif(1,0.2,0.8), p.b=rnorm(1), 
                        b=rnorm(nCovs),
                        sigma.b=array(runif(nCovs,0,2)))

# parameters to monitor
pars <- c("mu","trend","sd.site","sd.year","b","sigma","p.mean","p.b")


## running with rjags
#burnin <- jags.model(paste("bkt_trend_power_weather_cov_model.r", sep=""),
#                           dat, init, n.chains=3, n.adapt=30000)

# mcmc sample
#out <- coda.samples(burnin, pars, n.iter=30000, thin=10)
#summary(out)
#plot(out)

# gelman r value

#library(coda) 
#gelman.diag(out, multivariate=FALSE)

## parallel run in jagsUI
out2 <- jags(dat, init, pars, paste("bkt_trend_power_weather_cov_model.r", sep=""),
             n.chains=3, n.thin=10, n.iter=60000, n.burnin=30000, parallel=TRUE)
print(out2, dig=3)
