library(jagsUI)

xTrend<-seq(3,5,0.01)
x<-xTrend+rnorm(length(xTrend),0,1)

y<-x*2+rnorm(length(x),1)

jagsData<-list(x=x,
               y=y,
               N=length(y),
               year=1:length(y))
varsToMonitor=c("beta","sigma","trend")

out <- jags(
  data=jagsData,
  inits=NULL,
  model = "trendTest/testModel.R",
  parameters.to.save = varsToMonitor,
  n.chains=3,
  n.iter = 3000,
  n.thin = 1,
  n.burnin=2000,
  parallel=T)

