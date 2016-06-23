model{
  beta~dnorm(0,0.001)
  sigma~dunif(0,10)
  tau<-1/pow(sigma,2)
  
  trend~dnorm(0,0.001)
  
  for(i in 1:N){
    yExp[i]<-trend*year[i]+beta*x[i]
    y[i]~dnorm(yExp[i],tau)
  }
}